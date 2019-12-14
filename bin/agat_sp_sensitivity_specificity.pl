#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long;
use Pod::Usage;
use List::MoreUtils qw(uniq);
use Sort::Naturally;
use Bio::Tools::GFF;
use AGAT::Omniscient;

my $header = get_agat_header();
my $outfile = undef;
my $gff1 = undef;
my $gff2 = undef;
my $verbose = undef;
my $opt_help= 0;

my @copyARGV=@ARGV;
if ( !GetOptions(
    "help|h" => \$opt_help,
    "gff1=s" => \$gff1,
    "gff2=s" => \$gff2,
    "v!" => \$verbose,
    "output|outfile|out|o=s" => \$outfile))

{
    pod2usage( { -message => 'Failed to parse command line',
                 -verbose => 1,
                 -exitval => 1 } );
}

# Print Help and exit
if ($opt_help) {
    pod2usage( { -verbose => 99,
                 -exitval => 0,
                 -message => "$header\n" } );
}

if ( ! $gff1 or ! $gff2){
    pod2usage( {
           -message => "$header\nAt least 2 parameters are mandatory:\n".
                       "Input reference gff file1 (--gff1)\n".
                       "Input reference gff file2 (--gff2)\n\n",
           -verbose => 0,
           -exitval => 1 } );
}

######################
# Manage output file #
my $report = IO::File->new();
if ($outfile) {
  open($report, '>', $outfile) or die "Could not open file $outfile $!";
}
else{
  $report->fdopen( fileno(STDOUT), 'w' );
}

                #####################
                #     MAIN          #
                #####################


######################
### Parse GFF input #
my ($omniscient1, $hash_mRNAGeneLink1) = slurp_gff3_file_JD({ input => $gff1
                                                              });
my ($omniscient2, $hash_mRNAGeneLink2) = slurp_gff3_file_JD({ input => $gff2
                                                              });
print ("GFF3 files parsed\n");


my $sortBySeq1 = gather_and_sort_l1_location_by_seq_id_and_strand_chimere($omniscient1);
my $sortBySeq2 = gather_and_sort_l1_location_by_seq_id_and_strand_chimere($omniscient2);
print ("GFF3 files sorted\n");

#get top feature first
my $hash = get_levels_info(); # get from the file
my $top_features = $hash->{'other'}{'level'}{'topfeature'};

# ----- Remove $top_features ------
foreach my $sortBySeq ($sortBySeq1, $sortBySeq2){
  foreach my $locusID1 ( sort keys %{$sortBySeq}){ # tag_l1 = gene or repeat etc...

    # Skip to features
    foreach my $type_top_feature (keys %{$top_features}){
      if (exists_keys( $sortBySeq, ($locusID1, $type_top_feature) ) ){
        delete  $sortBySeq->{$locusID1}{$type_top_feature};
      }
    }
  }
}
# ----- END Remove $top_features ------

my %all; # <= Will contain all TP FN FP values

# ------------------------------------------------------------------------------
# ------------------------- STORE ALL LOCATIONS -------------------------
# ------------------------------------------------------------------------------
my $flattened_locations1 = {};
my $flattened_locations2 = {};
my $cpt = 0;
#use Data::Dumper; print "\n\n\n sortBySeq1 location1: ".Dumper($sortBySeq1) ;
#use Data::Dumper; print "\n\n\n sortBySeq2 : ".Dumper($sortBySeq2) ;
foreach my $sortBySeq ($sortBySeq1, $sortBySeq2){

  # select proper variable to work with (related to annotationA or annotationB)
  my $flattened_locations;
  my $omniscient;
  if(! $cpt){
    $cpt++;
    $flattened_locations = $flattened_locations1;
    $omniscient = $omniscient1;
  }
  else{
    $flattened_locations = $flattened_locations2;
    $omniscient = $omniscient2;
  }
  foreach my $locusID ( sort keys %{$sortBySeq}){ # tag_l1 = gene or repeat etc...
    foreach my $chimere_type_l1l2 ( sort keys %{$sortBySeq->{$locusID}} ) {
      #
      # Go through location from left to right ### !! if not empty
      #
      my $previous_location_l1 = undef;
      my $list_of_location_l2 = [];
      while ( my $location1 = shift @{$sortBySeq->{$locusID}{$chimere_type_l1l2}} ){

        # Define location l1
        my $current_location_l1 = [$location1->[1], $location1->[2]];
        my $l1_id = lc($location1->[0]);
        my $type_l1 = $location1->[3];

        push @{$flattened_locations->{$locusID}{$chimere_type_l1l2}{'level1'}{$type_l1}}, $current_location_l1 ;

        ################
        # Go to level2 #
        ################
        foreach my $type_l2 (keys %{$omniscient->{'level2'}}){
          if(exists_keys($omniscient,('level2', $type_l2, $l1_id) ) ){

            foreach my $l2_f ( @{$omniscient->{'level2'}{$type_l2}{$l1_id} } ){

              # Define location l2
              my $l2_id = lc($l2_f->_tag_value('ID'));
              my $current_location_l2 = [int($l2_f->start()), int($l2_f->end())];

              push @{$flattened_locations->{$locusID}{$chimere_type_l1l2}{'level2'}{$type_l2}}, $current_location_l2 ;

              ################
              # Go to level3 #
              ################
              foreach my $type_l3 (keys %{$omniscient->{'level3'}}){
                if(exists_keys($omniscient,('level3', $type_l3, $l2_id) ) ){
                  foreach my $l3_f ( @{$omniscient->{'level3'}{$type_l3}{$l2_id} } ){
                    # Define location l3
                    my $current_location_l3 = [int($l3_f->start()), int($l3_f->end()) ];
                    push @{$flattened_locations->{$locusID}{$chimere_type_l1l2}{'level3'}{$type_l3}}, $current_location_l3 ;
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}


# ------------------------------------------------------------------------------
# --------------------------- FIX OVERLAPPING LOCATIONS ------------- ----------
# ------------------------------------------------------------------------------
#use Data::Dumper; print "\n\n\n all locations1 sored: ".Dumper($flattened_locations1) ;
#use Data::Dumper; print "\n\n\n all locations2 sored: ".Dumper($flattened_locations2) ;
print "Now flattening the locations\n" if ($verbose);
foreach my $flattened_locations ( $flattened_locations1, $flattened_locations2 ){
  foreach my $locusID (  keys %{$flattened_locations} ){
    foreach my $chimere_type ( keys %{$flattened_locations->{$locusID}}){
      foreach my $level ( keys %{$flattened_locations->{$locusID}{$chimere_type}} ){
        foreach my $type ( keys %{$flattened_locations->{$locusID}{$chimere_type}{$level}} ){

          # Initialise all counter to 0. Useful later to compute the FN FP TP
          $all{$chimere_type}{$level}{$type}{'FN'}=0;
          $all{$chimere_type}{$level}{$type}{'FP'}=0;
          $all{$chimere_type}{$level}{$type}{'TP'}=0;

          print "investigate $type\n" if ($verbose);
          my @newlocations;
          my $previous_location = undef;
          foreach my $location ( sort {$a->[0] <=> $b->[0]} @{$flattened_locations->{$locusID}{$chimere_type}{$level}{$type}} ){
            print "investigate @$location\n" if ($verbose);
            # first round
            if (! $previous_location){
               push @newlocations, $location;
               $previous_location = $location;
            }
            # Not first round
            else{
              #  location A  -------------------------
              #  location B           -------------------------
              if ( ($previous_location->[0] <= $location->[1]) and ($previous_location->[1] >= $location->[0])){

                #  location A  -------------------------
                #  location B           -------------------------
                if($previous_location->[1] <= $location->[1]){
                  # take back last location pushed in the array
                  $previous_location = pop @newlocations;
                  $previous_location  = [$previous_location->[0], $location->[1]];
                  # push back into the array the previous location that has been modified
                  push @newlocations, $previous_location ;
                }
              }
              else{
                push @newlocations, $location ;
                $previous_location = $location;
              }
            }
          }
          @{$flattened_locations->{$locusID}{$chimere_type}{$level}{$type}} = @newlocations ;
        }
      }
    }
  }
}


# ------------------------------------------------------------------------------
# ---- NOW COMPARE FLATENED LOCATIONS OF THE TWO ANNOTAITON sorted by chimere name (l1l2)
# ------------------------------------------------------------------------------
#use Data::Dumper; print "\n\n\n flattened_locations1: ".Dumper($flattened_locations1) ;
#use Data::Dumper; print "\n\n\n flattened_locations2: ".Dumper($flattened_locations2) ;
print "COMPARE FLATENED LOCATIONS\n" if ($verbose);
foreach my $locusID (  keys %{$flattened_locations1} ){
  foreach my $chimere_type ( keys %{$flattened_locations1->{$locusID}} ){
    foreach my $level ( keys %{$flattened_locations1->{$locusID}{$chimere_type}} ){
      foreach my $type ( keys %{$flattened_locations1->{$locusID}{$chimere_type}{$level}} ){

        print "\nGENERAL loop over $level $type\n"if ($verbose);
        if ( exists_keys ($flattened_locations1, ($locusID,$chimere_type,$level,$type) ) ){ # We have to remove the locations2 to check at the end the FP that are remaining (only prenent in annotationB)

          my $previous_overlap = 0;
          my $previous_FP_right_2 = 0 ;
          if ($verbose) { print "list of location1 $level $type: "; foreach my $array ( @{$flattened_locations1->{$locusID}{$chimere_type}{$level}{$type}}){print "@{$array} - "; } print "\n";}
          while ( my $location1 = shift  @{$flattened_locations1->{$locusID}{$chimere_type}{$level}{$type}} ){ # here the location are supposed to be sorted
            print "location1 investigated:  @$location1\n" if ($verbose);


            if ( exists_keys ($flattened_locations2, ($locusID,$chimere_type,$level,$type) ) and
                scalar @{$flattened_locations2->{$locusID}{$chimere_type}{$level}{$type}} != 0 ){ # and


              while ( scalar @{$flattened_locations2->{$locusID}{$chimere_type}{$level}{$type}} != 0 ){
                if ($verbose) { print " list of location2 $level $type: "; foreach my $array ( @{$flattened_locations2->{$locusID}{$chimere_type}{$level}{$type}}){print "@{$array} - "; } print "\n";}
                my $shift_it = 1;
                my $location2 = $flattened_locations2->{$locusID}{$chimere_type}{$level}{$type}->[0];
                print " location2 investigated:  @$location2\n" if ($verbose);

                #  location A                         ----------------
                #  location B  ---------------
                if ($location2->[1] < $location1->[0]){
                  print " shift location2 because before location A!\n" if ($verbose);
                  if ($previous_overlap){
                    print " the previous location overlapped!\n" if ($verbose);
                    # If last locationB here we should take into accout locationA otherwise is lost
                    if (! @{$flattened_locations2->{$locusID}{$chimere_type}{$level}{$type}}){
                      my $FN = $location1->[1] - $location1->[0] + 1; #size
                      $all{$chimere_type}{$level}{$type}{'FN'} += $FN;
                      print " last locationB must take into account the current location A! FN => $FN\n" if ($verbose);
                    }
                  }
                  else{
                    my $FP = $location2->[1] - $location2->[0] + 1; #size
                    $all{$chimere_type}{$level}{$type}{'FP'} += $FP;
                    print " FP => $FP\n" if ($verbose);
                  }
                  $previous_overlap = 0;
                  $previous_FP_right_2 = 0;
                }

                # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
                #      ------------ OVERLAP -----------
                ## # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
                elsif( ($location1->[0] <= $location2->[1]) and ($location1->[1] >= $location2->[0])){
                  print " location @$location1 and @$location2 overlap !!!!\n" if $verbose;
                  my ($FN, $FP, $TP) = get_snsp_for_overlaps ($location1, $location2);
                  print " FN=$FN, FP=$FP, TP=$TP\n" if $verbose;

                  # get FP right           vvvv
                  #  location A    -------
                  #  location B  --------------
                  if($FP){
                    if ($location2->[1] > $location1->[1] ){
                      $previous_FP_right_2 = $location2->[1] - $location1->[1] + 1;
                      $shift_it = undef;
                    }
                  }

                  if (! $previous_overlap){
                    $all{$chimere_type}{$level}{$type}{'FN'} += $FN;
                    $all{$chimere_type}{$level}{$type}{'FP'} += $FP;
                    $all{$chimere_type}{$level}{$type}{'TP'} += $TP;

                  }
                  #  location A          -------------
                  #  location B  -----------   --  ------------
                  else{
                    # TP must always be added
                    $all{$chimere_type}{$level}{$type}{'TP'} += $TP;
                    $all{$chimere_type}{$level}{$type}{'FN'} -= $TP;
                    $all{$chimere_type}{$level}{$type}{'FP'} += $FP;
                  }

                  # From previous locationA analysis
                  if($previous_FP_right_2){
                    $all{$chimere_type}{$level}{$type}{'FP'} -= $TP;
                  }

                  # At the end of this foreach we should remove the locationB if it is last element
                  # and nothing is left in locationA list
                  # if yes we must remove it to not count it as FP
                  #if ( ! @{$flattened_locations1->{$locusID}{$chimere_type}{$level}{$type}} and  scalar @{$flattened_locations2->{$locusID}{$chimere_type}{$level}{$type}} == 1){
                  #    print "remove $level $type from locations2\n" if ($verbose);
                  #}
                  $previous_overlap = 1;
                }

                #  location A  -------------------------
                #  location B                                     -------------------------
                else{
                  print " last because location2 after\n" if ($verbose);
                  if ( ! $previous_overlap){
                    # We will take another location A, the current one is not yet taken into account if was not overlaping
                    my $FN = $location1->[1] - $location1->[0] + 1; #size
                    $all{$chimere_type}{$level}{$type}{'FN'} += $FN;
                    print " Take into account the current locationA! FN = $FN;\n" if ($verbose);
                    $shift_it = undef;
                  }
                  else{$previous_overlap = 0;}

                  last; # Go back to the list of locationA
                }

                if ($shift_it){
                  my $tothrow = shift  @{$flattened_locations2->{$locusID}{$chimere_type}{$level}{$type}};# Throw location B
                  print " Remove location2: @$tothrow\n" if ($verbose);
                }
              }# END WHILE until location B is after A
            }


            # The list of locationB is empty now
            else{
              print " LocationA only => FN\n" if ($verbose);
              my $FN += $location1->[1] - $location1->[0] + 1; #size
              $all{$chimere_type}{$level}{$type}{'FN'} += $FN;
            }
          }
        }

        # No such $type in annotationB, so it is specific to annotationA
        else{
          print "Specific to annotationA  => FN\n" if ($verbose);
          my $FN=0;
          foreach my $location ( @{$flattened_locations1->{$locusID}{$chimere_type}{$level}{$type}} ){ # here the location are supposed to be sorted
            $FN += $location->[1] - $location->[0] + 1; #size
          }
          $all{$chimere_type}{$level}{$type}{'FN'} += $FN;
        }
      }
    }
  }
}
#print "is already out??\n";exit;
if ($verbose) { use Data::Dumper; print "\n\n\nFlatenned location1: ".Dumper($flattened_locations1) ;}
if ($verbose) { use Data::Dumper; print "\n\n\nFlatenned location2: ".Dumper($flattened_locations2) ;}
if ($verbose) { use Data::Dumper; print "The all hash: ".Dumper(\%all); }
# ---- NOw deal with what is remaining in annotationB => FP
# Gather False positive => seq only annotated in annotationB, or type of feature annotated only in annotationB that was missing in annotatoin A
foreach my $locusID (  keys %{$flattened_locations2} ){
  foreach my $chimere_type ( keys %{$flattened_locations2->{$locusID}}){
    foreach my $level ( keys %{$flattened_locations2->{$locusID}{$chimere_type}} ){
      foreach my $type ( keys %{$flattened_locations2->{$locusID}{$chimere_type}{$level}} ){

        if ( exists_keys ($flattened_locations2, ($locusID,$chimere_type,$level,$type) ) ){ # We have to remove the locations2 to check at the end the FP that are remaining (only prenent in annotationB)
          while ( my $location2 = shift @{$flattened_locations2->{$locusID}{$chimere_type}{$level}{$type}} ){ # here the location are supposed to be sorted
            my $FP = $location2->[1] - $location2->[0] + 1; #size
            $all{$chimere_type}{$level}{$type}{'FP'} += $FP;
          }
        }
      }
    }
  }
}

# ------------------------------------------------------------------------------
# ------------------------- Now compute Sn Sp -------------------------
# ------------------------------------------------------------------------------
my %sensitivity;
my %specificity;
foreach my $chimere_type ( keys %all ){
  foreach my $level ( keys %{$all{$chimere_type}} ){
    foreach my $type ( keys %{$all{$chimere_type}{$level}} ){

      if ( exists_keys (\%all, ($chimere_type,$level,$type) ) ){ # We have to remove the locations2 to check at the end the FP that are remaining (only prenent in annotationB)
        my $FN=$all{$chimere_type}{$level}{$type}{'FN'};
        my $FP=$all{$chimere_type}{$level}{$type}{'FP'};
        my $TP=$all{$chimere_type}{$level}{$type}{'TP'};
        if($TP){
          $sensitivity{$chimere_type}{$level}{$type} = sprintf("%.2f", $TP / ($TP + $FN) );
          $specificity{$chimere_type}{$level}{$type} = sprintf("%.2f", $TP / ($TP + $FP) );
        }
        else{
          $sensitivity{$chimere_type}{$level}{$type} = 0;
          $specificity{$chimere_type}{$level}{$type} = 0;
        }
      }
    }
  }
}


# ------------------------------------------------------------------------------
# ------------------------- Now compute the opposite -------------------------
# ------------------------------------------------------------------------------
 if ($verbose) {use Data::Dumper; print "The sensitivity hash: ".Dumper(\%sensitivity)."\nThe specificity hash: ".Dumper(\%specificity);}

# ------------------------------------------------------------------------------
# ------------------------- Now print the Results -------------------------
# ------------------------------------------------------------------------------
my $string_to_print = "usage: $0 @copyARGV\nResults:\n\n";
$string_to_print .=  join('', '-') x 64;
$string_to_print .= "\n|".sizedPrint("Feature type",20)."|".sizedPrint("Sensitivity",20)."|".sizedPrint("Specificity",20)."|\n";
foreach my $chimere_type ( sort keys %all ){
  if ( exists_keys ( \%all, ( $chimere_type, 'level1') ) ){
    foreach my $type ( sort keys %{$all{$chimere_type}{'level1'}} ){
      $string_to_print .=  join('', '-') x 64;
      $string_to_print .= "\n|".sizedPrint($type,20)."|".sizedPrint($sensitivity{$chimere_type}{'level1'}{$type},20)."|".sizedPrint($specificity{$chimere_type}{'level1'}{$type}, 20)."|\n";
      if ( exists_keys ( \%all, ( $chimere_type, 'level2') ) ){
        foreach my $type ( sort keys %{$all{$chimere_type}{'level2'}} ){
          $string_to_print .= "|".sizedPrint($type,20)."|".sizedPrint($sensitivity{$chimere_type}{'level2'}{$type},20)."|".sizedPrint($specificity{$chimere_type}{'level2'}{$type}, 20)."|\n";
          if ( exists_keys ( \%all, ( $chimere_type, 'level3') ) ){
            foreach my $type ( sort keys %{$all{$chimere_type}{'level3'}} ){
              $string_to_print .= "|".sizedPrint($type,20)."|".sizedPrint( $sensitivity{$chimere_type}{'level3'}{$type},20)."|".sizedPrint( $specificity{$chimere_type}{'level3'}{$type}, 20)."|\n";
            }
          }
        }
      }
    }
  }
}
$string_to_print .=  join('', '-') x 64;
$string_to_print .= "\n";
print $string_to_print;
print "Bye Bye.\n";
#######################################################################################################################
        ####################
         #     METHODS    #
          ################
           ##############
            ############
             ##########
              ########
               ######
                ####
                 ##

 # @Purpose: Create a hash of level1 location (location = [level1ID,start,end]) sorted by feature type and localisation. A localisation is the sequence_id appended by the strand
 # @input: 1 => hash omniscient
 # @output: 1 => hash => LocusID->typeFeatureChimere =[ID,start,end, type]
 # TagChimere allows to divide the L1 into relatedted l2 type. (e.g like do split by level 2 feature)
 sub gather_and_sort_l1_location_by_seq_id_and_strand_chimere{
 	my ($omniscient) = @_;

 	my %hash_sortBySeq;

 	foreach my $tag_level1 (keys %{$omniscient->{'level1'}} ){
   	foreach my $level1_id (keys %{$omniscient->{'level1'}{$tag_level1}} ){
	    my $level1_feature = $omniscient->{'level1'}{$tag_level1}{$level1_id};

      my $tag_chimere = $tag_level1;
      my $l2_type_keep=undef;
      foreach my $l2_type (keys %{$omniscient->{'level2'}} ){
        if( exists_keys($omniscient,('level2', $l2_type, $level1_id ) ) ) {
          $l2_type_keep = $l2_type;
        }
      }
      if($l2_type_keep){$tag_chimere .= $l2_type_keep;}

    	my $ID = $level1_feature->_tag_value('ID');
	    my $strand="+";
	    if($level1_feature->strand != 1){$strand = "-";}
	    my $position_l1=$level1_feature->seq_id."".$strand;
	    push ( @{$hash_sortBySeq{$position_l1}{$tag_chimere}}, [$ID, int($level1_feature->start), int($level1_feature->end), $tag_level1] );
    }

    foreach my $position_l1 (keys %hash_sortBySeq){
      foreach my $tag_chimere (keys %{$hash_sortBySeq{$position_l1}} ){
        @{$hash_sortBySeq{$position_l1}{$tag_chimere}} = sort { ncmp ( $a->[1], $b->[1] ) } @{$hash_sortBySeq{$position_l1}{$tag_chimere}};
      }
    }
 	}
 	return \%hash_sortBySeq;
}

# location provided must overlap
sub get_snsp_for_overlaps{
  my ($location1, $location2)=@_;

  my $TP = 0;
  my $FN = 0;
  my $FP = 0;

  #  ---- SAME START ---
  if( $location1->[0] == $location2->[0]){
  #  location A  -------------------------
  #  location B  -------------------------
    if ($location1->[1] == $location2->[1]){
      $TP = $location1->[1] - $location1->[0] + 1; #size

    }
    #  location A  -------------------------
    #  location B  --------------
    elsif( $location1->[1] > $location2->[1] ){
      $TP = $location2->[1] - $location2->[0] + 1; #size
      $FN = $location1->[1] - $location2->[1] + 1; #size
    }
    #  location A  --------------
    #  location B  -------------------------
    #( $location1->[1] < $location2->[1] )
    else{
      $TP = $location1->[1] - $location1->[0] + 1; #size
      $FP = $location2->[1] - $location1->[1] + 1; #size
    }
  }
  # ---- SAME END ---
  elsif( $location1->[1] == $location2->[1]){
    #  location A  -------------------------
    #  location B             --------------
    if( $location1->[0] < $location2->[0] ){
      $TP = $location2->[1] - $location2->[0] + 1; #size
      $FN = $location2->[0] - $location1->[0] + 1; #size
    }
    #  location A             --------------
    #  location B  -------------------------
    #( $location1->[0] > $location2->[0] )
    else{
      $TP = $location1->[1] - $location1->[0] + 1; #size
      $FP = $location1->[0] - $location2->[0] + 1; #size
    }
  }
  #  ---- LOCATION A START BEFORE --- Not same start/end
  elsif( $location1->[0] < $location2->[0]){
    #  location A  -------------------------
    #  location B          -----------
    if($location1->[1] > $location2->[1]){
      $TP = $location2->[1] - $location2->[0] + 1; #size
      $FN = $location2->[0] - $location1->[0] + 1; #size
      $FN += $location1->[1] - $location2->[1] + 1; #size
    }
    #  location A  -------------------------
    #  location B             ---------------------
    # ( $location1->[1] < $location2->[1] )
    else{
      $TP = $location1->[1] - $location2->[0] + 1; #size
      $FP = $location2->[1] - $location1->[1] + 1; #size
      $FN = $location2->[0] - $location1->[0] + 1; #size
    }
  }
  #  ---- LOCATION A START AFTER ---  Not same start/end
  # ( $location1->[0] > $location2->[0] )
  else{
    #  location A       -----------
    #  location B  -------------------------
    if($location1->[1] < $location2->[1]){
      $TP = $location1->[1] - $location1->[0] + 1; #size
      $FP = $location1->[0] - $location2->[0] + 1; #size
      $FP += $location2->[1] - $location1->[1] + 1; #size
    }
    #  location A         -------------------------
    #  location B  ---------------------
    # ( $location1->[1] > $location2->[1] )
    else{
      $TP = $location2->[1] - $location1->[0] + 1; #size
      $FP = $location1->[0] - $location2->[0] + 1; #size
      $FN = $location1->[1] - $location2->[1] + 1; #size
    }
  }
  return  $FN, $FP, $TP;
}

__END__

=head1 NAME

agat_sp_filter_by_locus_distance.pl

=head1 DESCRIPTION

The script aims to remove or flag loci that are too close to each other.
Close loci are important to remove when training abinitio tools in order
to train intergenic region properly. Indeed if intergenic region
(surrouneded part of a locus) contain part of another locus,
the training on intergenic part will be biased.

=head1 SYNOPSIS

    agat_sp_filter_by_locus_distance.pl -gff infile.gff [ -o outfile ]
    agat_sp_filter_by_locus_distance.pl --help

=head1 OPTIONS

=over 8

=item B<-gff>

Input GTF/GFF file.

=item B<--dist> or B<-d>

The minimum inter-loci distance to allow.  No default (will not apply
filter by default).

=item B<--add> or B<--add_flag>

Instead of filter the result into two output files, write only one and add the flag <low_dist> in the gff.(tag = Lvalue or tag = Rvalue  where L is left and R right and the value is the distance with accordingle the left or right locus)

=item B<-o> , B<--output> , B<--out> or B<--outfile>

Output GFF file.  If no output file is specified, the output will be
written to STDOUT.

=item B<-v>

Verbose option, make it easier to follow what is going on for debugging purpose.

=item B<-h> or B<--help>

Display this helpful text.

=back

=head1 FEEDBACK

=head2 Did you find a bug?

Do not hesitate to report bugs to help us keep track of the bugs and their
resolution. Please use the GitHub issue tracking system available at this
address:

            https://github.com/NBISweden/AGAT/issues

 Ensure that the bug was not already reported by searching under Issues.
 If you're unable to find an (open) issue addressing the problem, open a new one.
 Try as much as possible to include in the issue when relevant:
 - a clear description,
 - as much relevant information as possible,
 - the command used,
 - a data sample,
 - an explanation of the expected behaviour that is not occurring.

=head2 Do you want to contribute?

You are very welcome, visit this address for the Contributing guidelines:
https://github.com/NBISweden/AGAT/blob/master/CONTRIBUTING.md

=cut

AUTHOR - Jacques Dainat
