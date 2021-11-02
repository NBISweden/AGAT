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
    "help|h"      => \$opt_help,
    "gff1=s"      => \$gff1,
    "gff2=s"      => \$gff2,
    "v!"          => \$verbose,
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
print ("Parsing $gff1\n");
my ($omniscient1, $hash_mRNAGeneLink1) = slurp_gff3_file_JD({ input => $gff1
                                                              });
print ("\n\nParsing $gff2\n");                                                              
my ($omniscient2, $hash_mRNAGeneLink2) = slurp_gff3_file_JD({ input => $gff2
                                                              });
print ("-- Files parsed --\n");


my $sortBySeq1 = gather_and_sort_l1_location_by_seq_id_and_strand_chimere($omniscient1);
my $sortBySeq2 = gather_and_sort_l1_location_by_seq_id_and_strand_chimere($omniscient2);
print ("GFF3 files sorted\n");

#get top feature first
my $top_features = get_feature_type_by_agat_value($omniscient1, 'level1', 'topfeature');

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
# Will merge same location types that overlap
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
foreach my $locusID ( sort  keys %{$flattened_locations1} ){
  foreach my $chimere_type ( sort keys %{$flattened_locations1->{$locusID}} ){
    foreach my $level ( sort keys %{$flattened_locations1->{$locusID}{$chimere_type}} ){
      foreach my $type ( sort keys %{$flattened_locations1->{$locusID}{$chimere_type}{$level}} ){

        print "\n========================================================\nGENERAL loop over $locusID $chimere_type $level <<$type>>\n"if ($verbose);
        if ( exists_keys ($flattened_locations1, ($locusID,$chimere_type,$level,$type) ) ){ # We have to remove the locations2 to check at the end the FP that are remaining (only prenent in annotationB)

          if ($verbose) { print "list of location1 $level $type: "; foreach my $array ( @{$flattened_locations1->{$locusID}{$chimere_type}{$level}{$type}}){print "@{$array} - "; } print "\n";}
          while ( my $location1 = shift  @{$flattened_locations1->{$locusID}{$chimere_type}{$level}{$type}} ){ # here the location are supposed to be sorted
            print "location1 investigated:  @$location1\n" if ($verbose);

            # keep track last locationA
            my $last_locationA = undef;
            $last_locationA = 1 if (scalar @{$flattened_locations1->{$locusID}{$chimere_type}{$level}{$type}} == 0);
            print "Last LocationA !!\n" if ( $last_locationA and $verbose);

            if ( exists_keys ($flattened_locations2, ($locusID,$chimere_type,$level,$type) ) and
                scalar @{$flattened_locations2->{$locusID}{$chimere_type}{$level}{$type}} != 0 ){ # and

              my $previous_overlap = 0;
              while ( scalar @{$flattened_locations2->{$locusID}{$chimere_type}{$level}{$type}} != 0 ){
                if ($verbose) { print " list of location2 $level $type: "; foreach my $array ( @{$flattened_locations2->{$locusID}{$chimere_type}{$level}{$type}}){print "@{$array} - "; } print "\n";}
                my $shift_it = 1;
                my $location2 = $flattened_locations2->{$locusID}{$chimere_type}{$level}{$type}->[0];
                print " location2 investigated:  @$location2\n" if ($verbose);
                print " Original TP: ".$all{$chimere_type}{$level}{$type}{'TP'}."\n" if $verbose;
                print " Original FN: ".$all{$chimere_type}{$level}{$type}{'FN'}."\n" if $verbose;
                print " Original FP: ".$all{$chimere_type}{$level}{$type}{'FP'}."\n" if $verbose;

                # keep track last locationA
                my $last_locationB = undef;
                $last_locationB = 1 if (scalar @{$flattened_locations2->{$locusID}{$chimere_type}{$level}{$type}} == 1);
                print " Last LocationB !!\n" if ($last_locationB and $verbose);

                #  location A                         ----------------
                #  location B  ---------------
                if ($location2->[1] < $location1->[0]){
                  print " shift location2 because before location A!\n" if ($verbose);
                  my $FP = $location2->[1] - $location2->[0] + 1; #size
                  $all{$chimere_type}{$level}{$type}{'FP'} += $FP;
                  print " FP => $FP\n" if ($verbose);
                  $previous_overlap = 0;
                }

                # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
                #      ------------ OVERLAP -----------
                ## # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
                elsif( ($location1->[0] <= $location2->[1]) and ($location1->[1] >= $location2->[0])){
                  $previous_overlap = 1;
                  print " location @$location1 and @$location2 overlap !!!!\n" if $verbose;
                  my ($FN, $FP, $TP) = get_snsp_for_overlaps ($location1, $location2);
                  print " FN=$FN, FP=$FP, TP=$TP\n" if $verbose;

                  # check                  vvvv
                  #  location A    -------
                  #  location B  --------------
                  # but only if not the last location
                  my $locationB_remain = 0 ;
                  my $locationA_remain = 0 ;
                  my $clean_cut = 0 ;
                  if( $location2->[1] > $location1->[1] ){
                    $locationB_remain = $location2->[1] - $location1->[1] + 1;
                    $shift_it = undef;
                  }
                  # check                       vvvv
                  #  location A          -----------
                  #  location B  --------------
                  elsif( $location1->[1] > $location2->[1] ){
                    $locationA_remain = 1;
                  }
                  #                            v
                  #  location A        -------- <
                  #  location B  -------------- <
                  else{
                    $clean_cut = 1;
                  }

                  #  location A          -------------
                  #  location B  -----------   --  ------------
                  print "locationA_remain $locationA_remain \n" if $verbose;
                  print "locationB_remain $locationB_remain \n" if $verbose;
                  print "clean_cut $clean_cut \n" if $verbose;
                  if ($locationA_remain and !$last_locationA){
                    # TP must always be added
                    $all{$chimere_type}{$level}{$type}{'TP'} += $TP;
                    $all{$chimere_type}{$level}{$type}{'FN'} -= $TP;
                    $all{$chimere_type}{$level}{$type}{'FP'} += $FP;
                    print " FN: removing ".$TP."\n" if $verbose;
                    print " FP: ADDING ".$FP."\n" if $verbose;
                    print " TP: ADDING ".$TP."\n" if $verbose;
                  }
                  elsif ($locationB_remain and !$last_locationB){
                    # TP must always be added
                    $all{$chimere_type}{$level}{$type}{'TP'} += $TP;
                    $all{$chimere_type}{$level}{$type}{'FN'} += $FN;
                    $all{$chimere_type}{$level}{$type}{'FP'} -= $TP;
                    print " FN: removing ".$TP."\n" if $verbose;
                    print " FP: ADDING ".$FP."\n" if $verbose;
                    print " TP: ADDING ".$TP."\n" if $verbose;
                  }
                  else{
                    $all{$chimere_type}{$level}{$type}{'FN'} += $FN;
                    $all{$chimere_type}{$level}{$type}{'FP'} += $FP;
                    $all{$chimere_type}{$level}{$type}{'TP'} += $TP;
                    print " FN: ADDING ".$FN."\n" if $verbose;
                    print " FP: ADDING ".$FP."\n" if $verbose;
                    print " TP: ADDING ".$TP."\n" if $verbose;
                  }

                  # From previous locationA analysis
                  if($locationB_remain){
                    # No more locationA
                    if( $last_locationA ){
                      print " No more location A 1\n" if $verbose;
                      # No more locationB
                      if( $last_locationB ){
                        print " No more location B\n" if $verbose;
                      }
                      print " Next location B\n" if $verbose;
                      $shift_it = 1;
                    }
                    else{
                      print " Next location A\n" if $verbose;
                      last; # next locationA
                    }
                  }
                  elsif( $clean_cut ){
                    print " Clean cut !!! Removing LocationB and next Location A\n" if $verbose;
                    my $tothrow = shift  @{$flattened_locations2->{$locusID}{$chimere_type}{$level}{$type}};# Throw location B
                    print " End3 TP: ".$all{$chimere_type}{$level}{$type}{'TP'}."\n" if $verbose;
                    print " End3 FN: ".$all{$chimere_type}{$level}{$type}{'FN'}."\n" if $verbose;
                    print " End3 FP: ".$all{$chimere_type}{$level}{$type}{'FP'}."\n\n" if $verbose;
                    last; # next locationA
                  }
                }

                #  location A  -------------------------
                #  location B                                     -------------------------
                else{
                  print " last because location2 after\n" if ($verbose);

                  if($previous_overlap){
                    my $FP = $location2->[1] - $location2->[0] + 1; #size
                    $all{$chimere_type}{$level}{$type}{'FP'} += $FP;
                    print " Take into account the current locationB! +FP FP;\n" if ($verbose);
                  }
                  else{
                    my $FN = $location1->[1] - $location1->[0] + 1; #size
                    $all{$chimere_type}{$level}{$type}{'FN'} += $FN;
                    print " Take into account the current locationA! +FN:$FN;\n" if ($verbose);
                    $shift_it = undef;
                  }

                  # If it was overlaping then we do not count it because it has been already taken into account
                  print "End2 TP: ".$all{$chimere_type}{$level}{$type}{'TP'}."\n" if $verbose;
                  print "End2 FN: ".$all{$chimere_type}{$level}{$type}{'FN'}."\n" if $verbose;
                  print "End2 FP: ".$all{$chimere_type}{$level}{$type}{'FP'}."\n\n" if $verbose;


                  # no more location1
                  if(scalar @{$flattened_locations1->{$locusID}{$chimere_type}{$level}{$type}} == 0 ){
                    my $tothrow = shift  @{$flattened_locations2->{$locusID}{$chimere_type}{$level}{$type}};# Throw location B
                    print "No more location A 2\n"  if $verbose;
                    next;
                  }
                  else{
                    $previous_overlap = 0;
                    last; # next locationA
                  }
                }

                #
                # - END -
                #
                if ($shift_it){
                  my $tothrow = shift  @{$flattened_locations2->{$locusID}{$chimere_type}{$level}{$type}};# Throw location B
                  print " Removing location2: @$tothrow\n" if ($verbose);
                }
                print "End1 TP: ".$all{$chimere_type}{$level}{$type}{'TP'}."\n" if $verbose;
                print "End1 FN: ".$all{$chimere_type}{$level}{$type}{'FN'}."\n" if $verbose;
                print "End1 FP: ".$all{$chimere_type}{$level}{$type}{'FP'}."\n\n" if $verbose;
              }# END WHILE until location B is after A
            }

            # The list of locationB is empty now
            else{
              my $FN += $location1->[1] - $location1->[0] + 1; #size
              print " LocationA only => +FN:$FN\n" if ($verbose);
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

# ---- NOw deal with what is remaining in annotationB => FP
# Gather False positive => seq only annotated in annotationB, or type of feature annotated only in annotationB that was missing in annotatoin A
foreach my $locusID (  keys %{$flattened_locations2} ){
  foreach my $chimere_type ( keys %{$flattened_locations2->{$locusID}}){
    foreach my $level ( keys %{$flattened_locations2->{$locusID}{$chimere_type}} ){
      foreach my $type ( keys %{$flattened_locations2->{$locusID}{$chimere_type}{$level}} ){

        if ( exists_keys ($flattened_locations2, ($locusID,$chimere_type,$level,$type) ) ){ # We have to remove the locations2 to check at the end the FP that are remaining (only prenent in annotationB)
          while ( my $location2 = shift @{$flattened_locations2->{$locusID}{$chimere_type}{$level}{$type}} ){ # here the location are supposed to be sorted
            my $FP = $location2->[1] - $location2->[0] + 1; #size
            print "remaining $chimere_type $level $type - location: ".$location2->[0]." ".$location2->[1]."  -  +FP $FP\n" if ($verbose);
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
        print "chimere_type:$chimere_type level:$level type/$type TP:$TP FN:$FN FP:$FP\n"  if $verbose;
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
# should it be implemented?

# ------------------------------------------------------------------------------
# ------------------------- Now print the Results -------------------------
# ------------------------------------------------------------------------------
#if ($verbose) {use Data::Dumper; print "The sensitivity hash: ".Dumper(\%sensitivity)."\nThe specificity hash: ".Dumper(\%specificity);}
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

if ($outfile){
  print $report $string_to_print;
}
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
      $FN = $location1->[1] - $location2->[1]; #size
    }
    #  location A  --------------
    #  location B  -------------------------
    #( $location1->[1] < $location2->[1] )
    else{
      $TP = $location1->[1] - $location1->[0] + 1; #size
      $FP = $location2->[1] - $location1->[1]; #size
    }
  }
  # ---- SAME END ---
  elsif( $location1->[1] == $location2->[1]){
    #  location A  -------------------------
    #  location B             --------------
    if( $location1->[0] < $location2->[0] ){
      $TP = $location2->[1] - $location2->[0] + 1; #size
      $FN = $location2->[0] - $location1->[0]; #size
    }
    #  location A             --------------
    #  location B  -------------------------
    #( $location1->[0] > $location2->[0] )
    else{
      $TP = $location1->[1] - $location1->[0] + 1; #size
      $FP = $location1->[0] - $location2->[0]; #size
    }
  }
  #  ---- LOCATION A START BEFORE --- Not same start/end
  elsif( $location1->[0] < $location2->[0]){
    #  location A  -------------------------
    #  location B          -----------
    if($location1->[1] > $location2->[1]){
      $TP = $location2->[1] - $location2->[0] + 1; #size
      $FN = $location2->[0] - $location1->[0]; #size
      $FN += $location1->[1] - $location2->[1]; #size
    }
    #  location A  -------------------------
    #  location B             ---------------------
    # ( $location1->[1] < $location2->[1] )
    else{
      $TP = $location1->[1] - $location2->[0] + 1; #size
      $FP = $location2->[1] - $location1->[1]; #size
      $FN = $location2->[0] - $location1->[0]; #size
    }
  }
  #  ---- LOCATION A START AFTER ---  Not same start/end
  # ( $location1->[0] > $location2->[0] )
  else{
    #  location A       -----------
    #  location B  -------------------------
    if($location1->[1] < $location2->[1]){
      $TP = $location1->[1] - $location1->[0] + 1; #size
      $FP = $location1->[0] - $location2->[0]; #size
      $FP += $location2->[1] - $location1->[1]; #size
    }
    #  location A         -------------------------
    #  location B  ---------------------
    # ( $location1->[1] > $location2->[1] )
    else{
      $TP = $location2->[1] - $location1->[0] + 1; #size
      $FP = $location1->[0] - $location2->[0]; #size
      $FN = $location1->[1] - $location2->[1]; #size
    }
  }
  return  $FN, $FP, $TP;
}

__END__

=head1 NAME

agat_sp_sensitivity_specificity.pl

=head1 DESCRIPTION

The script aims to compute the Sensitivity and Specificity in order to assess the quality
of an annotation according to a reference (that is supposed to be true high-quality annotation).
The Sensitivity (Sn) is the proportion of true predictions compared to the total number of correct genes (including missed predictions)
Sn = TP / TP+FN
The Specificity (Sp) is the proportion of true predictions among all predicted genes (including incorrectly predicted ones)
Sp = TP / TP+FP

reference annotation:     -------------
prediction          :           ------------
                            FN     TP    FP    TN

Sensitivity and Specificity will be computed for each feature types.
(and computed independentaly if part of different Level2 type. i.e. exons Sn Sp
for tRNA will not be mixed up with the exon Sn Sp of mRNA exons)

=head1 SYNOPSIS

    agat_sp_sensitivity_specificity.pl --gff1 infile1.gff --gff2 infile2.gff  [ -o outfile ]
    agat_sp_sensitivity_specificity.pl --help

=head1 OPTIONS

=over 8

=item B<-gff1>

Input GTF/GFF file 1.

=item B<-gff2>

Input GTF/GFF file 2.

=item B<-o> , B<--output> , B<--out> or B<--outfile>

Output GFF file.  If no output file is specified, the output will be
written to STDOUT.

=item B<-v>

Verbose option for debug purposes.

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
