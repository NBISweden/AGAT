#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long;
use Pod::Usage;
use List::MoreUtils qw(uniq);
use Sort::Naturally;
use AGAT::AGAT;

start_script();
my $header = get_agat_header();
# ------------------------------- LOAD OPTIONS --------------------------------
my $outfile = undef;
my $gff1 = undef;
my $gff2 = undef;
my $opt_help= 0;

# OPTION MANAGEMENT: split shared vs script options
my ($shared_argv, $script_argv) = split_argv_shared_vs_script(\@ARGV);

# Parse script-specific options
my $script_parser = Getopt::Long::Parser->new;
$script_parser->configure('bundling','no_auto_abbrev');
if ( ! $script_parser->getoptionsfromarray(
  $script_argv,
  'h|help!'                     => \$opt_help,
  'gff1=s'                      => \$gff1,
  'gff2=s'                      => \$gff2,
  'output|outfile|out|o=s'      => \$outfile,
  ) )
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

# --- Manage config ---
my $shared_opts = parse_shared_options($shared_argv);
initialize_agat({ config_file_in => ( $shared_opts->{config} ), input => $gff1, shared_opts => $shared_opts });

######################
# Manage output file #
my $report = prepare_fileout($outfile);

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>     MAIN     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

######################
### Parse GFF input #
my ($omniscient1) = slurp_gff3_file_JD({ input => $gff1 });

my $log = create_log_file({input => $gff2});
$LOGGING->{'log'} = $log ;
my ($omniscient2) = slurp_gff3_file_JD({ input => $gff2 });

my $sortBySeq1 = gather_and_sort_l1_location_by_seq_id_and_strand_chimere($omniscient1);
my $sortBySeq2 = gather_and_sort_l1_location_by_seq_id_and_strand_chimere($omniscient2);
dual_print1 "GFF3 files sorted\n";

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
dual_print2 "Now flattening the locations\n";
foreach my $flattened_locations ( $flattened_locations1, $flattened_locations2 ){
  foreach my $locusID (  keys %{$flattened_locations} ){
    foreach my $chimere_type ( keys %{$flattened_locations->{$locusID}}){
      foreach my $level ( keys %{$flattened_locations->{$locusID}{$chimere_type}} ){
        foreach my $type ( keys %{$flattened_locations->{$locusID}{$chimere_type}{$level}} ){

          # Initialise all counter to 0. Useful later to compute the FN FP TP
          $all{$chimere_type}{$level}{$type}{'FN'}=0;
          $all{$chimere_type}{$level}{$type}{'FP'}=0;
          $all{$chimere_type}{$level}{$type}{'TP'}=0;

          dual_print2 "investigate $type\n";
          my @newlocations;
          my $previous_location = undef;
          foreach my $location ( sort {$a->[0] <=> $b->[0]} @{$flattened_locations->{$locusID}{$chimere_type}{$level}{$type}} ){
            dual_print2 "investigate @$location\n";
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
dual_print2 "COMPARE FLATENED LOCATIONS\n";
foreach my $locusID ( sort  keys %{$flattened_locations1} ){
  foreach my $chimere_type ( sort keys %{$flattened_locations1->{$locusID}} ){
    foreach my $level ( sort keys %{$flattened_locations1->{$locusID}{$chimere_type}} ){
      foreach my $type ( sort keys %{$flattened_locations1->{$locusID}{$chimere_type}{$level}} ){

        dual_print2 "\n========================================================\nGENERAL loop over $locusID $chimere_type $level <<$type>>\n";
        if ( exists_keys ($flattened_locations1, ($locusID,$chimere_type,$level,$type) ) ){ # We have to remove the locations2 to check at the end the FP that are remaining (only prenent in annotationB)

          if ($CONFIG->{verbose}) { dual_print2 "list of location1 $level $type: "; foreach my $array ( @{$flattened_locations1->{$locusID}{$chimere_type}{$level}{$type}}){ dual_print2 "@{$array} - "; } dual_print2 "\n";}
          while ( my $location1 = shift  @{$flattened_locations1->{$locusID}{$chimere_type}{$level}{$type}} ){ # here the location are supposed to be sorted
            dual_print2 "location1 investigated:  @$location1\n";

            # keep track last locationA
            my $last_locationA = undef;
            $last_locationA = 1 if (scalar @{$flattened_locations1->{$locusID}{$chimere_type}{$level}{$type}} == 0);
            if ($last_locationA){ dual_print2 "Lets go for last LocationA !!\n"; }

            if ( exists_keys ($flattened_locations2, ($locusID,$chimere_type,$level,$type) ) and
                scalar @{$flattened_locations2->{$locusID}{$chimere_type}{$level}{$type}} != 0 ){ # and

              while ( scalar @{$flattened_locations2->{$locusID}{$chimere_type}{$level}{$type}} != 0 ){
                if ($CONFIG->{verbose} >1) { dual_print2 " list of location2 $level $type: "; foreach my $array ( @{$flattened_locations2->{$locusID}{$chimere_type}{$level}{$type}}){ dual_print2 "@{$array} - "; } dual_print2 "\n";}

                my $location2 = $flattened_locations2->{$locusID}{$chimere_type}{$level}{$type}->[0];
                dual_print2 " location2 investigated:  @$location2\n";
                dual_print2 " Original TP: ".$all{$chimere_type}{$level}{$type}{'TP'}."\n";
                dual_print2 " Original FN: ".$all{$chimere_type}{$level}{$type}{'FN'}."\n";
                dual_print2 " Original FP: ".$all{$chimere_type}{$level}{$type}{'FP'}."\n";

                # keep track last locationB
                my $last_locationB = undef;
                $last_locationB = 1 if (scalar @{$flattened_locations2->{$locusID}{$chimere_type}{$level}{$type}} == 1);
                if ($last_locationB){ dual_print2 " Lets go for last LocationB !!\n"; }

                # ===================== CASE 1 =====================
                #  location A                         ----------------
                #  location B  ---------------
                if ($location2->[1] < $location1->[0]){
                  my $FP = $location2->[1] - $location2->[0] + 1; #size
                  dual_print2 " +FP => $FP\n";
                  $all{$chimere_type}{$level}{$type}{'FP'} += $FP;

                  dual_print2 "End1 TP: ".$all{$chimere_type}{$level}{$type}{'TP'}."\n";
                  dual_print2 "End1 FN: ".$all{$chimere_type}{$level}{$type}{'FN'}."\n";
                  dual_print2 "End1 FP: ".$all{$chimere_type}{$level}{$type}{'FP'}."\n";
                  dual_print2 " Case1 - Next location2!\n\n";
                  my $tothrow = shift  @{$flattened_locations2->{$locusID}{$chimere_type}{$level}{$type}};# Throw location B
                }


                # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
                #      ------------ OVERLAP -----------
                ## # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
                elsif( ($location1->[0] <= $location2->[1]) and ($location1->[1] >= $location2->[0])){
                  dual_print2 " location @$location1 and @$location2 overlap !!!!\n";
                  my ($FN, $FP, $TP) = get_snsp_for_overlaps ($location1, $location2);
                  dual_print2 " FN=$FN, FP=$FP, TP=$TP\n";


                  my $locationB_remain = 0 ;
                  my $locationA_remain = 0 ;
                  # check                  vvvv
                  #  location A    -------
                  #  location B  --------------
                  # but only if not the last location
                  if( $location2->[1] > $location1->[1] ){
                    $locationB_remain = $location2->[1] - $location1->[1] + 1;
                  }
                  # check                       vvvv
                  #  location A          -----------
                  #  location B  --------------
                  elsif( $location1->[1] > $location2->[1] ){
                    $locationA_remain = 1;
                  }
                  #  location A          -------------
                  #  location B  -----------   --  ------------

                  dual_print2 "locationA_remain $locationA_remain \n";
                  dual_print2 "locationB_remain $locationB_remain \n";

                  if ($locationA_remain and !$last_locationA and !$last_locationB){
                    # TP must always be added
                    $all{$chimere_type}{$level}{$type}{'TP'} += $TP;
                    $all{$chimere_type}{$level}{$type}{'FN'} -= $TP;
                    $all{$chimere_type}{$level}{$type}{'FP'} += $FP;
                    dual_print2 " TP: ADDING ".$TP."\n";
                    dual_print2 " FN: removing ".$TP."\n";
                    dual_print2 " FP: ADDING ".$FP."\n";
                    dual_print2 " Case2 A - Next location B \n";
                    my $tothrow = shift  @{$flattened_locations2->{$locusID}{$chimere_type}{$level}{$type}};# Throw location B
                  }

                  elsif ($locationB_remain and !$last_locationA  and !$last_locationB){
                    # TP must always be added
                    $all{$chimere_type}{$level}{$type}{'TP'} += $TP;
                    $all{$chimere_type}{$level}{$type}{'FN'} += $FN;
                    $all{$chimere_type}{$level}{$type}{'FP'} -= $TP;

                    dual_print2 " TP: ADDING ".$TP."\n";
                    dual_print2 " FN: ADDING ".$FN."\n";
                    dual_print2 " FP: removing ".$TP."\n";
                    dual_print2 " Case2 B - Next location A\n";
                    last;
                  }

                  else{ #clean cut or end of one type of location (1 or 2)
                    $all{$chimere_type}{$level}{$type}{'FN'} += $FN;
                    $all{$chimere_type}{$level}{$type}{'FP'} += $FP;
                    $all{$chimere_type}{$level}{$type}{'TP'} += $TP;
                    dual_print2 " TP: ADDING ".$TP."\n";
                    dual_print2 " FN: ADDING ".$FN."\n";
                    dual_print2 " FP: ADDING ".$FP."\n";
                    dual_print2 " Case2 C ------\n";
                    dual_print2 " End TP: ".$all{$chimere_type}{$level}{$type}{'TP'}."\n";
                    dual_print2 " End FN: ".$all{$chimere_type}{$level}{$type}{'FN'}."\n";
                    dual_print2 " End FP: ".$all{$chimere_type}{$level}{$type}{'FP'}."\n\n";

                      if( $last_locationB and !$last_locationA){
                        dual_print2 " Case2 C2 - Remove last location B\n";
                        my $tothrow = shift  @{$flattened_locations2->{$locusID}{$chimere_type}{$level}{$type}};# Throw location B
                      }
                      elsif ($last_locationA and !$last_locationB){
                        dual_print2 " Case2 C3 - No more location A - LAST\n";
                        last;
                      }
                      elsif ($last_locationA and $last_locationB){
                        dual_print2 " Case2 C4 - No more locationA neither locationB. Removing locationB and LAST.\n";
                        my $tothrow = shift  @{$flattened_locations2->{$locusID}{$chimere_type}{$level}{$type}};# Throw location B
                        last;
                      }
                      #  clean cut                 v
                      #  location A        -------- <
                      #  location B  -------------- <
                      # No more locationA
                      elsif(!$locationA_remain and !$locationB_remain){
                        dual_print2 " Clean cut !!! Removing LocationB and next Location A\n";
                        my $tothrow = shift  @{$flattened_locations2->{$locusID}{$chimere_type}{$level}{$type}};# Throw location B
                        last; # next locationA
                      }
                  }
                }

                # ===================== CASE 3 =====================
                #  location A  -------------------------
                #  location B                                     -------------------------
                else{
                  dual_print2 " last because location2 after\n";

                  my $FN = $location1->[1] - $location1->[0] + 1; #size
                  $all{$chimere_type}{$level}{$type}{'FN'} += $FN;
                  dual_print2 " Take into account the current locationA! +FN: $FN;\n";

                  dual_print2 " End2 TP: ".$all{$chimere_type}{$level}{$type}{'TP'}."\n";
                  dual_print2 " End2 FN: ".$all{$chimere_type}{$level}{$type}{'FN'}."\n";
                  dual_print2 " End2 FP: ".$all{$chimere_type}{$level}{$type}{'FP'}."\n";
                  dual_print2 " Case3 - Next location A \n\n";
                  last; # next locationA
                }
              }# END WHILE until location B is after A
            }

            # ===================== CASE 4 =====================
            # The list of locationB is empty now
            else{
              my $FN += $location1->[1] - $location1->[0] + 1; #size
              dual_print2 " LocationA only => +FN:$FN\n";
              $all{$chimere_type}{$level}{$type}{'FN'} += $FN;
              dual_print2 " Case4 - Next location A \n\n";
            }
          }
        }

        # No such $type in annotationB, so it is specific to annotationA
        else{
          dual_print2 "Specific to annotationA  => FN\n";
          my $FN=0;
          foreach my $location ( @{$flattened_locations1->{$locusID}{$chimere_type}{$level}{$type}} ){ # here the location are supposed to be sorted
            $FN += $location->[1] - $location->[0] + 1; #size
            dual_print2 " Case5 - Next location A \n";
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
            dual_print2 "remaining $chimere_type $level $type - location: ".$location2->[0]." ".$location2->[1]."  -  +FP $FP\n";
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
        dual_print2 "chimere_type:$chimere_type level:$level type/$type TP:$TP FN:$FN FP:$FP\n";
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
#if ($CONFIG->{verbose} >1) {use Data::Dumper; dual_print2 "The sensitivity hash: ".Dumper(\%sensitivity)."\nThe specificity hash: ".Dumper(\%specificity);} 
my $string_to_print =  join('', '-') x 64;
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
dual_print1 $string_to_print;

# --- final messages ---
end_script();
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

=item B<-h> or B<--help>

Display this helpful text.

=back

=head1 SHARED OPTIONS

Shared options are defined in the AGAT configuration file and can be overridden via the command line for this script only.
Common shared options are listed below; for the full list, please refer to the AGAT agat_config.yaml.

=over 8

=item B<--config>

String - Path to a custom AGAT configuration file.  
By default, AGAT uses `agat_config.yaml` from the working directory if present, otherwise the default file shipped with AGAT
(available locally via `agat config --expose`).

=item B<--cpu>, B<--core>, B<--job> or B<--thread>

Integer - Number of parallel processes to use for file input parsing (via forking).

=item B<-v> or B<--verbose>

Integer - Verbosity, choice are 0,1,2,3,4. 0 is quiet, 1 is normal, 2,3,4 is more verbose. Default 1.

=back

=head1 FEEDBACK

For questions, suggestions, or general discussions about AGAT, please use the AGAT community forum:
https://github.com/NBISweden/AGAT/discussions

=head1 BUG REPORTING

Bug reports should be submitted through the AGAT GitHub issue tracker:
https://github.com/NBISweden/AGAT/issues

=cut

AUTHOR - Jacques Dainat
