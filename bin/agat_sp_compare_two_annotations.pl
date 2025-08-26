#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use File::Basename;
use List::MoreUtils qw(uniq);
use Sort::Naturally;
use AGAT::AGAT;

my $header = get_agat_header();
my @copyARGV = @ARGV;

my ( $opt, $usage, $config ) = AGAT::AGAT::describe_script_options( $header,
    [ 'gff1=s', 'Input reference gff file 1', { required => 1 } ],
    [ 'gff2=s', 'Input reference gff file 2', { required => 1 } ],
);

my $gff1       = $opt->gff1;
my $gff2       = $opt->gff2;
my $opt_output = $config->{output};
my $verbose    = $config->{verbose};
my $debug      = $config->{debug};

my $log;
if ( my $log_name = $config->{log_path} ) {
    open( $log, '>', $log_name )
      or die "Can not open $log_name for printing: $!";
    dual_print( $log, $header,  3 );
}

######################
# Manage output folder #

if ( !$opt_output ) {
  dual_print( $log, "Default output name: split_result\n");
  $opt_output = 'comparison_result';
}

if ( -d $opt_output ) {
  warn "The output directory choosen already exists. Please give me another Name.\n" if $verbose;
  exit();
}
mkdir $opt_output;

######################
# Manage output file #
my $report = prepare_fileout("$opt_output/report.txt");

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>     MAIN     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# Activate verbose when debug active
$verbose=1 if ($debug);

######################
### Parse GFF input #
dual_print( $log, "Parsing $gff1\n");
my ($omniscient1, $hash_mRNAGeneLink1) = slurp_gff3_file_JD({ input => $gff1,
                                                              config => $config
                                                              });
dual_print( $log, "\n\nParsing $gff2\n");
my ($omniscient2, $hash_mRNAGeneLink2) = slurp_gff3_file_JD({ input => $gff2,
                                                              config => $config
                                                              });
dual_print( $log, "-- Files parsed --\n");


my $sortBySeq1 = gather_and_sort_l1_location_by_seq_id_and_strand_chimere($omniscient1);
my $sortBySeq2 = gather_and_sort_l1_location_by_seq_id_and_strand_chimere($omniscient2);
dual_print( $log, "GFF3 files sorted\n");

my %overlap_info; # <= Will contain the overlap information
# ------------------------------------------------------------------------------
# ---------------------------- COLLECT ALL LOCATIONS ---------------------------
# ------------------------------------------------------------------------------
dual_print( $log, "COLLECT LOCATIONS\n");
my $bucket_locations1 = {};
my $bucket_locations2 = {};
my $cpt = 0;

foreach my $sortBySeq ($sortBySeq1, $sortBySeq2){

  # select proper variable to work with (related to annotationA or annotationB)
  my $bucket_locations;
  my $omniscient;
  if(! $cpt){
    $cpt++;
    $bucket_locations = $bucket_locations1;
    $omniscient = $omniscient1;
  }
  else{
    $bucket_locations = $bucket_locations2;
    $omniscient = $omniscient2;
  }
	# create bucket location. A location=[start, end, type, id]
  foreach my $locusID ( sort keys %{$sortBySeq}){ # tag_l1 = gene or repeat etc...
    foreach my $chimere_type_l1l2 ( sort keys %{$sortBySeq->{$locusID}} ) {
      #
      # Go through location from left to right ### !! if not empty
      #
      my $previous_location_l1 = undef;
      my $list_of_location_l2 = [];
      while ( my $location1 = shift @{$sortBySeq->{$locusID}{$chimere_type_l1l2}} ){

        # Define location l1
        my $l1_id = lc($location1->[0]);
        my $type_l1 = $location1->[3];
				my $current_location_l1 = [$location1->[1], $location1->[2], $type_l1, $l1_id];

        push @{$bucket_locations->{$locusID}{$chimere_type_l1l2}{$l1_id}{'level1'}}, $current_location_l1 ;

        ################
        # Go to level2 #
        ################
        foreach my $type_l2 (keys %{$omniscient->{'level2'}}){
          if(exists_keys($omniscient,('level2', $type_l2, $l1_id) ) ){

            foreach my $l2_f ( @{$omniscient->{'level2'}{$type_l2}{$l1_id} } ){

              # Define location l2
              my $l2_id = lc($l2_f->_tag_value('ID'));
              my $current_location_l2 = [int($l2_f->start()), int($l2_f->end()), $type_l1."@".$type_l2, $l1_id];

              push @{$bucket_locations->{$locusID}{$chimere_type_l1l2}{$l1_id}{'level2'}}, $current_location_l2 ;

              ################
              # Go to level3 #
              ################
              foreach my $type_l3 (keys %{$omniscient->{'level3'}}){
                if(exists_keys($omniscient,('level3', $type_l3, $l2_id) ) ){
                  foreach my $l3_f ( @{$omniscient->{'level3'}{$type_l3}{$l2_id} } ){
                    # Define location l3
                    my $current_location_l3 = [int($l3_f->start()), int($l3_f->end()), $type_l1."@".$type_l2."@".$type_l3, $l1_id ];
                    push @{$bucket_locations->{$locusID}{$chimere_type_l1l2}{$l1_id}{'level3'}{$type_l3}}, $current_location_l3 ;
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
# ------------------------- FLATTEN OVERLAPPING LOCATIONS ----------------------
# ------------------------------------------------------------------------------
# Will merge same location types that overlap
dual_print( $log, "FLATTEN LOCATIONS\n");
my $flattened_locations_clean;
my $flattened_locations1_clean = {};
my $flattened_locations2_clean = {};
my $cpt_fl;
foreach my $bucket_locations ( $bucket_locations1, $bucket_locations2 ){

  	if (! $cpt_fl){
		$flattened_locations_clean = $flattened_locations1_clean;
		$cpt_fl++;
  	}
  	else {
		$flattened_locations_clean = $flattened_locations2_clean;
  	}

  	foreach my $locusID (  keys %{$bucket_locations} ){
    	foreach my $chimere_type ( keys %{$bucket_locations->{$locusID}}){
      		foreach my $l1_id ( keys %{$bucket_locations->{$locusID}{$chimere_type}} ){
				# Try from l3 to l1 and take the first existing
				foreach my $level ( ("level3","level2","level1") ){
		  			if ( exists_keys ($bucket_locations, ($locusID,$chimere_type, $l1_id, $level) ) ){

						my @newlocations;
						# merge L3 locations from a same locus. CDS locations if any or exon if any or the rest if any, in this order
						if  ($level eq "level3"){

			  				# first get type of l3 to use
			  				my @types;
			  				if ( exists_keys ($bucket_locations, ($locusID,$chimere_type, $l1_id, $level, "cds") ) ){
			  					push @types, "cds";
			  				}
			  				elsif ( exists_keys ($bucket_locations, ($locusID,$chimere_type, $l1_id, $level, "exon") ) ){
								push @types, "exon";
			  				}
			  				else{
			    				foreach my $type ( keys %{$bucket_locations->{$locusID}{$chimere_type}{$l1_id}{$level}} ){
				  					push @types, $type;
			    				}
			  				}
			  
			 				foreach my $type (@types){
								my $previous_location = undef;

								foreach my $location ( sort {$a->[0] <=> $b->[0]} @{$bucket_locations->{$locusID}{$chimere_type}{$l1_id}{$level}{$type}} ){
                                                                    dual_print( $log, "investigate location of $l1_id from $chimere_type from $locusID at $level \n", $debug );
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

											if($previous_location->[1] <= $location->[1]){
												# take back last location pushed in the array
												$previous_location = pop @newlocations;
												$previous_location->[1]  = $location->[1];
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
							}
						}
						# merge L2 or l1 locations from a same locus
						else {
							my $previous_location = undef;
							foreach my $location ( sort {$a->[0] <=> $b->[0]} @{$bucket_locations->{$locusID}{$chimere_type}{$l1_id}{$level}} ){
                                                            dual_print( $log, "investigate location of $l1_id from $chimere_type from $locusID at $level \n", $debug );
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
										if($previous_location->[1] <= $location->[1]){
										# take back last location pushed in the array
										$previous_location = pop @newlocations;
										$previous_location->[1]  = $location->[1];
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
						}
						# Save flattened locations either from the deepest level => L3 L2 or L1 in this order
						push @{$flattened_locations_clean->{$locusID}{$chimere_type}}, [@newlocations] ; # each set of locations must be  kept as separate list
						last;
					}
	    		}
	  		}
		}
 	}
}

my $flattened_locations1_clean_sorted;
foreach my $locusID ( sort  keys %{$flattened_locations1_clean} ){
  foreach my $chimere_type ( sort keys %{$flattened_locations1_clean->{$locusID}} ){
	@{$flattened_locations1_clean_sorted->{$locusID}{$chimere_type}}	= sort {$a->[0][0] <=> $b->[0][0]} @{$flattened_locations1_clean->{$locusID}{$chimere_type}};
  }
}
my $flattened_locations2_clean_sorted;
foreach my $locusID ( sort  keys %{$flattened_locations2_clean} ){
	foreach my $chimere_type ( sort keys %{$flattened_locations2_clean->{$locusID}} ){
		@{$flattened_locations2_clean_sorted->{$locusID}{$chimere_type}}	= sort {$a->[0][0] <=> $b->[0][0]} @{$flattened_locations2_clean->{$locusID}{$chimere_type}};
	}
}
# ------------------------------------------------------------------------------
# ---- NOW COMPARE FLATENED LOCATIONS OF THE TWO ANNOTAITON sorted by chimere name (l1l2)
# ------------------------------------------------------------------------------

# When already checked the location is removed from flattened_locationsX_clean_sorted hash
dual_print( $log, "COMPARE LOCATIONS\n");
my %seen1;
my %seen2;
foreach my $locusID ( sort  keys %{$flattened_locations1_clean_sorted} ){
  foreach my $chimere_type ( sort keys %{$flattened_locations1_clean_sorted->{$locusID}} ){
		if ( exists_keys ($flattened_locations2_clean_sorted, ($locusID,$chimere_type) ) ){

			#/!\ MUST SORT LOCATIONS
			my @copy1 = @{$flattened_locations1_clean_sorted->{$locusID}{$chimere_type}};
                    dual_print( $log, scalar(@copy1)." locations 1\n", $debug );

			# With foreach the last location is properly handled automatically
			foreach my $locations1 (@copy1){
				if(! exists_keys(\%seen1, ($locations1->[0][3]) ) ) {

					if ( exists_keys ($flattened_locations2_clean_sorted, ($locusID,$chimere_type) ) ) {
						my @copy2 = @{$flattened_locations2_clean_sorted->{$locusID}{$chimere_type}};
                                            dual_print( $log, scalar(@copy2)." locations 2\n", $debug );

						# With foreach the last location is properly handled automatically
						foreach my $locations2 (@copy2){
							if(! exists_keys(\%seen2, ($locations2->[0][3] ) ) ){

               	# --------------------- OVERLAP --------------------------------
								if ( locations_overlap($locations1, $locations2) ){
                                                                    dual_print( $log, "overlap ".$locations1->[0][3]." ".$locations2->[0][3]."\n", $debug );

									remove_loc_by_id($flattened_locations1_clean_sorted, $locusID, $chimere_type, $locations1->[0][3]);
									$seen1{$locations1->[0][3]}++;
									remove_loc_by_id($flattened_locations2_clean_sorted, $locusID, $chimere_type, $locations2->[0][3]);
									$seen2{$locations2->[0][3]}++;

									my $flat_overlap_1 = $locations1;
									my $flat_overlap_2 = $locations2;
									my $overlap_A=1;
									my @overlap_A_id=($locations1->[0][3]);
									my $overlap_B=1;
									my @overlap_B_id=($locations2->[0][3]);
									my $loop="top";
									my $flip=2;
									my $current_locs;
									my $current_hash;
									my $current_seen;
									my $current_flat=();
									my $current_flat_oppo=();

									while ($flip){
                                                                            dual_print( $log, "  loop $loop\n", $debug );
										if ($loop eq "top"){
											$current_locs = $flattened_locations1_clean_sorted->{$locusID}{$chimere_type};
											$current_hash = $flattened_locations1_clean_sorted;
											$current_seen = \%seen1;
											$current_flat = $flat_overlap_1;
											$current_flat_oppo = $flat_overlap_2;
										} elsif ($loop eq "bot"){
											$current_locs = $flattened_locations2_clean_sorted->{$locusID}{$chimere_type};
											$current_hash = $flattened_locations2_clean_sorted;
											$current_seen = \%seen2;
											$current_flat = $flat_overlap_2;
											$current_flat_oppo = $flat_overlap_1;
										}

										my @overlap_loc;
										foreach my $locations (@{$current_locs}){
                                                                                    dual_print( $log, "  ".scalar(@{$current_locs})." locations $loop\n", $debug );

											#  location   -------------------------
											#  flat                                        -------------------------
											if ( $locations->[scalar(@{$locations})-1][1] < $current_flat_oppo->[0][0] ){
                                                                                            dual_print( $log, "  caseX\n", $debug );
												next;
											}
											#  location A                         ----------------
											#  location B  ---------------
											elsif ($current_flat_oppo->[scalar(@{$current_flat_oppo})-1][1] < $locations->[0][0] ){
                                                                                            dual_print( $log, "  caseY\n", $debug );
												last;
											}
											elsif ( locations_overlap($current_flat_oppo, $locations) ){
                                                                                            dual_print( $log, "  overlap in overlap ".$current_flat_oppo->[0][3]." ".$locations->[0][3]."\n", $debug );
												# keep track of ID that overlap to remove locations later out of the loop
												push @overlap_loc, $locations->[0][3];
												$current_flat = flatten_locations_and_merge($current_flat, $locations);

												if($loop eq "top"){
													$overlap_A++;
													push @overlap_A_id, $locations->[0][3];
												}
												elsif($loop eq "bot"){
													$overlap_B++;
													push @overlap_B_id, $locations->[0][3];
												}
											}
											else{
                                                                                            dual_print( $log, "  caseZ overlap in intron\n", $debug );
											}
										}

										#flip
										if ($loop and $loop eq "top"){
											$loop="bot";
										} elsif ($loop and $loop eq "bot") {
											$loop="top";
										}

										if(! @overlap_loc){
											$flip--;
										} else {
											foreach my $id (@overlap_loc){
												remove_loc_by_id($current_hash, $locusID, $chimere_type, $id);
												$current_seen->{$id}++;
											}
										}
									}
									push @{$overlap_info{$locations1->[0][2]}{$overlap_A}{$overlap_B}}, [ [@overlap_A_id], [@overlap_B_id] ];
								}
		# ----------------------------------- CASE 1 -----------------------------------
								#  location A                         ----------------
								#  location B  ---------------
								elsif ($locations2->[scalar(@{$locations2})-1][1] < $locations1->[0][0] ){

									push @{$overlap_info{$locations1->[0][2]}{0}{1}}, [[undef], [$locations2->[0][3]]];
                                                                    dual_print( $log, "Case1 notoverlap !\n\n", $debug );
									# throw loc2
									remove_loc_by_id($flattened_locations2_clean_sorted, $locusID, $chimere_type, $locations2->[0][3]);
									$seen2{$locations2->[0][3]}++;
									next;
								}

		# ----------------------------------- CASE 2 -----------------------------------
								#  location A  -------------------------
								#  location B                                     -------------------------
								elsif ($locations1->[scalar(@{$locations1})-1][1] < $locations2->[0][0] ){
									my $id1 = $locations1->[0][3];
                                                                    dual_print( $log, "Case2 notoverlap !\n", $debug );
									if(! exists_keys(\%seen1, ( $id1 ) ) ){ # else it has been dealed by overlap case
										push @{$overlap_info{$locations1->[0][2]}{1}{0}}, [ [$locations1->[0][3]], [undef] ];
										# throw loc1
										remove_loc_by_id($flattened_locations1_clean_sorted, $locusID, $chimere_type, $id1);
										$seen1{$id1}++;
									}
									last; # next location1
								}
							}
						}
					}
				}
			}
		}
  }
}

# ---- NOw deal with what is remaining in annotationA
# Gather False positive => seq only annotated in annotationB, or type of feature annotated only in annotationB that was missing in annotatoin A
dual_print( $log, "\nLook now what is specific to annotationA \n");
foreach my $locusID (  keys %{$flattened_locations1_clean_sorted} ){
  foreach my $chimere_type ( keys %{$flattened_locations1_clean_sorted->{$locusID}}){
    foreach my $locations1 ( @{$flattened_locations1_clean_sorted->{$locusID}{$chimere_type}} ){
			push @{$overlap_info{$locations1->[0][2]}{1}{0}}, [ [$locations1->[0][3]], [undef]];
                    dual_print( $log, " Case3 !\n", $debug );
    }
  }
}

# ---- NOw deal with what is remaining in annotationB-
# Gather False positive => seq only annotated in annotationB, or type of feature annotated only in annotationB that was missing in annotatoin A
dual_print( $log, "\nLook now what is specific to annotationB \n");
foreach my $locusID (  keys %{$flattened_locations2_clean_sorted} ){
  foreach my $chimere_type ( keys %{$flattened_locations2_clean_sorted->{$locusID}}){
    foreach my $locations2 ( @{$flattened_locations2_clean_sorted->{$locusID}{$chimere_type}} ){
			push @{$overlap_info{$locations2->[0][2]}{0}{1}},  [ [undef], [$locations2->[0][3]] ];
                    dual_print( $log, " Case4 !\n", $debug );
    }
  }
}

##############
# STATISTICS #
if($verbose){
  dual_print( $log, "Compute statistics for $gff1:\n");
	print_omniscient_statistics ({ input => $omniscient1 });

    dual_print( $log, "Compute statistics for $gff2:\n");
	print_omniscient_statistics ({ input => $omniscient2 });
}

# ------------------------------------------------------------------------------
# ------------------------- Now print the Results -------------------------
# ------------------------------------------------------------------------------
my %total;
my $separator_table = join('', '-') x 94;
$separator_table .= "\n";
my ($filename1,$path1,$ext1) = fileparse($gff1,qr/\.[^.]*/);
my ($filename2,$path2,$ext2) = fileparse($gff2,qr/\.[^.]*/);

my $string_to_print = "usage: $0 @copyARGV\nResults of number of genes from file1 that overlap genes from file2:\n\n";

my %file_handler;
foreach my $type_l1 ( sort keys %overlap_info ){

    $string_to_print .= "$separator_table|".sizedPrint("$type_l1",92)."|\n";
    $string_to_print .= "$separator_table|".sizedPrint($filename1.$ext1,30)."|".sizedPrint($filename2.$ext2,30)."|".sizedPrint("Number of cases",30)."|\n$separator_table";

	# Create file handlers
	
	foreach my $value1 ( sort {$a <=> $b} keys %{$overlap_info{$type_l1}} ){
      foreach my $value2 ( sort {$a <=> $b}  keys %{$overlap_info{$type_l1}{$value1}} ){
		# report ids file_handler
		my $report_ids = prepare_fileout("$opt_output/$type_l1"."_"."$value1"."_".$value2."_id_list.txt");
		$file_handler{$type_l1."_".$value1."_".$value2}=$report_ids;
	  }
	}

	$total{$type_l1}{'A'}=0;
	$total{$type_l1}{'B'}=0;
    foreach my $value1 ( sort {$a <=> $b} keys %{$overlap_info{$type_l1}} ){
      foreach my $value2 ( sort {$a <=> $b}  keys %{$overlap_info{$type_l1}{$value1}} ){
        $string_to_print .= "|".sizedPrint($value1,30)."|".sizedPrint($value2,30)."|".sizedPrint(scalar(@{$overlap_info{$type_l1}{$value1}{$value2}}),30)."|\n";
        
		# Use proper file handler for outputing IDs
		my  $report_ids = $file_handler{$type_l1."_".$value1."_".$value2};
		
		# file1
		foreach my $array ( @{$overlap_info{$type_l1}{$value1}{$value2}}) {
			my $cpt=0;
			my $last = scalar(@{$array->[0]});
			foreach my $value ( @{$array->[0]} ) { # array0 is id overlarpA
				$cpt++;
				if(! $value){
					print $report_ids "-";
				} else {
					if ($last == $cpt){
						print $report_ids $value;
					} else {
						print $report_ids $value.", ";
					}
				}
			}
			print $report_ids " | ";
			my $cpt2=0;
			my $last2 = scalar(@{$array->[1]});
			foreach my $value ( @{$array->[1]} ) { # array1 is id overlarpB
				$cpt2++;
				if(! $value){
					print $report_ids "-\n";
				} else {
					if ($last2 == $cpt2){
						print $report_ids "$value"
					} else {
						print $report_ids "$value, "
					}
				}
			}
			print $report_ids "\n";
		}


		if ($value1 != 0){
          $total{$type_l1}{'A'} += $value1 * scalar(@{$overlap_info{$type_l1}{$value1}{$value2}});
        }
        if ($value2 != 0){
          $total{$type_l1}{'B'} += $value2 * scalar(@{$overlap_info{$type_l1}{$value1}{$value2}});
        }
      }
    }
    $string_to_print .=  $separator_table;
    $string_to_print .= "Number gene in $filename1: $total{$type_l1}{'A'}\nNumber gene in $filename2: $total{$type_l1}{'B'}\n";
    $string_to_print .= "\n\n"

}
$string_to_print .= "\n";

if ($opt_output){
  print $report $string_to_print;
}
dual_print( $log, $string_to_print);
dual_print( $log, "Bye Bye.\n");
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

# locations are sorted by start position
 # return t1 is location overlap
sub remove_loc_by_id{
	my($list_location, $locusID, $chimere_type, $id)=@_;
    dual_print( $log, "Remove $id from locations \n");
	my @new_list;
	foreach my $locations ( @{$list_location->{$locusID}{$chimere_type}} ){
		if ( lc($locations->[0][3]) ne lc($id) ) {
			push @new_list, $locations;
		}
	}
	@{$list_location->{$locusID}{$chimere_type}} = @new_list;
}

# locations are sorted by start position
# return t1 is location overlap
sub locations_overlap{
	my($list_location1, $list_location2)=@_;
	my $overlap = undef;

	foreach my $location1 ( @{$list_location1} ){
		foreach my $location2 ( @{$list_location2} ){

			if (($location1->[0] <= $location2->[1]) and ($location1->[1] >= $location2->[0])){
				$overlap = 1;
				return $overlap;
			}
			elsif ($location1->[1] < $location2->[0] ){
				last;
			}
			elsif ($location2->[1] < $location1->[0] ){
				next;
			}
		}
	}
	return $overlap;
}

# merge overlapping locations in a list of locations
# lcation [[int, int,id],[int, int,id],[int, int,id]]
sub flatten_locations_and_merge{
  my ($locations1, $locations2, $verbose) = @_;

  my $locations = [@$locations1, @$locations2];
  my @newlocations;
  my $previous_location = undef;
  foreach my $location ( sort {$a->[0] <=> $b->[0]} @{$locations} ){

    # first round
    if (! $previous_location){
       push @newlocations, $location;
       $previous_location = $location;
    }
    # Not first round
    else{
      #  OVERLAP
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
  return \@newlocations ;
}

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

__END__

=head1 NAME

agat_sp_compare_two_annotations.pl

=head1 DESCRIPTION

The script aims to compare two annotation of the same assembly. It provided
information about split/fusion of genes between the two annotations.
The most common cases are:
1 => 0 ( gene uniq to file1)
0 => 1 ( gene uniq to file2)
1 => 1 ( 1 gene from file 1 overlaps only 1 gene from file2)
1 => <many> ( 1 gene from file 1 overlaps <many> genes from file2) => split case (with file 1 as reference)
<many> => 1 ( <many> genes from file 1 overlap only 1 gene from file2) => fusion case (with file 1 as reference)

Then you can get more complex cases:
<many> => <many>  (<many> genes from file 1 overlap <many> genes from file2)

The script output a folder containing a report of number of different cases as well as a file
per case type listing per line the gene feature's ID involved in each case.

=head1 SYNOPSIS

    agat_sp_compare_two_annotations.pl -gff1 infile1.gff -gff2 infile2.gff  [ -o outfile ]
    agat_sp_compare_two_annotations.pl --help

=head1 OPTIONS

=over 8

=item B<--gff1>

Input GTF/GFF file1.

=item B<--gff2>

Input GTF/GFF file2.

=item B<-o> , B<--output> or B<--out> 

Output folder.  It contains a report that resume the type and number of cases, as well as a file per case type 
containing one case per line with the list of gene feature's ID (or other type of feature level1) from file1 then file2 separated by a |.

=item  B<--debug> or B<-d>

Debug option, make it easier to follow what is going on for debugging purpose.

=item  B<--verbose> or B<-v>

Verbose option, make it easier to follow what is going on.

=item B<-c> or B<--config>

String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any, 
otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

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
