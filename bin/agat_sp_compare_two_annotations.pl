#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use List::MoreUtils qw(uniq);
use Sort::Naturally;
use AGAT::AGAT;

my $header = get_agat_header();
my $config;
my $opt_output = undef;
my $gff1 = undef;
my $gff2 = undef;
my $verbose = undef;
my $debug = undef;
my $opt_help= 0;

my @copyARGV=@ARGV;
if ( !GetOptions(
    'c|config=s'               => \$config,
    "h|help"         => \$opt_help,
    "gff1=s"         => \$gff1,
    "gff2=s"         => \$gff2,
	"debug|d!"       => \$debug,
    "verbose|v!"     => \$verbose,
    "output|out|o=s" => \$opt_output))

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
my $report = prepare_fileout($opt_output);

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>     MAIN     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# Activate verbose when debug active
$verbose=1 if ($debug);

######################
### Parse GFF input #
print ("Parsing $gff1\n");
my ($omniscient1, $hash_mRNAGeneLink1) = slurp_gff3_file_JD({ input => $gff1,
                                                              config => $config
                                                              });
print ("\n\nParsing $gff2\n");
my ($omniscient2, $hash_mRNAGeneLink2) = slurp_gff3_file_JD({ input => $gff2,
                                                              config => $config
                                                              });
print ("-- Files parsed --\n");


my $sortBySeq1 = gather_and_sort_l1_location_by_seq_id_and_strand_chimere($omniscient1);
my $sortBySeq2 = gather_and_sort_l1_location_by_seq_id_and_strand_chimere($omniscient2);
print ("GFF3 files sorted\n");

my %overlap_info; # <= Will contain the overlap information
# ------------------------------------------------------------------------------
# ---------------------------- COLLECT ALL LOCATIONS ---------------------------
# ------------------------------------------------------------------------------
print "COLLECT LOCATIONS\n" if ($verbose);
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
print "FLATTEN LOCATIONS\n" if ($verbose);
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
									print "investigate location of $l1_id from $chimere_type from $locusID at $level \n" if ($debug);
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
								print "investigate location of $l1_id from $chimere_type from $locusID at $level \n" if ($debug);
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
print "COMPARE LOCATIONS\n" if ($verbose);
my %seen1;
my %seen2;
foreach my $locusID ( sort  keys %{$flattened_locations1_clean_sorted} ){
  foreach my $chimere_type ( sort keys %{$flattened_locations1_clean_sorted->{$locusID}} ){
		if ( exists_keys ($flattened_locations2_clean_sorted, ($locusID,$chimere_type) ) ){

			#/!\ MUST SORT LOCATIONS
			my @copy1 = @{$flattened_locations1_clean_sorted->{$locusID}{$chimere_type}};
			print scalar(@copy1)." locations 1\n" if ($debug);

			# With foreach the last location is properly handled automatically
			foreach my $locations1 (@copy1){
				if(! exists_keys(\%seen1, ($locations1->[0][3]) ) ) {

					if ( exists_keys ($flattened_locations2_clean_sorted, ($locusID,$chimere_type) ) ) {
						my @copy2 = @{$flattened_locations2_clean_sorted->{$locusID}{$chimere_type}};
						print scalar(@copy2)." locations 2\n" if ($debug);

						# With foreach the last location is properly handled automatically
						foreach my $locations2 (@copy2){
							if(! exists_keys(\%seen2, ($locations2->[0][3] ) ) ){

               	# --------------------- OVERLAP --------------------------------
								if ( locations_overlap($locations1, $locations2) ){
									print "overlap ".$locations1->[0][3]." ".$locations2->[0][3]."\n" if ($debug);

									remove_loc_by_id($flattened_locations1_clean_sorted, $locusID, $chimere_type, $locations1->[0][3]);
									$seen1{$locations1->[0][3]}++;
									remove_loc_by_id($flattened_locations2_clean_sorted, $locusID, $chimere_type, $locations2->[0][3]);
									$seen2{$locations2->[0][3]}++;

									my $flat_overlap_1 = $locations1;
									my $flat_overlap_2 = $locations2;
									my $overlap_A=1;
									my $overlap_B=1;
									my $loop="top";
									my $flip=2;
									my $current_locs;
									my $current_hash;
									my $current_seen;
									my $current_flat=();
									my $current_flat_oppo=();

									while ($flip){
										print "  loop $loop\n" if ($debug);
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
											print "  ".scalar(@{$current_locs})." locations $loop\n" if ($debug);

											#  location   -------------------------
											#  flat                                        -------------------------
											if ( $locations->[scalar(@{$locations})-1][1] < $current_flat_oppo->[0][0] ){
												print "  caseX\n" if ($debug);
												next;
											}
											#  location A                         ----------------
											#  location B  ---------------
											elsif ($current_flat_oppo->[scalar(@{$current_flat_oppo})-1][1] < $locations->[0][0] ){
												print "  caseY\n" if ($debug);
												last;
											}
											elsif ( locations_overlap($current_flat_oppo, $locations) ){
												print "  overlap in overlap ".$current_flat_oppo->[0][3]." ".$locations->[0][3]."\n" if ($debug);
												# keep track of ID that overlap to remove locations later out of the loop
												push @overlap_loc, $locations->[0][3];
												$current_flat = flatten_locations_and_merge($current_flat, $locations);

												if($loop eq "top"){
													$overlap_A++;
												}
												elsif($loop eq "bot"){
													$overlap_B++;
												}
											}
											else{
												print "  caseZ overlap in intron\n" if ($debug);
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
									$overlap_info{$locations1->[0][2]}{$overlap_A}{$overlap_B}++;
								}
		# ----------------------------------- CASE 1 -----------------------------------
								#  location A                         ----------------
								#  location B  ---------------
								elsif ($locations2->[scalar(@{$locations2})-1][1] < $locations1->[0][0] ){

									$overlap_info{$locations1->[0][2]}{0}{1}++;
									print "Case1 notoverlap !\n\n" if ($debug);
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
									print "Case2 notoverlap !\n" if ($debug);
									if(! exists_keys(\%seen1, ( $id1 ) ) ){ # else it has been dealed by overlap case
										$overlap_info{$locations1->[0][2]}{1}{0}++;
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
print "\nLook now what is specific to annotationA \n" if ($verbose);
foreach my $locusID (  keys %{$flattened_locations1_clean_sorted} ){
  foreach my $chimere_type ( keys %{$flattened_locations1_clean_sorted->{$locusID}}){
    foreach my $locations1 ( @{$flattened_locations1_clean_sorted->{$locusID}{$chimere_type}} ){
			$overlap_info{$locations1->[0][2]}{1}{0}++;
			print " Case3 !\n" if ($debug);
    }
  }
}

# ---- NOw deal with what is remaining in annotationB-
# Gather False positive => seq only annotated in annotationB, or type of feature annotated only in annotationB that was missing in annotatoin A
print "\nLook now what is specific to annotationB \n" if ($verbose);
foreach my $locusID (  keys %{$flattened_locations2_clean_sorted} ){
  foreach my $chimere_type ( keys %{$flattened_locations2_clean_sorted->{$locusID}}){
    foreach my $locations2 ( @{$flattened_locations2_clean_sorted->{$locusID}{$chimere_type}} ){
			$overlap_info{$locations2->[0][2]}{0}{1}++;
			print " Case4 !\n" if ($debug);
    }
  }
}

##############
# STATISTICS #
if($verbose){
  print "Compute statistics for $gff1:\n";
	print_omniscient_statistics ({ input => $omniscient1 });

	print "Compute statistics for $gff2:\n";
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


foreach my $type_l1 ( sort keys %overlap_info ){

    $string_to_print .= "$separator_table|".sizedPrint("$type_l1",92)."|\n";
    $string_to_print .= "$separator_table|".sizedPrint($filename1.$ext1,30)."|".sizedPrint($filename2.$ext2,30)."|".sizedPrint("Number of cases",30)."|\n$separator_table";

		$total{$type_l1}{'A'}=0;
		$total{$type_l1}{'B'}=0;
    foreach my $value1 ( sort {$a <=> $b} keys %{$overlap_info{$type_l1}} ){
      foreach my $value2 ( sort {$a <=> $b}  keys %{$overlap_info{$type_l1}{$value1}} ){
        $string_to_print .= "|".sizedPrint($value1,30)."|".sizedPrint($value2,30)."|".sizedPrint($overlap_info{$type_l1}{$value1}{$value2},30)."|\n";
        if ($value1 != 0){
          $total{$type_l1}{'A'} += $value1 * $overlap_info{$type_l1}{$value1}{$value2};
        }
        if ($value2 != 0){
          $total{$type_l1}{'B'} += $value2 * $overlap_info{$type_l1}{$value1}{$value2};
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

# locations are sorted by start position
 # return t1 is location overlap
sub remove_loc_by_id{
	my($list_location, $locusID, $chimere_type, $id)=@_;
	print "Remove $id from locations \n" if ($verbose);
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
The most common case are:
1 => 0 ( gene uniq to file1)
0 => 1 ( gene uniq to file2)
1 => 1 ( 1 gene from file 1 overlaps only 1 gene from file2)
1 => <many> ( 1 gene from file 1 overlaps <many> genes from file2) => split case (with file 1 as reference)
<many> => 1 ( <many> genes from file 1 overlap only 1 gene from file2) => fusion case (with file 1 as reference)

Then you can get more complex cases:
<many> => <many>  (<many> genes from file 1 overlap <many> genes from file2)

=head1 SYNOPSIS

    agat_sp_compare_two_annotations.pl -gff1 infile.gff [ -o outfile ]
    agat_sp_compare_two_annotations.pl --help

=head1 OPTIONS

=over 8

=item B<-gff1>

Input GTF/GFF file1.

=item B<-gff2>

Input GTF/GFF file2.

=item B<-o> , B<--output> or B<--out> 

Output GFF file.  If no output file is specified, the output will be
written to STDOUT.

=item  B<--debug> or B<-d>

Debug option, make it easier to follow what is going on for debugging purpose.

=item  B<--verbose> or B<-v>

Verbose option, make it easier to follow what is going on.

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
