#!/usr/bin/perl -w

package AGAT::OmniscientTool;

use strict;
use warnings;
use Bio::Tools::GFF;
use Bio::Seq;
use Clone 'clone';
use Sort::Naturally;
use Exporter;
use AGAT::Utilities;
use AGAT::OmniscientJson;


our @ISA = qw(Exporter);
our @EXPORT = qw(exists_undef_value is_single_exon_gene get_most_right_left_cds_positions l2_has_cds
l1_has_l3_type check_record_positions l2_identical group_l1IDs_from_omniscient
complement_omniscients rename_ID_existing_in_omniscient keep_only_uniq_from_list2
check_gene_overlap_at_CDSthenEXON location_overlap_update location_overlap nb_feature_level1
check_gene_positions gather_and_sort_l1_location_by_seq_id gather_and_sort_l1_location_by_seq_id_and_strand
gather_and_sort_l1_by_seq_id gather_and_sort_l1_by_seq_id_and_strand extract_cds_sequence group_l1features_from_omniscient
create_omniscient_from_idlevel2list get_feature_l2_from_id_l2_l1 remove_omniscient_elements_from_level2_feature_list
remove_omniscient_elements_from_level2_ID_list featuresList_identik group_features_from_omniscient featuresList_overlap
check_level1_positions check_level2_positions info_omniscient fil_cds_frame
check_all_level1_positions check_all_level2_positions remove_element_from_omniscient
append_omniscient merge_omniscients remove_omniscient_elements_from_level1_id_list
fill_omniscient_from_other_omniscient_level1_id subsample_omniscient_from_level1_id_list check_if_feature_overlap
remove_tuple_from_omniscient create_or_replace_tag remove_element_from_omniscient_attributeValueBased
remove_shortest_isoforms check_gene_overlap_at_level3 gather_and_sort_l1_by_seq_id_for_l2type
gather_and_sort_l1_by_seq_id_for_l1type collect_l1_info_sorted_by_seqid_and_location
remove_l1_and_relatives remove_l2_and_relatives remove_l3_and_relatives
check_mrna_positions check_features_overlap);

sub import {
  AGAT::OmniscientTool->export_to_level(1, @_); # to be able to load the EXPORT functions when direct call; (normal case)
  AGAT::OmniscientTool->export_to_level(2, @_); # to be able to load the EXPORT functions when called from one level up;
}

=head1 SYNOPSIS

This is the code to handle data store in Omniscient.

=head1 DESCRIPTION

 Bench of funciton to fill, modify, create, etc... omniscient data structure.

=head1 AUTHOR

  Jacques Dainat - jacques.dainat@nbis.se

=cut

#				   +------------------------------------------------------+
#				   |+----------------------------------------------------+|
#				   || 			HANDLE OMNISCIENT => Fill / Modify		 ||
#				   |+----------------------------------------------------+|
#				   +------------------------------------------------------+



# omniscient is a hash containing a whole gXf file in memory sorted in a specific way (3 levels)
# If a feature/record already exists in omniscient_to_append, it will be replaced by the new one
# (If the new one content less features, the surnumerary ones are actually erased/removed)
sub fill_omniscient_from_other_omniscient_level1_id {

	my ($level_id_list, $hash_omniscient, $omniscient_to_append)=@_;

	#################
	# == LEVEL 1 == #
	#################
	foreach my $primary_tag_key_level1 (keys %{$hash_omniscient->{'level1'}}){ # primary_tag_key_level1 = gene or repeat etc...
		foreach my $id_tag_key_level1 (@$level_id_list){ #select Id of the list
			if( exists_keys($hash_omniscient, ('level1', $primary_tag_key_level1, $id_tag_key_level1)) ){
				$omniscient_to_append->{'level1'}{$primary_tag_key_level1}{$id_tag_key_level1} = $hash_omniscient->{'level1'}{$primary_tag_key_level1}{$id_tag_key_level1}; # print feature

				#################
				# == LEVEL 2 == #
				#################
				foreach my $primary_tag_key_level2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
					if( exists_keys($hash_omniscient, ('level2', $primary_tag_key_level2, $id_tag_key_level1)) ){
						@{$omniscient_to_append->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1}} = @{$hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1}};

						foreach my $feature_level2 ( @{$hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1}}) {

							#################
							# == LEVEL 3 == #
							#################
							my $level2_ID = lc($feature_level2->_tag_value('ID'));

							# remove feature from omniscient_to_append related to level2_ID. Like that those that are not anymore in hash_omniscient will not appear anymore at the end of the next step (#Now add the new L3 feature of level2_ID)
							foreach my $tag_l3 (keys %{$omniscient_to_append->{'level3'}}){
                if( exists_keys($omniscient_to_append, ('level3', $tag_l3, $level2_ID)) ){
                  delete $omniscient_to_append->{'level3'}{$tag_l3}{$level2_ID} ;
								}
              }
							#Now add the new L3 feature of level2_ID
							foreach my $primary_tag_key_level3 (keys %{$hash_omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...
								if( exists_keys($hash_omniscient, ('level3', $primary_tag_key_level3, $level2_ID)) ){
									@{$omniscient_to_append->{'level3'}{$primary_tag_key_level3}{$level2_ID}} = @{$hash_omniscient->{'level3'}{$primary_tag_key_level3}{$level2_ID}};
								}
							}
						}
					}
				}
			}
		}
	}
}

# @Purpose: Append hash1 with gene and subfeature of hash2 if:
# -> Have CDS and do not overlap at CDS level (if overlap at UTR it's fine) and CDS size is over threshold.
# -> Do not have CDS and do not overlap at exon level
# @input: 4 =>  omniscient1 hash, omniscient2 hash, int, int
# @output 1 =>  omniscient1
sub complement_omniscients {
	my ($omniscient1, $omniscient2, $size_min, $verbose)=@_;

	my %add_omniscient;

	$size_min=0 if (! $size_min);

	if(! $verbose){$verbose=0;}
	my $omniscient1_sorted = gather_and_sort_l1_location_by_seq_id_and_strand($omniscient1);
	my $omniscient2_sorted = gather_and_sort_l1_location_by_seq_id_and_strand($omniscient2);

	foreach my $locusID ( keys %{$omniscient2_sorted}){ # tag_l1 = gene or repeat etc...

		foreach my $tag_l1 ( keys %{$omniscient2_sorted->{$locusID}} ) {

			# Go through location from left to right ### !!
			foreach my $location ( @{$omniscient2_sorted->{$locusID}{$tag_l1}} ) {
				my $id1_l1 = lc($location->[0]);
				print "\nlets look at $id1_l1.\n" if ($verbose >= 3);
				my $take_it=1;

				if( exists_keys($omniscient1_sorted, ($locusID,$tag_l1) ) ) {

					foreach my $location2 ( @{$omniscient1_sorted->{$locusID}{$tag_l1}} ) {
						my $id2_l1 = lc($location2->[0]);

						#If location_to_check start if over the end of the reference location, we stop
						if($location2->[1] > $location->[2]) {last;}
						print "location_to_check start if over the end of the reference location.\n" if ($verbose >= 3);

						#If location_to_check end if inferior to the start of the reference location, we continue next
						if($location2->[2] < $location->[1]) {next;}
						print "location_to_check start if inferior to the start of the reference location.\n" if ($verbose >= 3);

						# Let's check at Gene LEVEL
						if( location_overlap($location, $location2) ){ #location overlap at gene level check now level3
							#let's check at CDS level (/!\ id1_l1 is corresponding to id from $omniscient2)
							if(check_gene_overlap_at_CDSthenEXON($omniscient2, $omniscient1, $id1_l1, $id2_l1)){ #If contains CDS it has to overlap at CDS level, otherwise any type of feature level3 overlaping is sufficient to decide that they overlap
								print "$id2_l1 overlaps $id1_l1, we skip it.\n" if ($verbose >= 3);
								$take_it=undef; last;
							}
							print "$id2_l1 overlaps $id1_l1 overlap but not at CDS level.\n" if ($verbose >= 3);
						}
						else{
							print "$id2_l1 DO NOT OVERLAP $id1_l1.\n" if ($verbose >= 3);
						}
					}
				}

				# We keep it because is not overlaping
				if($take_it){
					print "We take it : $id1_l1\n" if ($verbose >= 3);

					#look at size
					my $still_take_it=undef;
					foreach my $tag_l2 (keys %{$omniscient2->{'level2'}} ){
						if(exists_keys($omniscient2,('level2', $tag_l2, $id1_l1))){
							foreach my $feature_l2 ( @{$omniscient2->{'level2'}{$tag_l2}{$id1_l1}} ){
								my $id_l2 = $feature_l2->_tag_value('ID');

								if(exists_keys($omniscient2,('level3', 'cds', lc($id_l2)))){
									my $cds_size=0;
									foreach my $feature_l3 ( @{$omniscient2->{'level3'}{'cds'}{lc($id_l2)}} ){
										my $size=$feature_l3->end - $feature_l3->start +1;
										$cds_size += $size;
									}
									if($cds_size >= $size_min){
										$still_take_it=1;
										last;
									}
								}
								else{
									$still_take_it=1;
								}

								last if $still_take_it;
							}
							last if $still_take_it;
						}
					}
					# We keep it because has size over threshold
					if ($still_take_it){
						#save level1
						$add_omniscient{'level1'}{$tag_l1}{$id1_l1} = $omniscient2->{'level1'}{$tag_l1}{$id1_l1};
						#save level2
						foreach my $tag_l2 (keys %{$omniscient2->{'level2'}} ){
							if(exists_keys($omniscient2,('level2', $tag_l2, $id1_l1))){
								# Add the level2 list data
								$add_omniscient{'level2'}{$tag_l2}{$id1_l1} = $omniscient2->{'level2'}{$tag_l2}{$id1_l1};
								# for each level2 get the level3 subfeatures
								foreach my $feature_l2 ( @{$omniscient2->{'level2'}{$tag_l2}{$id1_l1}} ){
									my $id_l2 = $feature_l2->_tag_value('ID');
									#save level3
									foreach my $tag_l3 (keys %{$omniscient2->{'level3'}} ){
										if(exists_keys($omniscient2,('level3', $tag_l3, lc($id_l2)))){
											$add_omniscient{'level3'}{$tag_l3}{lc($id_l2)} = $omniscient2->{'level3'}{$tag_l3}{lc($id_l2)};
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

	#Now populate hash1 with data from hash2
	merge_omniscients($omniscient1, \%add_omniscient);

	undef %add_omniscient;

	return $omniscient1;
}

# omniscient is a hash containing a whole gXf file in memory sorted in a specific way (3 levels)
# rename ID in hash_omniscient2 that already exist in hash_omniscient1
sub rename_ID_existing_in_omniscient {

	my ($hash_omniscient1, $hash_omniscient2, $verbose)=@_;

	if(! $verbose){$verbose=1;}

	my $hash_whole_IDs = get_all_IDs($hash_omniscient1);
	my $hash2_whole_IDs = get_all_IDs($hash_omniscient2);

	my %hash_miscCount;
	my $miscCount = \%hash_miscCount;
	my $resume_case=undef;

	#################
	# == LEVEL 1 == #
	#################
	foreach my $tag_l1 (keys %{$hash_omniscient2->{'level1'}}){ # tag_l1 = gene or repeat etc...
		foreach my $id_l1 (keys %{$hash_omniscient2->{'level1'}{$tag_l1}}){
			my $new_parent=undef;
			my $uID = $hash_omniscient2->{'level1'}{$tag_l1}{$id_l1}->_tag_value('ID');

			if ( exists ( $hash_whole_IDs->{$id_l1} ) ){
				$resume_case++;
				my $feature = $hash_omniscient2->{'level1'}{$tag_l1}{$id_l1};
				$uID = replace_by_uniq_ID( $feature, $hash_whole_IDs,  $hash2_whole_IDs, $miscCount);
				$hash_omniscient2->{'level1'}{$tag_l1}{lc($uID)} = delete $hash_omniscient2->{'level1'}{$tag_l1}{$id_l1}; # save feature
				$new_parent=1;
			}
			#################
			# == LEVEL 2 == #
			#################
			foreach my $tag_l2 (keys %{$hash_omniscient2->{'level2'}}){ # tag_l2 = mrna or mirna or ncrna or trna etc...

				if (exists_keys ($hash_omniscient2, ('level2', $tag_l2, $id_l1) ) ){ #Non present in hash2, we create a list with one element

					foreach my $feature_l2 ( @{$hash_omniscient2->{'level2'}{$tag_l2}{$id_l1}}) {

						my $new_parent_l2=undef;

						if($new_parent){
							create_or_replace_tag($feature_l2, 'Parent', $uID);
						}

						my $uID_l2 = $feature_l2->_tag_value('ID');
						my $id_l2 = lc($uID_l2);

						if ( exists ( $hash_whole_IDs->{$id_l2} ) ){

							$resume_case++;
							$uID_l2 = replace_by_uniq_ID($feature_l2, $hash_whole_IDs,  $hash2_whole_IDs, $miscCount);
							$new_parent_l2=1;
						}

						#################
						# == LEVEL 3 == #
						#################
						foreach my $tag_l3 (keys %{$hash_omniscient2->{'level3'}}){

							if (exists_keys ($hash_omniscient2, ('level3', $tag_l3, $id_l2) ) ){

								foreach my $feature_l3 ( @{$hash_omniscient2->{'level3'}{$tag_l3}{$id_l2}}) {

									if($new_parent_l2){
										create_or_replace_tag($feature_l3, 'Parent', $uID_l2);
									}

									my $uID_l3 = $feature_l3->_tag_value('ID');
									my $id_l3 = lc($uID_l3);

									if ( exists ( $hash_whole_IDs->{$id_l2} ) ){
										$resume_case++;
										$uID_l3 = replace_by_uniq_ID($feature_l3, $hash_whole_IDs,  $hash2_whole_IDs, $miscCount);

									}
								}
								#save list feature level3
								if($new_parent_l2){
									$hash_omniscient2->{'level3'}{$tag_l3}{lc($uID_l2)} = delete $hash_omniscient2->{'level3'}{$tag_l3}{lc($id_l2)} ;
								}
							}
						}
					}
					#save list feature level2
					if($new_parent){
						$hash_omniscient2->{'level2'}{$tag_l2}{lc($uID)} = delete $hash_omniscient2->{'level2'}{$tag_l2}{lc($id_l1)};
					}
				}
			}
		}
	}
	print "we renamed $resume_case cases\n" if($verbose and $resume_case);

	return $hash_omniscient2;
}

# put data from hash_omniscient2 in hash_omniscient1
# Features are added even if they are identical. If they have similar name, new name will be given too.
sub merge_omniscients {
	# $hash_omniscient1 = omniscient to append !!!
	my ($hash_omniscient1, $hash_omniscient2, $hash_whole_IDs)=@_;

	if (! $hash_whole_IDs){
		$hash_whole_IDs = get_all_IDs($hash_omniscient1);
	}
	my $hash2_whole_IDs = get_all_IDs($hash_omniscient2);

	my %hash_miscCount;
	my $miscCount = \%hash_miscCount;


  #################
  # ==  HEADER == #
  #################
  if ( exists_keys($hash_omniscient2,('other') ) ){
    foreach my $thing (keys %{$hash_omniscient2->{'other'}}){
      # append new header lines
      if ($thing eq 'header'){
        foreach my $value ( @{$hash_omniscient2->{'other'}{'header'}} ){
          if ( !( grep { $_ eq $value }  @{ $hash_omniscient1->{'other'}{'header'} } ) ){
            push @{$hash_omniscient1->{'other'}{'header'}}, $value; #add value which is new
          }
        }
      }
      # For other thing we take only if no values/key
      else{
        if(! exists_keys($hash_omniscient1,('other', $thing) ) ) {
          $hash_omniscient1->{'other'}{$thing} = clone($hash_omniscient2->{'other'}{$thing});
        }
      }
    }
  }

	#################
	# == LEVEL 1 == #
	#################
	foreach my $tag_l1 (keys %{$hash_omniscient2->{'level1'}}){ # tag_l1 = gene or repeat etc...
		foreach my $id_l1 (keys %{$hash_omniscient2->{'level1'}{$tag_l1}}){

			my $new_parent=undef;
			my $uID = $hash_omniscient2->{'level1'}{$tag_l1}{$id_l1}->_tag_value('ID');

			if ( ! exists_keys ( $hash_whole_IDs,($id_l1) ) ){
					$hash_omniscient1->{'level1'}{$tag_l1}{$id_l1} = delete $hash_omniscient2->{'level1'}{$tag_l1}{$id_l1}; # save feature level1
					$hash_whole_IDs->{$id_l1}++;
			}
			else{
				#print "INFO level1:  Parent $id_l1 already exist. We generate a new one to avoid collision !\n";
				my $feature = $hash_omniscient2->{'level1'}{$tag_l1}{$id_l1};
				$uID = replace_by_uniq_ID( $feature, $hash_whole_IDs,  $hash2_whole_IDs, $miscCount);
				$hash_omniscient1->{'level1'}{$tag_l1}{lc($uID)} = $hash_omniscient2->{'level1'}{$tag_l1}{$id_l1}; # save feature level1
				$new_parent=1;
			}

			#################
			# == LEVEL 2 == #
			#################
			foreach my $tag_l2 (keys %{$hash_omniscient2->{'level2'}}){ # tag_l2 = mrna or mirna or ncrna or trna etc...

				if (exists_keys ($hash_omniscient2, ('level2', $tag_l2, $id_l1) ) ){ #Non present in hash2, we create a list with one element

					foreach my $feature_l2 ( @{$hash_omniscient2->{'level2'}{$tag_l2}{$id_l1}}) {

						my $new_parent_l2=undef;
						if($new_parent){
							create_or_replace_tag($feature_l2, 'Parent', $hash_omniscient1->{'level1'}{$tag_l1}{lc($uID)}->_tag_value('ID'));
						}

						my $uID_l2 = $feature_l2->_tag_value('ID');
						my $id_l2 = lc($uID_l2);

						if ( exists_keys ( $hash_whole_IDs,($id_l2) ) ){

							#print "INFO level2:  Parent $id_l2 already exist. We generate a new one to avoid collision !\n";
							$uID_l2 = replace_by_uniq_ID($feature_l2, $hash_whole_IDs,  $hash2_whole_IDs, $miscCount);
							$new_parent_l2=1;
						}
						else{$hash_whole_IDs->{$id_l2}++;}

						#################
						# == LEVEL 3 == #
						#################
						foreach my $tag_l3 (keys %{$hash_omniscient2->{'level3'}}){

							if (exists_keys ($hash_omniscient2, ('level3', $tag_l3, $id_l2) ) ){

								foreach my $feature_l3 ( @{$hash_omniscient2->{'level3'}{$tag_l3}{$id_l2}}) {

									if($new_parent_l2){
										create_or_replace_tag($feature_l3, 'Parent', $uID_l2);
									}

									my $uID_l3 = $feature_l3->_tag_value('ID');
									my $id_l3 = lc($uID_l3);

									if ( exists_keys ( $hash_whole_IDs,($id_l3) ) ){
									#	print "INFO level3:  Parent $id_l3 already exist. We generate a new one to avoid collision !\n";
										$uID_l3 = replace_by_uniq_ID($feature_l3, $hash_whole_IDs,  $hash2_whole_IDs, $miscCount);
									}
									else{$hash_whole_IDs->{$id_l3}++;}
								}
								#save list feature level3
								$hash_omniscient1->{'level3'}{$tag_l3}{lc($uID_l2)}  = delete $hash_omniscient2->{'level3'}{$tag_l3}{$id_l2} ;
							}
						}
					}
					#save list feature level2
					$hash_omniscient1->{'level2'}{$tag_l2}{lc($uID)} = delete $hash_omniscient2->{'level2'}{$tag_l2}{$id_l1} ;
				}
			}
		}
	}
	return $hash_omniscient1, $hash_whole_IDs;
}

# The hash of reference will be hash_omniscient2. We will keep name from this one.
# When an overlap is found, the ID/parent are fixed
# Check id are not used twice
# sub merge_omniscients_properly {
# 	# $hash_omniscient1 = omniscient to append !!!
# 	my ($hash_omniscient1, $hash_omniscient2)=@_;

# 	my $hash_whole_IDs1 = get_all_IDs($hash_omniscient1);
# 	my $hash_sortBySeq1 = gather_and_sort_l1_location_by_seq_id_and_strand($hash_omniscient1);

# 	my $hash_whole_IDs2 = get_all_IDs($hash_omniscient2);
# 	my $hash_sortBySeq2 = gather_and_sort_l1_location_by_seq_id_and_strand($hash_omniscient2);

# 	my %hash_miscCount;
# 	my $miscCount = \%hash_miscCount;

# 	#################
# 	# == LEVEL 1 == #
# 	#################
# 	foreach my $position1 ( keys %{$hash_sortBySeq1} ){
# 		foreach my $tag1 (keys %{$hash_sortBySeq1->{$position1}}){
# 			if (! exists_keys($hash_sortBySeq2, ($position1, $tag1))){
# 				foreach my $location1 @{$hash_sortBySeq1->{$position1}{tag1}}{
# 					my $l1_id1 = $location1->[0];
# 					my $l1_id_to_use = check_record_ids($l1_id1, $hash_omniscient1, $hash_whole_IDs1, $hash_whole_IDs2, $miscCount);
# 					my @level_id_list = ($l1_id_to_use);
# 					fill_omniscient_from_other_omniscient_level1_id(\@level_id_list, $hash_omniscient1, $hash_omniscient2);
# 				}
# 			}
# 			else{ #check if overlap
# 				# Go through location from left to right ### !!
#       			for (my $i = 0; $i < scalar @{$sortBySeq->{$locusID}{$tag_l1}}; $i++) {
# 			        my $location = shift  @{$sortBySeq->{$locusID}{$tag_l1}};# This location will be updated on the fly
# 			        my $id_l1 = $location->[0];


# 				foreach my $location1 @{$hash_sortBySeq1->{$position1}{tag1}}{
# 					if(check_gene_overlap_at_CDSthenEXON($hash_omniscient1, $hash_omniscient2 , lc($l1_feature1->_tag_value('ID')), lc($id2_l1) )){ #OVERLAP

# 					}
# 					else{ # feature do not overlap
# 						my $l1_id1 = $l1_feature1->_tag_value('ID');
# 						my $l1_id_to_use = check_record_ids($l1_id1, $hash_omniscient1, $hash_whole_IDs1, $hash_whole_IDs2, $miscCount);
# 						my @level_id_list = ($l1_id_to_use);
# 						fill_omniscient_from_other_omniscient_level1_id(\@level_id_list, $hash_omniscient1, $hash_omniscient2);
# 					}

# 				}

# 			}


# 	return $hash_omniscient2;
# }

sub append_omniscient {

	my ($omniscient, $level1,$level2,$level3)=@_;

	foreach my $feature (@$level1){
		my $primaryTag = lc($feature->primary_tag);
		my $id  = lc($feature->_tag_value('ID'));

		if( ! exists_keys($omniscient, ('level1', $primaryTag, $id)) ){
			$omniscient->{"level1"}{$primaryTag}{$id}=$feature;
		}
	}
	foreach my $feature (@$level2){ # if exist, try to append the list
		my $primaryTag = lc($feature->primary_tag);
		my $parent_id  = lc($feature->_tag_value('Parent'));

		if( ! exists_keys($omniscient, ('level2', $primaryTag, $parent_id)) ){
			push(@{$omniscient->{"level2"}{$primaryTag}{$parent_id}}, $feature);###
		}
		else{ # append element in the list if not existing
			my $exist_in_list="no";
			my $id = lc($feature->_tag_value('ID'));

			foreach my $feature_original (@{$omniscient->{"level2"}{$primaryTag}{$parent_id}}){
				my $original_id = lc($feature_original->_tag_value('ID'));
				if ($original_id eq $id){
					$exist_in_list="yes"; last;
				}
			}
			if($exist_in_list eq "no"){ # feature doesnt exist in the feature list already present. So, we append it.
				push(@{$omniscient->{"level2"}{$primaryTag}{$parent_id}}, $feature)
			}
		}
	}
	foreach my $feature (@$level3){
		my $primaryTag = lc($feature->primary_tag);
		my $parent_id = lc($feature->_tag_value('Parent'));

		if( ! exists_keys($omniscient, ('level3', $primaryTag, $parent_id)) ){
			push(@{$omniscient->{"level3"}{$primaryTag}{$parent_id}}, $feature);
		}
		else{ # append element in the list if not existing
			my $exist_in_list="no";
			my $id = lc($feature->_tag_value('ID'));

			foreach my $feature_original (@{$omniscient->{"level3"}{$primaryTag}{$parent_id}}){

				my $original_id = lc($feature_original->get_tag_values('ID'));

				if ($original_id eq $id){
					$exist_in_list="yes"; last;
				}
			}
			if($exist_in_list eq "no"){ # feature doesnt exist in the feature list already present. So, we append it.
				push(@{$omniscient->{"level3"}{$primaryTag}{$parent_id}}, $feature)
			}
		}
	}
}

#				   +------------------------------------------------------+
#				   |+----------------------------------------------------+|
#				   || 			HANDLE OMNISCIENT => Remove				 ||
#				   |+----------------------------------------------------+|
#				   +------------------------------------------------------+

# omniscient is a hash containing a whole gXf file in memory sorted in a specific way (3 levels)
# Input: list of level1 id
#        omniscient
#
sub remove_omniscient_elements_from_level1_id_list {

	my ($hash_omniscient, $level_id_list) = @_  ;

	#################
	# == LEVEL 1 == #
	#################
	foreach my $primary_tag_key_level1 (keys %{$hash_omniscient->{'level1'}}){ # primary_tag_key_level1 = gene or repeat etc...
		foreach my $id_tag_key_level1 (keys %{$hash_omniscient->{'level1'}{$primary_tag_key_level1}}){
			foreach my $level_id (@$level_id_list){
				if($id_tag_key_level1 eq lc($level_id)){

					#################
					# == LEVEL 2 == #
					#################
					foreach my $primary_tag_key_level2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...

						if ( exists ($hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1} ) ){
							foreach my $feature_level2 ( @{$hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1}}) {

								#################
								# == LEVEL 3 == #
								#################
								my $level2_ID = lc($feature_level2->_tag_value('ID'));

								foreach my $primary_tag_key_level3 (keys %{$hash_omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...
									if ( exists ($hash_omniscient->{'level3'}{$primary_tag_key_level3}{$level2_ID} ) ){
										delete $hash_omniscient->{'level3'}{$primary_tag_key_level3}{$level2_ID} # delete level3
									}
								}
							}
							delete $hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1} # delete level2
						}
					}
				delete $hash_omniscient->{'level1'}{$primary_tag_key_level1}{$id_tag_key_level1}; # delete level1
				}
			}
		}
	}
}

# /!\XXX Has to be improved, we should loop over the feature list and extract the id_tag_key_level1 before to loop over the hash
# Input: list of level2 id
#        omniscient
#
sub remove_omniscient_elements_from_level2_feature_list {

	my ($hash_omniscient, $feature_list) = @_  ;

	#################
	# == LEVEL 2 == #
	#################
	foreach my $primary_tag_key_level2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
		foreach my $id_tag_key_level1 (keys %{$hash_omniscient->{'level2'}{$primary_tag_key_level2}}){
			if( exists_keys($hash_omniscient, ('level2', $primary_tag_key_level2, $id_tag_key_level1)) ){
				foreach my $feature_level2 ( @{$hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1}}) {
					my $level2_ID= lc($feature_level2->_tag_value('ID'));

					foreach my $feature (@$feature_list){
						my $feature_ID = lc($feature->_tag_value('ID'));
						my $feature_Parent_ID = lc($feature->_tag_value('Parent'));

						if($level2_ID eq $feature_ID){

							#################
							# == LEVEL 3 == #
							#################
							foreach my $primary_tag_key_level3 (keys %{$hash_omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...
								if( exists_keys($hash_omniscient, ('level3', $primary_tag_key_level3, $level2_ID)) ){
									delete $hash_omniscient->{'level3'}{$primary_tag_key_level3}{$level2_ID} # delete level3
								}
							}
						my @id_concern_list=($feature_Parent_ID);
						my @id_list_to_remove=($feature_ID);
						my @list_tag_key=('all');
						remove_element_from_omniscient(\@id_concern_list, \@id_list_to_remove, $hash_omniscient, 'level2','false', \@list_tag_key);

							if( ! exists_keys($hash_omniscient, ('level2', $primary_tag_key_level2, $id_tag_key_level1)) ){
								#New new list was empty so l2 has been removed, we can now remove l1
								foreach my $primary_tag_key_level1 (keys %{$hash_omniscient->{'level1'}}){
									if( exists_keys($hash_omniscient, ('level1', $primary_tag_key_level1, $feature_Parent_ID)) ){
										delete $hash_omniscient->{'level1'}{$primary_tag_key_level1}{$feature_Parent_ID}
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

# After cleaning if nothing left attache to level 2 we removed it, and the same for level1
sub remove_omniscient_elements_from_level2_ID_list {

	my ($hash_omniscient, $ID_l2_list) = @_  ;

	#################
	# == LEVEL 2 == #
	#################
	foreach my $primary_tag_key_level2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
		foreach my $id_tag_key_level1 (keys %{$hash_omniscient->{'level2'}{$primary_tag_key_level2}}){
			if( exists_keys($hash_omniscient, ('level2', $primary_tag_key_level2, $id_tag_key_level1)) ){
				foreach my $feature_level2 ( @{$hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1}}) {
					my $level2_ID= lc($feature_level2->_tag_value('ID'));
					my $feature_Parent_ID = lc($feature_level2->_tag_value('Parent'));

					foreach my $feature_ID (@$ID_l2_list){
						$feature_ID = lc($feature_ID);

						if($level2_ID eq $feature_ID){

							#################
							# == LEVEL 3 == #
							#################
							foreach my $primary_tag_key_level3 (keys %{$hash_omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...
								if( exists_keys($hash_omniscient, ('level3', $primary_tag_key_level3, $level2_ID)) ){
									delete $hash_omniscient->{'level3'}{$primary_tag_key_level3}{$level2_ID} # delete level3
								}
							}

							my @id_concern_list=($feature_Parent_ID);
							my @id_list_to_remove=($feature_ID);
							my @list_tag_key=('all');
							remove_element_from_omniscient(\@id_concern_list, \@id_list_to_remove, $hash_omniscient, 'level2','false', \@list_tag_key);

							if( ! exists_keys($hash_omniscient, ('level2', $primary_tag_key_level2, $id_tag_key_level1)) ){
								#New list was empty so l2 has been removed, we can now remove l1
								foreach my $primary_tag_key_level1 (keys %{$hash_omniscient->{'level1'}}){
									if( exists_keys($hash_omniscient, ('level1', $primary_tag_key_level1, $feature_Parent_ID)) ){
										delete $hash_omniscient->{'level1'}{$primary_tag_key_level1}{$feature_Parent_ID}
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

# remove value of hash from omniscient in level $level which have the tag incriminated
sub remove_tuple_from_omniscient {

	my ($id_list_to_remove, $hash_omniscient, $level, $bolean, $list_tag_key)=@_;

	# bolean true => we remove if in list_tag_key
	# bolean false => we remove if ti is not in list_tag_key
	my $remove;
  $level = lc($level);
	foreach my $tag_key  (keys %{$hash_omniscient->{$level}}){
		if($bolean eq 'true'){
			$remove="no";
		}else{$remove="yes";}
		foreach my $tag_key_to_match (@$list_tag_key){
			if ((lc($tag_key) eq  lc($tag_key_to_match)) and ($bolean eq 'true')){
				$remove="yes";
			}
			if ((lc($tag_key) eq  lc($tag_key_to_match)) and ($bolean eq 'false')){
				$remove="no";last;
			}
		}
		if ($remove eq 'yes'){
			foreach my $id_key  (keys %{$hash_omniscient->{$level}{$tag_key}}){
				foreach my $id_to_remove (@$id_list_to_remove){
					if(lc($id_to_remove) eq $id_key ){
						delete $hash_omniscient->{ $level }{ $tag_key }{lc($id_to_remove)}; #REMOVE THAT KEY-VALUE pair
					}
				}
			}
		}
	}
}

# from omniscient: remove feature from "feature list" of level2 or level3 with id present in $id_list_to_remove
# $id_concern = ID of parent we will check
sub remove_element_from_omniscient {

	my ($id_concern_list, $id_list_to_remove, $hash_omniscient, $level, $bolean, $list_tag_key)=@_;

	# bolean true => we remove if in list_tag_key
	# bolean false => we remove if is not in list_tag_key
	my $remove;
	#Check level and tag
	foreach my $tag_key  (keys %{$hash_omniscient->{$level}}){
		if($bolean eq 'true'){
			$remove="no";
		}else{$remove="yes";}
		foreach my $tag_key_to_match (@$list_tag_key){

			if ((lc($tag_key) eq  lc($tag_key_to_match)) and ($bolean eq 'true')){
				$remove="yes";
			}
			if ((lc($tag_key) eq  lc($tag_key_to_match)) and ($bolean eq 'false')){
				$remove="no";last;
			}
		}
		#Check feature id from list
		if ($remove eq 'yes'){
			foreach my $id_concern (@$id_concern_list){
				my $mustModifyList=undef;
				my @listok;

				if(exists_keys($hash_omniscient, ($level,$tag_key,lc($id_concern)))){
					foreach my $feature (@{$hash_omniscient->{$level}{$tag_key}{lc($id_concern)}}){
						my $id  = lc($feature->_tag_value('ID'));
						my $shouldremoveit=undef;

						foreach my $id_to_remove (@$id_list_to_remove){

							if(lc ($id_to_remove) eq $id){ # These feature is in list to remove
								$mustModifyList="yes"; $shouldremoveit="yes"; last;
							}
						}
						if(! $shouldremoveit){
							push(@listok, $feature);
						} # Feature not present in id_to_remove, we keep it in list.
					}
					if($mustModifyList){ # at least one feature has been removed from list. Save the new list
						if(@listok){
							@{$hash_omniscient->{$level}{$tag_key}{$id_concern}}=@listok;
						}
						else{ # The list is empty we could remove the key (otherwise we would have saved a emplty list)
							delete $hash_omniscient->{$level}{$tag_key}{$id_concern};
						}
					}
				}
			}
		}
	}
}

# $id_concern = ID of parent we will check
# Go trhough all the element (L1 or L2 list) and check if we find one with the specied tag attribute and value attribute. In that case we remove it from the list
sub remove_element_from_omniscient_attributeValueBased {

	my ($id_concern_list, $attributeValue, $attributeTag, $hash_omniscient, $level, $bolean, $list_tag_key)=@_;

	# bolean true => we remove if in list_tag_key
	# bolean false => we remove if is not in list_tag_key
	my $remove;
	#Check level and tag
	foreach my $tag_key  (keys %{$hash_omniscient->{$level}}){
		if($bolean eq 'true'){
			$remove="no";
		}else{$remove="yes";}
		foreach my $tag_key_to_match (@$list_tag_key){

			if ((lc($tag_key) eq  lc($tag_key_to_match)) and ($bolean eq 'true')){
				$remove="yes";
			}
			if ((lc($tag_key) eq  lc($tag_key_to_match)) and ($bolean eq 'false')){
				$remove="no";last;
			}
		}
		#Check feature id from list
		if ($remove eq 'yes'){
			foreach my $id_concern (@{$id_concern_list}){
				my $mustModifyList=undef;
				my @listok;

				if(exists_keys($hash_omniscient, ($level,$tag_key,lc($id_concern)))){
					foreach my $feature (@{$hash_omniscient->{$level}{$tag_key}{lc($id_concern)}}){
						my $id  = lc($feature->_tag_value('ID'));
						my $shouldremoveit=undef;

						if($feature->has_tag($attributeTag)){
							if( lc($feature->_tag_value($attributeTag)) eq lc($attributeValue) ){
								$mustModifyList="yes"; $shouldremoveit="yes";
							}
							else{push(@listok, $feature);}
						}
						else{push(@listok, $feature);}
					}
					if($mustModifyList){ # at least one feature has been removed from list. Save the new list
						@{$hash_omniscient->{$level}{$tag_key}{$id_concern}}=@listok;
					}
				}
			}
		}
	}
}

# @Purpose: remove from omniscient l1 feature and all subfeatures
# @input: 3 => hash(omniscient hash), feature L1, optional fh to write case removed
# @output: 1 => hash (nb feature removed)
sub remove_l1_and_relatives{
  my ($omniscient, $feature, $fh_removed)=@_;

	my %cases;
	my $cases_l1 = 0; my $cases_l2 = 0; my $cases_l3 = 0;
	my $tag_l1 = lc($feature->primary_tag);
	my $id_l1 = lc($feature->_tag_value('ID'));

  foreach my $ptag_l2 (keys %{$omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...

    if ( exists_keys( $omniscient, ('level2', $ptag_l2, $id_l1) ) ){
      foreach my $feature_l2 ( @{$omniscient->{'level2'}{$ptag_l2}{$id_l1}}) {

        my $level2_ID = lc($feature_l2->_tag_value('ID'));

        foreach my $ptag_l3 (keys %{$omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...
          if ( exists_keys( $omniscient, ('level3', $ptag_l3, $level2_ID) ) ){
            foreach my $feature_l3 ( @{$omniscient->{'level3'}{$ptag_l3}{$level2_ID}}) {
              $cases_l3++;
							print $fh_removed $feature_l3->_tag_value('ID')."\n" if ($fh_removed);
            }
            delete $omniscient->{'level3'}{$ptag_l3}{$level2_ID} # delete level3
          }
        }
        $cases_l2++;
				print $fh_removed $feature_l2->_tag_value('ID')."\n" if($fh_removed);
      }
      delete $omniscient->{'level2'}{$ptag_l2}{$id_l1} # delete level2
    }
  }
  print $fh_removed $feature->gff_string()."\n" if $fh_removed;
  delete $omniscient->{'level1'}{$tag_l1}{$id_l1}; # delete level1
	$cases_l1++;

	$cases{'l1'} = $cases_l1;
	$cases{'l2'} = $cases_l2;
	$cases{'l3'} = $cases_l3;
	$cases{'all'} = $cases_l1+$cases_l2+$cases_l3;

	return \%cases;
}

# @Purpose: remove from omniscient l2 feature and all subfeatures
# @input: 5 => hash(omniscient hash), featureL2,  primary tag l1, id l1, optional fh to write case removed
# @output: 1 => hash (nb feature removed)
sub remove_l2_and_relatives{
  my ($omniscient, $feature, $ptag_l1, $id_l1, $fh_removed)=@_;

  my %cases;
	my $cases_l1 = 0; my $cases_l2 = 0; my $cases_l3 = 0;
	my $ptag_l2 = lc($feature->primary_tag);
  my $level2_Parent_ID = lc($feature->_tag_value('Parent'));
  my $level2_ID = lc($feature->_tag_value('ID'));

  if ( exists_keys( $omniscient, ('level2', $ptag_l2, $id_l1) ) ){ # just extra security in case
    foreach my $feature_l2 ( @{$omniscient->{'level2'}{$ptag_l2}{$id_l1}}) {

      if($level2_ID eq lc($feature_l2->_tag_value('ID')) ){
        # let's delete all l3 subfeatures before to remove the l2
        foreach my $ptag_l3 (keys %{$omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...
          if ( exists_keys( $omniscient, ('level3', $ptag_l3, $level2_ID)  ) ){
            foreach my $feature_l3 ( @{$omniscient->{'level3'}{$ptag_l3}{$level2_ID}}) {
              $cases_l3++;
							print $fh_removed $feature_l3->_tag_value('ID')."\n" if ($fh_removed);
            }
            delete $omniscient->{'level3'}{$ptag_l3}{$level2_ID} # delete level3
          }
        }
      }
    }

    # delete level2 and the hash pointer if the list is empty (no isoform left)
    my @id_concern_list=($id_l1);
    my @id_list_to_remove=($level2_ID);
    my @list_tag_key=('all');
    $cases_l2++;
		print $fh_removed $feature->gff_string()."\n" if ($fh_removed);
    remove_element_from_omniscient(\@id_concern_list, \@id_list_to_remove, $omniscient, 'level2','false', \@list_tag_key);


    if( ! exists_keys($omniscient, ('level2', $ptag_l2, $id_l1) ) ){
      #The list was empty so l2 has been removed, we can now remove l1
      if( exists_keys($omniscient, ('level1', $ptag_l1, $id_l1) ) ){
        $cases_l1++;
				print $fh_removed $id_l1."\n" if ($fh_removed);
        delete $omniscient->{'level1'}{$ptag_l1}{$id_l1};
      }
    }
  }

	$cases{'l1'} = $cases_l1;
	$cases{'l2'} = $cases_l2;
	$cases{'l3'} = $cases_l3;
	$cases{'all'} = $cases_l1+$cases_l2+$cases_l3;

	return \%cases;
}

# @Purpose: remove from omniscient l3 feature and all related feature if needed
# @input: 7 => hash(omniscient hash), feature L3, primary tag l2, id l2, primary tag l2, id l2, optional fh to write case removed
# @output: 1 => hash (nb feature removed)
sub remove_l3_and_relatives{
  my ($omniscient, $feature, $ptag_l1, $id_l1, $ptag_l2, $id_l2, $fh_removed)=@_;

	my %cases;
	my $cases_l1 = 0; my $cases_l2 = 0; my $cases_l3 = 0;
  my $level3_Parent_ID = lc($feature->_tag_value('Parent'));
  my $id_l3 = lc($feature->_tag_value('ID'));
	my $ptag_l3 = lc($feature->primary_tag);

  if ( exists_keys( $omniscient, ('level3', $ptag_l3, $id_l2) ) ){ # just extra security in case
    foreach my $feature_l3 ( @{$omniscient->{'level3'}{$ptag_l3}{$id_l2}}) {

      if($id_l3 eq lc($feature_l3->_tag_value('ID')) ){
        #remove one feature and pointer if no more feature left in the list
        my @id_concern_list=($level3_Parent_ID);
        my @id_list_to_remove=($id_l3);
        my @list_tag_key=('all');
        $cases_l3++;
				print $fh_removed $feature->gff_string()."\n" if ($fh_removed);
        remove_element_from_omniscient(\@id_concern_list, \@id_list_to_remove, $omniscient, 'level3','false', \@list_tag_key);
      }
    }
  }

  # List empty check if we remove l2 or other l3 linked to it
  if( ! exists_keys($omniscient, ('level3', $ptag_l3, $id_l2)) ){
    foreach my $tag_l3 (keys %{$omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...
      if ( exists_keys( $omniscient, ('level3', $tag_l3, $id_l2) ) ){
				$cases{'l1'} = $cases_l1;
				$cases{'l2'} = $cases_l2;
				$cases{'l3'} = $cases_l3;
				$cases{'all'} = $cases_l1+$cases_l2+$cases_l3;

				return \%cases;
      }
    }

    # if we arrive here it means no more L3 feature is attached to L2
    # we remove the L2 parent properly (if isoforms they are kept)
    my @id_concern_list=($id_l1);
    my @id_list_to_remove=($id_l2);
    my @list_tag_key=('all');
    $cases_l2++;
		print $fh_removed $id_l2."\n" if ($fh_removed);
    remove_element_from_omniscient(\@id_concern_list, \@id_list_to_remove, $omniscient, 'level2','false', \@list_tag_key);

    # Should we remove l1 too?
    if( ! exists_keys($omniscient, ('level2', $ptag_l2, $id_l1) ) ){
      #The list was empty so l2 has been removed, we can now remove l1
      if( exists_keys($omniscient, ('level1', $ptag_l1, $id_l1) ) ){
        $cases_l1++;
				print $fh_removed $id_l1."\n" if($fh_removed);
        delete $omniscient->{'level1'}{$ptag_l1}{$id_l1}
      }
    }
  }

	$cases{'l1'} = $cases_l1;
	$cases{'l2'} = $cases_l2;
	$cases{'l3'} = $cases_l3;
	$cases{'all'} = $cases_l1+$cases_l2+$cases_l3;

	return \%cases;
}

#				   +------------------------------------------------------+
#				   |+----------------------------------------------------+|
#				   || 			HANDLE OMNISCIENT => CREATE				 ||
#				   |+----------------------------------------------------+|
#				   +------------------------------------------------------+


# @Purpose: Create an omniscient from list of feature L1,L2 and L3
# @input: 3 =>   list L1, List L2, List L3
# @output 1 =>  omniscient hash reference
sub create_omniscient {

	my ($level1,$level2,$level3)=@_;

	my $omniscient;

	foreach my $feature (@$level1){
		my $id = lc($feature->_tag_value('ID'));
		$omniscient->{"level1"}{lc($feature->primary_tag)}{$id}=$feature;
	}
	foreach my $feature (@$level2){
		my $id = lc($feature->_tag_value('Parent'));
		push(@{$omniscient->{"level2"}{lc($feature->primary_tag)}{$id}}, $feature);###
	}
	foreach my $feature (@$level3){
		my @parentList = lc( $feature->get_tag_values('Parent'));
		foreach my $id (@parentList){
			push(@{$omniscient->{"level3"}{lc($feature->primary_tag)}{$id}}, $feature);
		}
	}
	return $omniscient;
}

# This method will create a new omniscient from an omniscient of reference and a list of id level2
#$list_id_l2 has to be lower case
sub create_omniscient_from_idlevel2list{

	my ($omniscientref, $hash_mRNAGeneLink, $list_id_l2)=@_;

	my %omniscient_new;

	foreach my $id_l2 (@$list_id_l2){
		my  $id_l1 = lc($hash_mRNAGeneLink->{$id_l2});

		# ADD LEVEL1
		foreach my $tag_l1 (keys %{$omniscientref->{'level1'}}){
			if( exists_keys($omniscientref,('level1',$tag_l1,$id_l1) ) ){
				$omniscient_new{'level1'}{$tag_l1}{$id_l1}=$omniscientref->{'level1'}{$tag_l1}{$id_l1};
				last;
			}
		}
		# ADD LEVEL2
		foreach my $tag_l2 (keys %{$omniscientref->{'level2'}}){
			if( exists_keys($omniscientref,('level2',$tag_l2,$id_l1) ) ){
				foreach my $feature_l2 ( @{$omniscientref->{'level2'}{$tag_l2}{$id_l1}}){
					if(lc($feature_l2->_tag_value('ID')) eq $id_l2 ){
						push (@{$omniscient_new{'level2'}{$tag_l2}{$id_l1}}, $feature_l2);
						last;
					}
				}
			}
		}
		# ADD LEVEL3
		foreach my $tag_l3 (keys %{$omniscientref->{'level3'}}){
			if( exists_keys($omniscientref,('level3',$tag_l3,$id_l2) ) ){
				foreach my $feature_l3 ( @{$omniscientref->{'level3'}{$tag_l3}{$id_l2}}){
					push (@{$omniscient_new{'level3'}{$tag_l3}{$id_l2}}, $feature_l3);
				}
			}
		}
	}
	return \%omniscient_new;
}

# @Purpose: filter an omniscient to return a new omnicient containing only data related by the list of level1 IDs
# @input: 1 =>  omniscient hash reference
# @output 1 =>  omniscient hash reference
sub subsample_omniscient_from_level1_id_list {

	my ($hash_omniscient, $level_id_list) = @_  ;

	my %new_hash;

	#################
	# == LEVEL 1 == #
	#################
	foreach my $primary_tag_key_level1 (keys %{$hash_omniscient->{'level1'}}){ # primary_tag_key_level1 = gene or repeat etc...

		foreach my $id_tag_key_level1_raw (@$level_id_list){
			my $id_tag_key_level1 = lc($id_tag_key_level1_raw);
			if(exists ($hash_omniscient->{'level1'}{$primary_tag_key_level1}{$id_tag_key_level1})){

				$new_hash{'level1'}{$primary_tag_key_level1}{$id_tag_key_level1} = delete $hash_omniscient->{'level1'}{$primary_tag_key_level1}{$id_tag_key_level1};

				#################
				# == LEVEL 2 == #
				#################
				foreach my $primary_tag_key_level2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
					if ( exists ($hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1} ) ){
						foreach my $feature_level2 ( @{$hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1}}) {
							my $level2_ID = lc($feature_level2->_tag_value('ID'));

							#################
							# == LEVEL 3 == #
							#################
							foreach my $primary_tag_key_level3 (keys %{$hash_omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...
								if ( exists ($hash_omniscient->{'level3'}{$primary_tag_key_level3}{$level2_ID} ) ){
									$new_hash{'level3'}{$primary_tag_key_level3}{$level2_ID} = delete $hash_omniscient->{'level3'}{$primary_tag_key_level3}{$level2_ID};
								}
							}
						}
						$new_hash{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1} = delete $hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1};
					}
				}
			}
		}
	}
	return \%new_hash;
}

#				   +------------------------------------------------------+
#				   |+----------------------------------------------------+|
#				   || 					Miscenaleous					 ||
#				   |+----------------------------------------------------+|
#				   +------------------------------------------------------+


# INPUT: feature object, String tag, String or Array ref;
# Output: None
sub create_or_replace_tag{

	my ($feature, $tag, $value)=@_;

	if ($feature->has_tag($tag) ) {
			$feature->remove_tag($tag);
			if(ref($value) eq "ARRAY"){
				$feature->add_tag_value($tag,@{$value});
			}
			else{
        		$feature->add_tag_value($tag,$value);
        	}
	}
	else{
		if(ref($value) eq "ARRAY"){
			$feature->add_tag_value($tag,@{$value});
		}
		else{
			$feature->add_tag_value($tag,$value);
		}
	}
}

# frame explanation
# 0 indicates that the feature begins with a whole codon at the 5' most base.
# 1 means that there is one extra base (the third base of a codon) before the first whole codon
# 2 means that there are two extra bases (the second and third bases of the codon) before the first codon.
sub fil_cds_frame {

	my ($hash_omniscient, $db, $verbose)=@_;

	foreach my $primary_tag_key_level2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...

		foreach my $id_tag_key_level1 (keys %{$hash_omniscient->{'level2'}{$primary_tag_key_level2}}) {

			foreach my $feature_level2 ( @{$hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1}}) {

				my $level2_ID = lc($feature_level2->_tag_value('ID'));

				# == LEVEL 3 == #
				if ( exists_keys($hash_omniscient,('level3','cds',$level2_ID) ) ){

					my $strand=$feature_level2->strand;
					my @cds_list;
					if(($feature_level2->strand eq "+") or ($feature_level2->strand eq "1")){
						@cds_list=sort {$a->start <=> $b->start}  @{$hash_omniscient->{'level3'}{'cds'}{$level2_ID}};
					}else{
						@cds_list=sort {$b->start <=> $a->start}  @{$hash_omniscient->{'level3'}{'cds'}{$level2_ID}};
					}

					my $phase = _get_cds_start_phase( $db, $hash_omniscient->{'level3'}{'cds'}{$level2_ID} );

					# Particular case If no phase found and a phase does not exist in the CDS feature we set it to 0 to start
					if ( ! defined( $phase ) and  $cds_list[0]->frame eq "." ) {
						$phase = 0;
						warn "Particular case: No phase found for the CDS start (None in the feature and none can be determined looking at the ORFs)\n".
						"We will assume then to be in phase 0" if $verbose;
					}

					# If no phase found and a phase exists in the CDS feature we keep the original
					# otherwise we loop over CDS features to set the correct phase
          if ( defined( $phase ) ) {
            foreach my $cds_feature ( @cds_list) {
              my $original_phase = $cds_feature->frame;

              if ( ($original_phase eq ".") or ($original_phase != $phase) ){
                print "Original phase $original_phase replaced by $phase for ".$cds_feature->_tag_value("ID")."\n" if $verbose;
                $cds_feature->frame($phase);
              }
  						my $cds_length=$cds_feature->end-$cds_feature->start +1;
  						$phase=(3-(($cds_length-$phase)%3))%3; #second modulo allows to avoid the frame with 3. Instead we have 0.
  					}
          }
				}
			}
		}
	}
}

# @Purpose: get the proper phase of the start of a CDS by looking at the ORF of the different frames
# @input: 1 =>  $db of the fasta genome, list of CDS features, codon table
# @output 1 => integer (0,1,2) or undef
sub _get_cds_start_phase {
  my ($db, $cds_list, $codonTableId) = @_;

  if(! $codonTableId){$codonTableId = 0;}

  my $cds_dna_seq = undef;
  my @cds_list_sorted=sort {$a->start <=> $b->start}  @{$cds_list};
  foreach my $cds_feature ( @cds_list_sorted) {
    $cds_dna_seq .= $db->seq( $cds_feature->seq_id, $cds_feature->start, $cds_feature->end );
  }
  my $cds_obj = Bio::Seq->new(-seq => $cds_dna_seq, -alphabet => 'dna' );
  #Reverse the object depending on strand
  if ($cds_list->[0]->strand == -1 or $cds_list->[0]->strand eq "-"){
    $cds_obj = $cds_obj->revcom();
  }
  my $codonTable = Bio::Tools::CodonTable->new( -id => $codonTableId);
  # case we have start codon => phase 0
  if ($codonTable->is_start_codon(substr($cds_obj->seq, 0, 3)) ) {
    return 0;
  }# No start codon we have to check the phase
  else{
      #try wihtout offset
      my $protein_seq_obj = $cds_obj->translate();
      my $lastChar =  substr $protein_seq_obj->seq(),-1,1;
      my $count = () = $protein_seq_obj->seq() =~ /\*/g;
      if ($lastChar eq "*"){ # if last char is a stop we remove it

        if ($count == 1){
            #print "Missing start codon, phase 0, stop present\n";
            return 0;
        }
      }
      else{
        if ($count == 0){
          #print "Missing start codon, phase 0, missing stop codon\n";
          return 0;
        }
      }

      #try wiht offset (+2 nucleotide)
      $protein_seq_obj = $cds_obj->translate(-offset => 3); #remove 2 nucleotide at the beginning
      $lastChar =  substr $protein_seq_obj->seq(),-1,1;
      $count = () = $protein_seq_obj->seq() =~ /\*/g;
      if ($lastChar eq "*"){ # if last char is a stop we remove it
        if ($count == 1){
            #print "Missing start codon, phase +2, stop present\n";
            return 2;
        }
      }
      else{
        if ($count == 0){
          #print "Missing start codon, phase +2, missing stop codon\n";
          return 2;
        }
      }

      #try wiht offset (+1 nucleotide)
      $protein_seq_obj = $cds_obj->translate(-offset => 2); #remove 2 nucleotide at the beginning
      $lastChar =  substr $protein_seq_obj->seq(),-1,1;
      $count = () = $protein_seq_obj->seq() =~ /\*/g;
      if ($lastChar eq "*"){ # if last char is a stop we remove it
        if ($count == 1){
            #print "Missing start codon, phase +1, stop present\n";
            return 1;
        }
      }
      else{
        if ($count == 0){
          #print "Missing start codon, phase +1, missing stop codon\n";
          return 1;
        }
      }

      # always stop codon in the middle of the sequence... cannot determine correct phase, keep original phase and trow a warning !
      warn "WARNING OmniscientTools _get_cds_start_phase: No phase found for the CDS by looking at the ORFs. ".
      "All frames contain an internal stop codon, thus we cannot determine the correct phase. We will keep original stored phase information.\n";
      return undef;
  }
}

sub info_omniscient {

	my ($hash_omniscient)=@_;

	my %resu;

	foreach my $tag (keys %{$hash_omniscient->{'level1'}}){
    	my $nb=keys %{$hash_omniscient->{'level1'}{$tag}};
    	$resu{$tag}=$nb;
	}

	foreach my $level (keys %{$hash_omniscient}){
  		if ($level ne 'level1' and $level ne 'other'){
    			foreach my $tag (keys %{$hash_omniscient->{$level}}){
      				foreach my $id (keys %{$hash_omniscient->{$level}{$tag}}){
        				my $nb=$#{$hash_omniscient->{$level}{$tag}{$id}}+1;
						if(exists_keys(\%resu,($tag))){
						        $resu{$tag}=$resu{$tag}+$nb;
						}
						else{
							$resu{$tag}=$nb;
						}
     				}
    			}
  		}
	}
	foreach my $tag (keys %resu){
		print "There is $resu{$tag} $tag\n";
	}
}

# omniscient is a hash containing a whole gXf file in memory sorted in a specific way (3 levels)
# This method group all features of a seq_id together.
sub group_features_from_omniscient {

	my ($hash_omniscient) = @_  ;

	my %group;
	#################
	# == LEVEL 1 == #
	#################
	foreach my $primary_tag_key_level1 (keys %{$hash_omniscient->{'level1'}}){ # primary_tag_key_level1 = gene or repeat etc...
		my $key;
		foreach my $id_tag_key_level1 (keys %{$hash_omniscient->{'level1'}{$primary_tag_key_level1}}){
			my $feature_l1=$hash_omniscient->{'level1'}{$primary_tag_key_level1}{$id_tag_key_level1};
			my $seq_id=$feature_l1->seq_id;
			$key="$primary_tag_key_level1$id_tag_key_level1";
			push(@{$group{$seq_id}{$key}}, $hash_omniscient->{'level1'}{$primary_tag_key_level1}{$id_tag_key_level1});
			#################
			# == LEVEL 2 == #
			#################
			foreach my $primary_tag_key_level2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...

				if ( exists ($hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1} ) ){
					foreach my $feature_level2 ( @{$hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1}}) {
						push(@{$group{$seq_id}{$key}}, $feature_level2);
						#################
						# == LEVEL 3 == #
						#################
						my $level2_ID = lc($feature_level2->_tag_value('ID'));

						############
						# THEN ALL THE REST
						foreach my $primary_tag_key_level3 (keys %{$hash_omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...

							if ( exists ($hash_omniscient->{'level3'}{$primary_tag_key_level3}{$level2_ID} ) ){
								foreach my $feature_level3 ( @{$hash_omniscient->{'level3'}{$primary_tag_key_level3}{$level2_ID}}) {
									push(@{$group{$seq_id}{$key}}, $feature_level3);
								}
							}
						}
					}
				}
			}
		}
	}
	return \%group;
}

# omniscient is a hash containing a whole gXf file in memory sorted in a specific way (3 levels)
# This method group all level1 features of a seq_id together.
# hash{seq_id} = @(feature1, feature2 ...)
sub group_l1features_from_omniscient {

	my ($hash_omniscient) = @_  ;

	my %group;
	#################
	# == LEVEL 1 == #
	#################
	foreach my $primary_tag_key_level1 (keys %{$hash_omniscient->{'level1'}}){ # primary_tag_key_level1 = gene or repeat etc...
		foreach my $id_tag_key_level1 (keys %{$hash_omniscient->{'level1'}{$primary_tag_key_level1}}){
			my $feature_l1=$hash_omniscient->{'level1'}{$primary_tag_key_level1}{$id_tag_key_level1};
			my $seq_id=$feature_l1->seq_id;
			push(@{$group{$seq_id}}, $feature_l1);

		}
	}
	return \%group;
}

# omniscient is a hash containing a whole gXf file in memory sorted in a specific way (3 levels)
# This method group all level1 features of a seq_id together.
# hash{seq_id} = @(id1, id2 ...)
sub group_l1IDs_from_omniscient {

	my ($hash_omniscient) = @_  ;

	my %group;
	#################
	# == LEVEL 1 == #
	#################
	foreach my $primary_tag_key_level1 (keys %{$hash_omniscient->{'level1'}}){ # primary_tag_key_level1 = gene or repeat etc...
		foreach my $id_tag_key_level1 (keys %{$hash_omniscient->{'level1'}{$primary_tag_key_level1}}){
			my $feature_l1=$hash_omniscient->{'level1'}{$primary_tag_key_level1}{$id_tag_key_level1};
			my $seq_id=$feature_l1->seq_id;
			push(@{$group{$seq_id}}, lc($feature_l1->_tag_value('ID')));

		}
	}
	return \%group;
}

sub get_feature_l2_from_id_l2_l1 {
	my ($hash_omniscient, $id_l2, $id_l1) = @_  ;
	foreach my $tag_l2 (keys %{$hash_omniscient->{'level2'}}){
		if(exists ($hash_omniscient->{'level2'}{$tag_l2}{$id_l1})){
			foreach my $feature (@{$hash_omniscient->{'level2'}{$tag_l2}{$id_l1}}){
				if ( lc($feature->_tag_value('ID')) eq lc($id_l2) ) {
					return $feature
				}
			}
		}
		else{print "element level2 $tag_l2 $id_l1 doesnt exists in omniscient\n";}
	}
}

#extract sequences form list of cds features in a fasta db
# return a Bio::Seq object
sub extract_cds_sequence {
	my ($feature_list, $db)=@_;

	my $sequence="";
	foreach my $feature (sort {$a->start <=> $b->start} @$feature_list){
		$sequence .= $db->subseq($feature->seq_id,$feature->start,$feature->end);
	}
	my $seq  = Bio::Seq->new( '-format' => 'fasta' , -seq => $sequence);
	if($feature_list->[0]->strand eq "-1" or $feature_list->[0]->strand eq "-"){
		$seq=$seq->revcom;
	}
	return $seq ;
}

# @Purpose: from a omniscient and a gene_id, will get back the extrem value for start and end
# @input: 2 => hash(omniscient), string(gene identifier)
# @output: 2 => integer(extrem start position), integer(extrem end position)
sub get_longest_cds_start_end {
  my  ($hash_omniscient,$gene_id)=@_;
  my $resu_start=100000000000;
  my $resu_end=0;

  #check full CDS for each mRNA
  foreach my $mrna_feature (@{$hash_omniscient->{'level2'}{'mrna'}{lc($gene_id)}}){
    my $mrna_id = lc($mrna_feature->_tag_value('ID'));
    my $extrem_start=100000000000;
    my $extrem_end=0;

    #check all cds pieces
    foreach my $cds_feature (@{$hash_omniscient->{'level3'}{'cds'}{$mrna_id}}){
      if ($cds_feature->start < $extrem_start){
        $extrem_start=$cds_feature->start;
      }
      if($cds_feature->end > $extrem_end){
              $extrem_end=$cds_feature->end ;
      }
    }

    if($extrem_start < $resu_start){
        $resu_start=$extrem_start;
    }
    if($extrem_end > $resu_end){
      $resu_end=$extrem_end;
    }
  }
  return $resu_start,$resu_end;
}

# @Purpose: Filter the RNA to remove short isoforms. Based at CDS level if exists or exon if no CDS
# @input: 1 => hash(omniscient hash)
# @output: return cleaned hash
sub remove_shortest_isoforms{
  my ($hash_omniscient)= @_;

	my @list_to_remove;
	my $case_exon=0;
	my $case_cds=0;

  #################
  # == LEVEL 1 == #
  #################
  foreach my $primary_tag_l1 (keys %{$hash_omniscient->{'level1'}}){ # primary_tag_key_level1 = gene or repeat etc...
    foreach my $id_tag_l1 (keys %{$hash_omniscient->{'level1'}{$primary_tag_l1}}){

      #################
      # == LEVEL 2 == #
      #################
      foreach my $primary_tag_l2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
        if ( exists_keys ($hash_omniscient, ('level2', $primary_tag_l2, $id_tag_l1) ) ){

          #check if there is isoforms
          ###########################

          if ($#{$hash_omniscient->{'level2'}{$primary_tag_l2}{$id_tag_l1}} > 0){

            my $longestL2cds = undef;
						my $longestL2exon = undef;
            my $longestCDSsize = 0;
            my $longestEXONsize = 0;

            foreach my $feature_level2 ( @{$hash_omniscient->{'level2'}{$primary_tag_l2}{$id_tag_l1}}) {

	            my $level2_ID =   lc($feature_level2->_tag_value('ID') ) ;
	            if ( exists_keys( $hash_omniscient, ('level3','cds',$level2_ID ) ) ) {

									my $cdsSize=0;
	                foreach my $cds ( @{$hash_omniscient->{'level3'}{'cds'}{$level2_ID}} ) { # primary_tag_key_level3 = cds or exon or start_codon or utr etc...
	                  $cdsSize += ( $cds->end - $cds->start + 1 );
	                }

	                if($cdsSize > $longestCDSsize ){
										if($longestL2cds){ # we found a longest CDS. The previous is shortest
											$case_cds++;
											push @list_to_remove, [$longestL2cds, $primary_tag_l1, $id_tag_l1];
										}
	                  $longestL2cds = $feature_level2;
	                  $longestCDSsize = $cdsSize;
	                }
									else{ # we have a longest CDS. The current is shortest
										$case_cds++;
										push @list_to_remove, [$feature_level2, $primary_tag_l1, $id_tag_l1];
									}
	            }
	            elsif ( exists_keys( $hash_omniscient, ('level3','exon',$level2_ID ) ) ) {

								if ($longestL2cds){
									# We have a CDS for another isoform so we remove this one that do not have CDS
									push @list_to_remove, [ $feature_level2, $primary_tag_l1, $id_tag_l1];
									$case_exon++;
								}
								else{
									my $exonSize=0;
	                foreach my $exon ( @{$hash_omniscient->{'level3'}{'exon'}{$level2_ID}} ) { # primary_tag_key_level3 = cds or exon or start_codon or utr etc...
	                  $exonSize += ( $exon->end - $exon->start + 1 );
	                }
	                if($exonSize > $longestEXONsize ){
										if($longestL2exon){ # we found a longest exon. The previous is shortest
											push @list_to_remove, [$longestL2exon, $primary_tag_l1, $id_tag_l1];
											$case_exon++;
										}
	                  $longestL2exon = $feature_level2;
	                  $longestEXONsize = $exonSize;
	                }
									else{ # we have a longest exons. The current is shortest
										push @list_to_remove, [ $feature_level2, $primary_tag_l1, $id_tag_l1 ];
										$case_exon++;
									}
								}
            	}
            }
          }
        }
      }
    }
  }
	# remove listed l2
	foreach my $infos (@list_to_remove) {
		my $cases = remove_l2_and_relatives( $hash_omniscient, @$infos);
	}

  return $case_cds, $case_exon;
}

# @Purpose: Counter the number of feature level in an omniscient
# @input: 1 => hash(omniscient hash)
# @output: integer
sub nb_feature_level1 {

  my ($omniscient)=@_;
  my $resu=0;
	foreach my $tag_level1 (keys %{$omniscient->{'level1'}}){
		$resu += (keys %{$omniscient->{'level1'}{$tag_level1}})
	}
	return $resu;
}

# @Purpose: get all the ID present in an omniscient
# @input: 1 => hash(omniscient hash)
# @output: hash of the whole IDs
sub get_all_IDs{
	my ($omniscient)=@_;

	my %whole_IDs;

	#################
	# == LEVEL 1 == #
	#################
	foreach my $primary_tag_l1 (keys %{$omniscient->{'level1'}}){ # primary_tag_l1 = gene or repeat etc...
		foreach my $id_l1 (keys %{$omniscient->{'level1'}{$primary_tag_l1}}){
			$whole_IDs{$id_l1}++;
		}
	}
	#################
	# == LEVEL 2 == #
	#################
	foreach my $primary_tag_l2 (keys %{$omniscient->{'level2'}}){ # primary_tag_l2 = mrna or mirna or ncrna or trna etc...
		foreach my $id_l1 ( keys %{$omniscient->{'level2'}{$primary_tag_l2}}) {
			foreach my $feature_level2 ( @{$omniscient->{'level2'}{$primary_tag_l2}{$id_l1}}) {
				my $level2_ID  = lc($feature_level2->_tag_value('ID'));
				$whole_IDs{$level2_ID}++;
			}
		}
	}
	#################
	# == LEVEL 3 == #
	#################
	foreach my $primary_tag_l3 (keys %{$omniscient->{'level3'}}){ # primary_tag_l3 = cds or exon or start_codon or utr etc...
		foreach my $level2_ID ( keys %{$omniscient->{'level3'}{$primary_tag_l3}}) {
			foreach my $feature_level3 ( @{$omniscient->{'level3'}{$primary_tag_l3}{$level2_ID}}) {
				my $level3_ID  = lc($feature_level3->_tag_value('ID'));
				$whole_IDs{$level3_ID}++;
			}
		}
	}
	return \%whole_IDs;
}

# @Purpose: Replace ID by Uniq ID and modify all parent attribute of child feature to stay in line with the modification
# @input: 4 => feature objetc, hash of ids, hash of ids, hash of feature counted to give more rapidly a name
# @output: uniq ID
sub replace_by_uniq_ID{
	my ($feature, $hash_whole_IDs, $hash2_whole_IDs, $miscCount) = @_;

	my $id = $feature->_tag_value('ID');
	my $prefix = "IDmodified";
	my $key;

	if($prefix){
		$key=$prefix."-".lc($feature->primary_tag);
	}
	else{
		$key=lc($feature->primary_tag);
	}

	my $uID=$id;
	while( exists_keys($hash_whole_IDs, (lc($uID)) ) or exists_keys($hash2_whole_IDs, (lc($uID)) ) ){	 #loop until we found an uniq tag
		$miscCount->{$key}++;
		$uID = $key."-".$miscCount->{$key};
	}

	#push the new ID
	$hash_whole_IDs->{lc($uID)}=$id;

	# modify the feature ID with the correct one chosen
	create_or_replace_tag($feature,'ID', $uID); #modify ID to replace by parent value

	#Now repercute this modification to the subfeatures
	return $uID;
}

#				   +------------------------------------------------------+
#				   |+----------------------------------------------------+|
#				   || 			MANIPULATION AT OMNISCIENT LEVEL1/2/3	 ||
#				   |+----------------------------------------------------+|
#				   +------------------------------------------------------+

# @Purpose: Check 2 lists of feature L2 and remove the identical ones from the second list.
# @input: 4 =>  omniscient Hash reference, list1 reference of L2 features,  list2 reference of L2 features, verbose option for debug
# @output: list2 minus all the feature identical to one of the list1 feature
sub keep_only_uniq_from_list2{
	my ($omniscient, $list1_l2, $list2_l2, $verbose)= @_;

	my @new_list2;
	my $keep = 1;

	foreach my $feature2 ( @{$list2_l2} ){
		foreach my $feature1 ( @{$list1_l2} ){
			if(l2_identical($omniscient, $feature1, $feature2, $verbose )){
				$keep = undef; last;
			}
		}
		if($keep){
			push(@new_list2, $feature2);
		}
		else{ # We dont keep the l2 feature so we have to remove all related features
			remove_l2_related_feature($omniscient, $feature2, $verbose);
		}
	}
	return \@new_list2;
}

# check if l2 are identical. So look recursively at the level under.
# return 1 if identical
sub l2_identical{
	my ($omniscient, $feature1_l2, $feature2_l2, $verbose)= @_;
	my $result=1;

	my $id1_l2 = lc($feature1_l2->_tag_value('ID') );
	my $id2_l2 = lc($feature2_l2->_tag_value('ID') );

	foreach my $l3_type (keys %{$omniscient->{'level3'}} ){
		if(exists_keys($omniscient,('level3', $l3_type, $id1_l2))){
			if(exists_keys($omniscient,('level3', $l3_type, $id2_l2))){

				if(scalar @{$omniscient->{'level3'}{$l3_type}{$id1_l2}} ==  scalar @{$omniscient->{'level3'}{$l3_type}{$id2_l2}}){

					foreach my $feature1_level3 ( sort {$a->start <=> $b->start} @{$omniscient->{'level3'}{$l3_type}{$id1_l2}}) {

						my $identik = undef;
						foreach my $feature2_level3 ( sort {$a->start <=> $b->start} @{$omniscient->{'level3'}{$l3_type}{$id2_l2}}) {

							if( ($feature1_level3->start == $feature2_level3->start) and ($feature1_level3->end == $feature2_level3->end) ){
								$identik=1;
							}
						}
						if(! $identik){
							return undef;
						}
					}
				}
				else{return undef;} # Not same number of features. Cannot be identical
			}
			else{return undef;}
		}
		else{
			if(exists_keys($omniscient,('level3', $l3_type, $id2_l2))){ # $id1_l2 do not have l3 but $id2_l2 has !
				return undef;
			}
		}
	}
	print "The isoforms $id1_l2 and $id2_l2 are identical\n" if ($verbose >= 2 and $result);
    return $result;
}

#
#
#Remove everything related to l2. But not itself ... why ?
sub remove_l2_related_feature{
	my ($omniscient, $feature2, $verbose) = @_;

	my $l1_id = lc($feature2->_tag_value('Parent'));
	my $l2_id = lc($feature2->_tag_value('ID'));

	#remove level 1 feature
	foreach my $tag (keys %{$omniscient->{'level1'}}){
		if(exists_keys($omniscient, ('level1', $tag, $l1_id))){
			delete $omniscient->{'level1'}{$tag}{$l1_id};
			last;
		}
	}
	foreach my $tag (keys %{$omniscient->{'level3'}}){
		if(exists_keys($omniscient, ('level3', $tag, $l2_id))){
			delete $omniscient->{'level3'}{$tag}{$l2_id};
		}
	}
}




# @Purpose: Copy past an attribute and change its tag
# @input: 3 =>  omniscient Hash reference, String = attribute tag original, String = attribute tag new
# @output none =>  The hash itself is modified
sub create_attribute_from_existing_attribute{
	my ($omniscient, $original_attribute, $new_attribute) = @_;

	foreach my $primary_tag_l1 (keys %{$omniscient->{'level1'}}){ # primary_tag_l1 = gene or repeat etc...
		foreach my $id_l1 (keys %{$omniscient->{'level1'}{$primary_tag_l1}}){
			my $feature_l1 = $omniscient->{'level1'}{$primary_tag_l1}{$id_l1};

			if( $feature_l1->has_tag($original_attribute)  and ! $feature_l1->has_tag($new_attribute)  ) {
				create_or_replace_tag($feature_l1,$new_attribute, $feature_l1->get_tag_values($original_attribute));
			}

			#################
			# == LEVEL 2 == #
			#################
			foreach my $primary_tag_l2 (keys %{$omniscient->{'level2'}}){ # primary_tag_l2 = mrna or mirna or ncrna or trna etc...
				if ( exists ($omniscient->{'level2'}{$primary_tag_l2}{$id_l1} ) ){
					foreach my $feature_l2 ( @{$omniscient->{'level2'}{$primary_tag_l2}{$id_l1}}) {

						if( $feature_l2->has_tag($original_attribute)  and ! $feature_l2->has_tag($new_attribute) ){
							create_or_replace_tag($feature_l2,$new_attribute, $feature_l2->get_tag_values($original_attribute));
						}
						#################
						# == LEVEL 3 == #
						#################
						my $level2_ID  = lc($feature_l2->_tag_value('ID'));

						foreach my $primary_tag_l3 (keys %{$omniscient->{'level3'}}){ # primary_tag_l3 = cds or exon or start_codon or utr etc...
							if ( exists ($omniscient->{'level3'}{$primary_tag_l3}{$level2_ID} ) ){
								foreach my $feature_l3 ( @{$omniscient->{'level3'}{$primary_tag_l3}{$level2_ID}}) {
									if( $feature_l3->has_tag($original_attribute) and ! $feature_l3->has_tag($new_attribute) ){
										create_or_replace_tag($feature_l3, $new_attribute, $feature_l3->get_tag_values($original_attribute));
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

# @Purpose: Create a locus tag for all feature L1 and share it with all children features
# @input: 3 =>  omniscient Hash reference, String = attribute tag original to  use as locus tag, String = locus_tag attribute tag
# @output none =>  The hash itself is modified
sub create_locus_tag{
	my ($omniscient, $original_attribute, $locus_tag) = @_;

	foreach my $primary_tag_l1 (keys %{$omniscient->{'level1'}}){ # primary_tag_l1 = gene or repeat etc...
		foreach my $id_l1 (keys %{$omniscient->{'level1'}{$primary_tag_l1}}){
			my $feature_l1 = $omniscient->{'level1'}{$primary_tag_l1}{$id_l1};

			my $locus_tag_value=undef;
			if( $feature_l1->has_tag($original_attribute)  and ! $feature_l1->has_tag($locus_tag)  ) {
				$locus_tag_value =  $feature_l1->_tag_value($original_attribute);
				create_or_replace_tag($feature_l1, $locus_tag, $locus_tag_value);
			}
			else{
				$locus_tag_value =  $feature_l1->_tag_value($locus_tag);
			}

			#################
			# == LEVEL 2 == #
			#################
			foreach my $primary_tag_l2 (keys %{$omniscient->{'level2'}}){ # primary_tag_l2 = mrna or mirna or ncrna or trna etc...
				if ( exists ($omniscient->{'level2'}{$primary_tag_l2}{$id_l1} ) ){
					foreach my $feature_l2 ( @{$omniscient->{'level2'}{$primary_tag_l2}{$id_l1}}) {

						if( $feature_l2->has_tag($original_attribute)  and ! $feature_l2->has_tag($locus_tag) ){
							create_or_replace_tag($feature_l2,$locus_tag, $locus_tag_value);
						}
						#################
						# == LEVEL 3 == #
						#################
						my $level2_ID  = lc($feature_l2->_tag_value('ID'));

						foreach my $primary_tag_l3 (keys %{$omniscient->{'level3'}}){ # primary_tag_l3 = cds or exon or start_codon or utr etc...
							if ( exists ($omniscient->{'level3'}{$primary_tag_l3}{$level2_ID} ) ){
								foreach my $feature_l3 ( @{$omniscient->{'level3'}{$primary_tag_l3}{$level2_ID}}) {
									if( $feature_l3->has_tag($original_attribute) and ! $feature_l3->has_tag($locus_tag) ){
										create_or_replace_tag($feature_l3, $locus_tag, $locus_tag_value);
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




#				   +------------------------------------------------------+
#				   |+----------------------------------------------------+|
#				   || 			FEATURES LOCATION       				 ||
#				   |+----------------------------------------------------+|
#				   +------------------------------------------------------+

# looking the end and the start, the method check if two location overlap.
#A location os [Id, position1, position2]
# return t1 is location overlap
sub location_overlap{
	my($location1, $location2)=@_;
	my $overlap = undef;

	if (($location1->[1] <= $location2->[2]) and ($location1->[2] >= $location2->[1])){
		$overlap = 1;
	}

	return $overlap;
}

# looking the end and the start, the method check if two location overlap.
# A location iss [Id, position1, position2]
# return the intersect of locations
sub location_overlap_update{
	my($location1, $location2)=@_;
	my $location = $location1;
	my $overlap = undef;

	if (($location1->[1] <= $location2->[2]) and ($location1->[2] >= $location2->[1])){
		$overlap = 1;
		if($location2->[1] < $location1->[1]){
			$location->[1] = $location2->[1]
		}
		if($location2->[2] > $location1->[2]){
			$location->[2] = $location2->[2]
		}
	}

	return $location, $overlap;
}

# Check if two genes have at least one L2 isoform which overlap at level3 feature
# Ouput: 1 => String corresponding to the tag of the l3 feature
sub check_gene_overlap_at_level3{
  my  ($hash_omniscient, $hash_omniscient2, $gene_id, $gene_id2, $tag_l3)=@_;

	$tag_l3 = lc($tag_l3); # lowercase in case

	foreach my $l2_type (keys %{$hash_omniscient->{'level2'}} ){

		#check full CDS for each mRNA
		if(exists_keys($hash_omniscient,('level2', $l2_type, lc($gene_id)))){

			foreach my $mrna_feature (@{$hash_omniscient->{'level2'}{$l2_type}{lc($gene_id)}}){
				my $mrna_id1 = $mrna_feature->_tag_value('ID');

				if(exists_keys($hash_omniscient2,('level2', $l2_type, lc($gene_id2)))){

			    foreach my $mrna_feature2 (@{$hash_omniscient2->{'level2'}{$l2_type}{lc($gene_id2)}}){ # from here bothe feature level2 are the same type

						my $mrna_id2 = $mrna_feature2->_tag_value('ID');

				    #check l3 pieces
				    if(exists_keys($hash_omniscient,('level3', $tag_l3, lc($mrna_id1)))){
			      	if(exists_keys($hash_omniscient2,('level3', $tag_l3, lc($mrna_id2)))){
					    	foreach my $cds_feature1 (@{$hash_omniscient->{'level3'}{$tag_l3}{lc($mrna_id1)}}){
					        foreach my $cds_feature2 (@{$hash_omniscient2->{'level3'}{$tag_l3}{lc($mrna_id2)}}){

					        	if(($cds_feature2->start <= $cds_feature1->end) and ($cds_feature2->end >= $cds_feature1->start )){ # they overlap
					            return $tag_l3;
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
	return undef;
}

# Check if two genes have at least one L2 isoform which overlap at cds level.
# if no CDS we check if overlap at any other l3 feature.
# Check exon level only if no CDS exists !
sub check_gene_overlap_at_CDSthenEXON{
  my  ($hash_omniscient, $hash_omniscient2, $gene_id, $gene_id2)=@_;

	foreach my $l2_type (keys %{$hash_omniscient->{'level2'}} ){

		#check full CDS for each mRNA
		if(exists_keys($hash_omniscient,('level2', $l2_type, lc($gene_id)))){
			foreach my $mrna_feature (@{$hash_omniscient->{'level2'}{$l2_type}{lc($gene_id)}}){
				my $mrna_id1 = $mrna_feature->_tag_value('ID');

				if(exists_keys($hash_omniscient2,('level2', $l2_type, lc($gene_id2)))){

          # from here both feature level2 are the same type
			    foreach my $mrna_feature2 (@{$hash_omniscient2->{'level2'}{$l2_type}{lc($gene_id2)}}){
					  my $mrna_id2 = $mrna_feature2->_tag_value('ID');

				    #check all cds pieces - CDS against CDS
				    if(exists_keys($hash_omniscient,('level3', 'cds', lc($mrna_id1))) and
				      exists_keys($hash_omniscient2,('level3', 'cds', lc($mrna_id2))) ) {

						  foreach my $cds_feature1 (@{$hash_omniscient->{'level3'}{'cds'}{lc($mrna_id1)}}){
					      foreach my $cds_feature2 (@{$hash_omniscient2->{'level3'}{'cds'}{lc($mrna_id2)}}){
					        if(($cds_feature2->start <= $cds_feature1->end) and ($cds_feature2->end >= $cds_feature1->start )){ # they overlap
					          return "cds";
					        }
					      }
				      }

				    }# CDS not in both, check at CDS / exon / match level only if same level2 type
				    else{
				    	foreach my $tag_l3 (keys %{$hash_omniscient->{'level3'}}){

				    		if(exists_keys($hash_omniscient,('level3', $tag_l3, lc($mrna_id1)))){
				    			foreach my $feature1 (@{$hash_omniscient->{'level3'}{$tag_l3}{lc($mrna_id1)}}){

										if(exists_keys($hash_omniscient2,('level3', $tag_l3, lc($mrna_id2)))){
				    					foreach my $feature2 (@{$hash_omniscient2->{'level3'}{$tag_l3}{lc($mrna_id2)}}){

				    						if(($feature2->start <= $feature1->end) and ($feature2->end >= $feature1->start )){ # they overlap
					            		return "exon";
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
		}
	}
  return undef;
}

# @Purpose: Check the start and end of gene feature based on its mRNAs and eventualy fix it.
# @input: 2 => hash(omniscient hash), string(gene identifier)
# @output: none
sub check_gene_positions {

  my ($hash_omniscient, $gene_id)=@_;

  #####
  #Modify gene start-end (have to check size of each mRNA)
  my $geneExtremStart=1000000000000;
  my $geneExtremEnd=0;
  foreach my $primary_tag_l2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
  	if (exists_keys($hash_omniscient, ('level2', $primary_tag_l2, lc($gene_id) ) ) ){ # check if they have mRNA avoiding autovivifcation
	    foreach my $mrna_feature ( @{$hash_omniscient->{'level2'}{$primary_tag_l2}{lc($gene_id)}}) {
	      	my $start=$mrna_feature->start();
	      	my $end=$mrna_feature->end();

	      if ($start < $geneExtremStart){
	        $geneExtremStart=$start;
	      }
	      if($end > $geneExtremEnd){
	        $geneExtremEnd=$end;
	      }
	    }
	}
  }
  my $gene_feature=$hash_omniscient->{'level1'}{'gene'}{lc($gene_id)};
  if ($gene_feature->start != $geneExtremStart){
      $gene_feature->start($geneExtremStart);
   }
  if($gene_feature->end != $geneExtremEnd){
    $gene_feature->end($geneExtremEnd);
  }
}

# Check the start and end of level1 feature based on all features level2;
sub check_all_level1_positions {
	my ($args) = @_;

	my $resume_case=undef;
	# -------------- INPUT --------------
	# Check we receive a hash as ref
	if(ref($args) ne 'HASH'){ warn "Hash Arguments expected for check_all_level1_positions. Please check the call.\n";exit;	}
	# -- Declare all variables and fill them --
	my ($hash_omniscient, $verbose, $log);
	if( defined($args->{omniscient})) {$hash_omniscient = $args->{omniscient};} else{ print "Input omniscient mandatory to use check_all_level1_positions!"; exit; }
	if( defined($args->{verbose}) ) { $verbose = $args->{verbose}; } else { $verbose = 0;}
	if( defined($args->{log}) ) { $log = $args->{log}; }

	foreach my $tag_l1 (keys %{$hash_omniscient->{'level1'}}){ # primary_tag_key_level1 = gene or repeat etc...
		foreach my $id_l1 ( keys %{$hash_omniscient->{'level1'}{$tag_l1}} ) { #sort by position

			my $level1_feature = $hash_omniscient->{'level1'}{$tag_l1}{$id_l1};

			$resume_case++ if(check_level1_positions({ omniscient => $hash_omniscient,
																								 feature => $level1_feature,
																								 verbose => $verbose}));
		}
	}

	if($resume_case){
		dual_print($log, "We fixed $resume_case wrong level1 location cases\n", $verbose );
	}
	else{
		dual_print($log, "No problem found\n", $verbose );
	}
}

# Purpose: review all the feature L2 to adjust their start and stop according to the extrem start and stop from L3 sub features.
sub check_all_level2_positions{
	my ($args) = @_;
	my $resume_case=undef;
	# -------------- INPUT --------------
	# Check we receive a hash as ref
	if(ref($args) ne 'HASH'){ warn "Hash Arguments expected for check_all_level1_positions. Please check the call.\n";exit;	}
	# -- Declare all variables and fill them --
	my ($hash_omniscient, $verbose, $log);
	if( defined($args->{omniscient})) {$hash_omniscient = $args->{omniscient};} else{ print "Input omniscient mandatory to use check_all_level1_positions!"; exit; }
	if( defined($args->{verbose}) ) { $verbose = $args->{verbose}; } else { $verbose = 0;}
	if( defined($args->{log}) ) { $log = $args->{log}; }

	foreach my $tag_l1 (keys %{$hash_omniscient->{'level1'}}){ # primary_tag_key_level1 = gene or repeat etc...
		foreach my $id_l1 ( keys %{$hash_omniscient->{'level1'}{$tag_l1}} ) { #sort by position

			foreach my $tag_level2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
					if ( exists_keys ($hash_omniscient, ('level2', $tag_level2, $id_l1) ) ){

					foreach my $mRNA_feature ( @{$hash_omniscient->{'level2'}{$tag_level2}{$id_l1}}){
						my $level2_ID = lc($mRNA_feature->_tag_value('ID'));
						my @feature_list=();
						foreach my $primary_tag_l3 ( keys %{$hash_omniscient->{'level3'}}){ # primary_tag_l3 = cds or exon or start_codon or utr etc...

							if ( exists_keys( $hash_omniscient, ('level3', $primary_tag_l3, $level2_ID) ) ){
								push @feature_list, @{$hash_omniscient->{'level3'}{$primary_tag_l3}{$level2_ID}};
							}
						}
						if( @feature_list ){ #could be emtpy like in match match_part features, so avoid this cases
							 $resume_case++ if( check_mrna_positions({ l2_feature => $mRNA_feature,
							 																					exon_list => \@feature_list,
																												log => $log,
																												verbose => $verbose} ) );
						}
					}
				}
			}
		}
	}
	if($resume_case){
		dual_print($log, "We fixed $resume_case wrong level2 location cases\n", $verbose );
	}
	else{
		dual_print($log, "No problem found\n", $verbose );
	}
}

# Check the start and end of mRNA based a list of feature like list of exon;
sub check_mrna_positions{
	my ($args) = @_;
	my $result=undef;
	# -------------- INPUT --------------
	# Check we receive a hash as ref
	if(ref($args) ne 'HASH'){ warn "Hash Arguments expected for check_mrna_positions. Please check the call.\n";exit;	}
	# -- Declare all variables and fill them --
	my ($mRNA_feature, $exon_list, $verbose, $log);
	if( defined($args->{l2_feature})) {$mRNA_feature = $args->{l2_feature};} else{ print "Input l2_feature mandatory to use check_mrna_positions!"; exit; }
	if( defined($args->{exon_list})) {$exon_list = $args->{exon_list};} else{ print "Input exon_list mandatory to use check_mrna_positions!"; exit; }
	if( defined($args->{verbose}) ) { $verbose = $args->{verbose}; } else { $verbose = 0;}
	if( defined($args->{log}) ) { $log = $args->{log}; }

	my @exon_list_sorted = sort {$a->start <=> $b->start} @{$exon_list};
	my $exonStart=$exon_list_sorted[0]->start;

	@exon_list_sorted = sort {$a->end <=> $b->end} @exon_list_sorted;
	my $exonEnd=$exon_list_sorted[$#exon_list_sorted]->end;

	#check start
	if ($mRNA_feature->start != $exonStart){
		dual_print($log, "We modified the L2 LEFT extremity for the sanity the biological data!\n", 0); # print log only
		$mRNA_feature->start($exonStart);
		$result=1;
	}
	#check stop
	if($mRNA_feature->end != $exonEnd){
		dual_print($log, "We modified the L2 RIGHT extremity for the sanity the biological data!\n", 0); # print log only
		$mRNA_feature->end($exonEnd);
		$result=1;
	}

	return $result;
}

# Check the start and end of level1 feature based on all features level2;
#return 1 if something modified
sub check_level1_positions {
	my ($args) = @_;
	my $result=undef;

	# -------------- INPUT --------------
	# Check we receive a hash as ref
	if(ref($args) ne 'HASH'){ warn "Hash Arguments expected for check_level1_positions. Please check the call.\n";exit;	}

	my ($hash_omniscient, $feature_l1, $verbose, $log);
	if( defined($args->{omniscient})) {$hash_omniscient = $args->{omniscient};} else{ print "Input omniscient mandatory to use check_level1_positions!"; exit; }
	if( defined($args->{feature})) {$feature_l1 = $args->{feature};} else{ print "Input feature mandatory to use check_level1_positions!"; exit; }
	if( defined($args->{verbose}) ) { $verbose = $args->{verbose}; } else { $verbose = 0;}
	if( defined($args->{log}) ) { $log = $args->{log}; }


	my $extrem_start=1000000000000;
	my $extrem_end=0;
	my $check_existence_feature_l2=undef;
	my $id_l1 = lc($feature_l1->_tag_value('ID'));
	my $tag_l1 = lc($feature_l1->primary_tag());

	# Skip top and standalone features
	if (! exists_keys ($hash_omniscient, ('other', 'level', 'level1') ) ){ # Check info is present in $hash_omniscient
		get_levels_info({verbose => 0, omniscient => $hash_omniscient});
	}
	if ($hash_omniscient->{'other'}{'level'}{'level1'}{$tag_l1} eq 'standalone' or
				$hash_omniscient->{'other'}{'level'}{'level1'}{$tag_l1} eq 'topfeature'){
		return $result;
	}

	foreach my $tag_level2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
    	if ( exists_keys ($hash_omniscient, ('level2', $tag_level2, $id_l1) ) ){
    		$check_existence_feature_l2=1;

	    	my $extrem_start_A=1000000000000;
		  	my $extrem_end_A=0;
	   		foreach my $feature ( @{$hash_omniscient->{'level2'}{$tag_level2}{$id_l1}}) {
          if( $feature->seq_id eq $feature_l1->seq_id ){
	      		my $start=$feature->start();
	      		my $end=$feature->end();
	      		if ($start < $extrem_start_A){
	       			$extrem_start_A=$start;
	      		}
	      		if($end > $extrem_end_A){
	        		$extrem_end_A=$end;
	      		}
          }
	      }

	    	if ($extrem_start_A < $extrem_start){
	    		$extrem_start=$extrem_start_A;
	    	}
	    	if($extrem_end_A > $extrem_end){
	    		$extrem_end=$extrem_end_A;
	    	}
	    }
    }
    if(! $check_existence_feature_l2){
    	dual_print($log, "check_level1_positions: NO level2 feature to check positions of the level1 feature !".$feature_l1->gff_string()."\n", $verbose);
    }
    else{
	    # modify START if needed
	    if($feature_l1->start != $extrem_start){
	    	$feature_l1->start($extrem_start);
	    	$result=1;
	    	dual_print($log, "check_level1_positions: We modified the L1 LEFT extremity for the sanity the biological data!\n", 0); # print in log only
	    }

	    # modify END if needed
	    if($feature_l1->end != $extrem_end){
	    	$feature_l1->end($extrem_end);
	    	$result=1;
	    	dual_print($log, "check_level1_positions: We modified the L1 RIGHT extremity for the sanity the biological data!\n", 0); # print in log only
	    }
	}
	return $result;
}

# Check the start and end of level2 feature based on all features level3;
sub check_level2_positions {
	my ($hash_omniscient, $level2_feature)=@_;

	 my @values = $level2_feature->get_tag_values('ID');
     my $level2_feature_name = lc(shift @values) ;

	my $extrem_start=1000000000000;
	my $extrem_end=0;
	foreach my $tag_level3 (keys %{$hash_omniscient->{'level3'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
    	my $extrem_start_A=1000000000000;
		my $extrem_end_A=0;
		if( exists_keys ($hash_omniscient, ('level3', $tag_level3, $level2_feature_name) ) ){
	  		foreach my $feature ( @{$hash_omniscient->{'level3'}{$tag_level3}{$level2_feature_name}}) {
	   			my $start=$feature->start();
	   			my $end=$feature->end();
	   			if ($start < $extrem_start_A){
	   				$extrem_start_A=$start;
	   			}
	   			if($end > $extrem_end_A){
	      			$extrem_end_A=$end;
	   			}
	   		}
	   	}
   		if ($extrem_start_A < $extrem_start){
    		$extrem_start=$extrem_start_A;
    	}
   		if($extrem_end_A > $extrem_end){
    		$extrem_end=$extrem_end_A;
    	}
    }

    # modify START if needed
    if($level2_feature->start != $extrem_start){
    	$level2_feature->start($extrem_start);
    }

    # modify END if needed
    if($level2_feature->end != $extrem_end){
    	$level2_feature->end($extrem_end);
    }
}

# return 1 if feature overlaps
sub check_features_overlap {
  my ($feature1, $feature2)=@_;
  my $result = undef;

  if( ($feature2->start <= $feature1->end) and ($feature2->end >= $feature1->start ) ){ # they overlap
    $result = 1;
  }
  return $result;
}

#calcul the overlaping percentage betwwen 2 CDS list or 2 exon list etc...
# /!\ Be careful if you test the output, a overlaping gene can have a percentage overlap to 0. And if you test "if(featuresList_overlap)" and you have a 0, it will fail. So you have to check defined(featuresList_overlap)
sub featuresList_overlap {

	my ($listCDS1ref, $listCDS2ref)=@_;
	my $resu;

	###
	# sort the list
	my @listCDS1 = sort {$a->start <=> $b->start} @{$listCDS1ref};
	my @listCDS2 = sort {$a->start <=> $b->start} @{$listCDS2ref};
	# foreach my $t (@listCDS1){
	# 	print "list1: ".$t->start." ".$t->end."\n";
	# }
	# foreach my $t (@listCDS2){
	# 	print "list2: ".$t->start." ".$t->end."\n";
	# }
	my $size_overlap=0;
	my $cds1_size=0;
	foreach my $cds1 (@listCDS1){

		$cds1_size=$cds1_size+($cds1->end - $cds1->start)+1;
		my $starto;
		my $endo;

		foreach my $cds2 (@listCDS2){
			if($cds2->start > $cds1->end){ #we are after of the investigated cds.
				last;
			}
			elsif($cds2->end < $cds1->start){ # we are before investigated cds.
				next;
			}
			else{ #we are overlaping
				#check start
				if($cds1->start >= $cds2->start){
					$starto=$cds1->start;
				}
				else{$starto=$cds2->start;}
				#check end
				if($cds1->end >= $cds2->end){
					$endo=$cds2->end;
				}
				else{$endo=$cds1->end;}

				#calcul overlap;
				$size_overlap=$size_overlap+($endo - $starto + 1);
			}
		}
	}

	#Now calcul percentage overlap
	if($size_overlap != 0){
		$resu=($size_overlap*100)/$cds1_size;
		$resu = sprintf('%.0f', $resu);
		return $resu;
	}
	else{return undef;}
}


sub featuresList_identik {
	my ($list1, $list2)=@_;

	my @slist1 = sort {$a->start <=> $b->start} @{$list1};
	my @slist2 = sort {$a->start <=> $b->start} @{$list2};
	my $identik="true";

	if($#slist1 == $#slist2){
		my $cpt=0;

		while ($cpt <= $#slist1){

			my $feature1=$slist1[$cpt];
			my $feature2=$slist2[$cpt];
			#print $feature1->start." != ".$feature2->start." or  ".$feature1->end." != ".$feature2->end." or  ".$feature1->strand." ne ".$feature2->strand." or  ".$feature1->seq_id." ne ".$feature2->seq_id."\n";
			if( ($feature1->start != $feature2->start) or ($feature1->end != $feature2->end) or ($feature1->strand ne $feature2->strand) or ($feature1->seq_id ne $feature2->seq_id)){
				$identik=undef;last;
			}
			$cpt++;
		}
	}
	else{$identik=undef;}
	return $identik;
}

# @Purpose: Check the start and end of l1 l2 features from l3 features.
# @input: 2 => hash(omniscient hash), string(gene identifier)
# @output: none
sub check_record_positions {
	my ($hash_omniscient, $gene_id_raw)=@_;

	my $gene_id = lc($gene_id_raw);
	my $ExtremStart=1000000000000;
	my $ExtremEnd=0;

  	foreach my $primary_tag_l2 (keys %{$hash_omniscient->{'level2'}} ){
	  	if (exists_keys($hash_omniscient, ('level2', $primary_tag_l2, $gene_id ) ) ){
		    foreach my $mrna_feature ( @{$hash_omniscient->{'level2'}{$primary_tag_l2}{$gene_id}} ) {
		      	my $l2_id = lc($mrna_feature->_tag_value('ID'));
		      	my $l2_ExtremStart=1000000000000;
	  			my $l2_ExtremEnd=0;
		      	foreach my $tag_l3 ( keys %{$hash_omniscient->{'level3'}} ) {
	  				if ( exists_keys ( $hash_omniscient, ('level3', $tag_l3, $l2_id ) ) ){
		    			foreach my $feature_l3 ( @{$hash_omniscient->{'level3'}{$tag_l3}{$l2_id}} ) {

						    if ($feature_l3->start() < $l2_ExtremStart){
						       $l2_ExtremStart = $feature_l3->start();
						    }
						    if($feature_l3->end() > $l2_ExtremEnd){
						       $l2_ExtremEnd = $feature_l3->end();
			      			}
			      		}
		      		}
		      	}
		      	if ($mrna_feature->start != $l2_ExtremStart and $l2_ExtremStart != 1000000000000){
			      $mrna_feature->start($l2_ExtremStart);
			   	}
			  	if($mrna_feature->end != $l2_ExtremEnd and $l2_ExtremEnd != 0){
			   	 $mrna_feature->end($l2_ExtremEnd);
			  	}
			  	if ( $l2_ExtremStart < $ExtremStart ){
			  		$ExtremStart = $l2_ExtremStart;
			  	}
			  	if ($l2_ExtremEnd > $ExtremEnd  ){
			  		$ExtremEnd = $l2_ExtremEnd;
			  	}
		    }
		}


	  	my $gene_feature=$hash_omniscient->{'level1'}{'gene'}{$gene_id};
	  	if ($gene_feature->start != $ExtremStart and $ExtremStart != 1000000000000){
	      $gene_feature->start($ExtremStart);
	   	}
	  	if($gene_feature->end != $ExtremEnd and $ExtremEnd != 0){
	   	 $gene_feature->end($ExtremEnd);
	  	}
	}
}

# @Purpose: Sort by locusID and location
# @input: 2 => hash(omniscient hash), optional list or hash of id to filter
# @output 3: return 3 similar hashes:  LocusID->uniqLocationId = [id => X, tag => Y].
# One contains all features L1(wihtout topfeature), the second only standalone features, the third only top features
# The standalone feature hash is needed when we loop over omniscient from L2. We will miss stand alone featurs that
# only have L1 features.  It is way to still access them.

sub collect_l1_info_sorted_by_seqid_and_location{
	my ($omniscient, $filterid) = @_;

	my %hash_sortBySeq; # all classified feature L1 even standalone feature
	my %hash_topfeatures; # topfeatures
	my %hash_stdfeatures; # standalone features only

	#Check option filterid
	my $hash_filterid;
	if ($filterid){
		if( ref($filterid) eq 'ARRAY' ){
			$hash_filterid->{$_}++ for (@{$filterid});
		}
		elsif( ref($filterid) eq 'HASH' ){
			$hash_filterid = $filterid;
		}
		else{
			warn "optional filterid parameter need to be a list or hash ref\n";
		}
	}

	# get list of feature type that shoud appear at the very beginning of each sequence id
	my $top_features = get_feature_type_by_agat_value($omniscient, 'level1', 'topfeature');
	my $std_features = get_feature_type_by_agat_value($omniscient, 'level1', 'standalone');

	foreach my $tag_level1 ( keys %{$omniscient->{'level1'}}){
		foreach my $level1_id ( keys %{$omniscient->{'level1'}{$tag_level1}}){

				if ($hash_filterid and ! exists_keys ($hash_filterid, $level1_id) ){
					next;
				}

				my $seqid = $omniscient->{'level1'}{$tag_level1}{$level1_id}->seq_id;
				my $uniq = $omniscient->{'level1'}{$tag_level1}{$level1_id}->start."|".$omniscient->{'level1'}{$tag_level1}{$level1_id}->end.$tag_level1.$level1_id;

				if ( exists_keys($top_features, ($tag_level1) ) ) {
					$hash_topfeatures{$seqid}{$uniq} = { tag => $tag_level1, id => $level1_id };
				}
				else{
					if ( exists_keys($std_features, ($tag_level1) ) ) { #save in standalone is standalone
						$hash_stdfeatures{$seqid}{$uniq} = { tag => $tag_level1, id => $level1_id };
					}
					$hash_sortBySeq{$seqid}{$uniq} = { tag => $tag_level1, id => $level1_id };
				}
	  }
	}
	return \%hash_sortBySeq, \%hash_stdfeatures, \%hash_topfeatures;
}

# Sort by locusID
# LocusID->typeFeature = [feature, feature, feature]
#return a hash. Key is position,tag and value is list of feature l1. The list is sorted
sub gather_and_sort_l1_by_seq_id{
	my ($omniscient) = @_;

	my %hash_sortBySeq;
	foreach my $tag_level1 ( keys %{$omniscient->{'level1'}}){
		foreach my $level1_id ( keys %{$omniscient->{'level1'}{$tag_level1}}){
	    	my $position=$omniscient->{'level1'}{$tag_level1}{$level1_id}->seq_id;
	    	push (@{$hash_sortBySeq{$position}{$tag_level1}}, $omniscient->{'level1'}{$tag_level1}{$level1_id});
	  	}
	  	foreach my $position_l1 (keys %hash_sortBySeq){
        @{$hash_sortBySeq{$position_l1}{$tag_level1}} = sort { ncmp ($a->start.$a->end.$a->_tag_value('ID'), $b->start.$b->end.$b->_tag_value('ID') ) } @{$hash_sortBySeq{$position_l1}{$tag_level1}};
      }
	}
 	return \%hash_sortBySeq;
}


# Sort by locusID only for l2 of a specific type (tag from 3rd column)
# LocusID->typeFeature = [feature, feature, feature]
#return a hash. Key is position,tag and value is list of feature l1. The list is sorted
sub gather_and_sort_l1_by_seq_id_for_l2type{
	my ($omniscient, $l2_tag) = @_;

	my %hash_sortBySeq;
	foreach my $tag_level1 ( keys %{$omniscient->{'level1'}}){

		foreach my $level1_id ( keys %{$omniscient->{'level1'}{$tag_level1}}){
			if (exists_keys ($omniscient, ('level2', $l2_tag, $level1_id) ) ){
		    my $position=$omniscient->{'level1'}{$tag_level1}{$level1_id}->seq_id;
		    push (@{$hash_sortBySeq{$position}{$tag_level1}}, $omniscient->{'level1'}{$tag_level1}{$level1_id});
			}
	  }
	  foreach my $position_l1 (keys %hash_sortBySeq){
      @{$hash_sortBySeq{$position_l1}{$tag_level1}} = sort { ncmp ($a->start.$a->end.$a->_tag_value('ID'), $b->start.$b->end.$b->_tag_value('ID') ) } @{$hash_sortBySeq{$position_l1}{$tag_level1}};
    }
	}
 	return \%hash_sortBySeq;
}

# Sort by locusID only for l1 of a specific type (tag from 3rd column)
# LocusID->typeFeature = [feature, feature, feature]
#return a hash. Key is position,tag and value is list of feature l1. The list is sorted
sub gather_and_sort_l1_by_seq_id_for_l1type{
	my ($omniscient, $tag_level1) = @_;

	my %hash_sortBySeq;
	if (exists_keys ($omniscient, ('level1', $tag_level1) ) ){

		foreach my $level1_id ( keys %{$omniscient->{'level1'}{$tag_level1}}){
		  my $position=$omniscient->{'level1'}{$tag_level1}{$level1_id}->seq_id;
		  push (@{$hash_sortBySeq{$position}{$tag_level1}}, $omniscient->{'level1'}{$tag_level1}{$level1_id});
		}
	  foreach my $position_l1 (keys %hash_sortBySeq){
      @{$hash_sortBySeq{$position_l1}{$tag_level1}} = sort { ncmp ($a->start.$a->end.$a->_tag_value('ID'), $b->start.$b->end.$b->_tag_value('ID') ) } @{$hash_sortBySeq{$position_l1}{$tag_level1}};
    }
	}

 	return \%hash_sortBySeq;
}

# Sort by locusID and strand
# LocusID_strand->typeFeature = [feature, feature, feature]
# return a hash. Key is position,tag and value is list of feature l1. The list is sorted
sub gather_and_sort_l1_by_seq_id_and_strand{
  my ($omniscient) = @_;

  my %hash_sortBySeq;
    foreach my $tag_level1 (keys %{$omniscient->{'level1'}}){
      	foreach my $level1_id (keys %{$omniscient->{'level1'}{$tag_level1}}){
        	my $level1_feature = $omniscient->{'level1'}{$tag_level1}{$level1_id};
        	my $position_l1=$level1_feature->seq_id.$level1_feature->strand;
        	push (@{$hash_sortBySeq{$position_l1}{$tag_level1}}, $level1_feature);
        }
        foreach my $position_l1 (keys %hash_sortBySeq){
        	@{$hash_sortBySeq{$position_l1}{$tag_level1}} = sort { ncmp ($a->start.$a->end.$a->_tag_value('ID'), $b->start.$b->end.$b->_tag_value('ID') ) } @{$hash_sortBySeq{$position_l1}{$tag_level1}};
        }
    }
  return \%hash_sortBySeq;
}

# @Purpose: Create a hash of level1 location (location = [level1ID,start,end]) sorted by feature type and localisation. A localisation is the sequence_id appended by the strand
# @input: 1 => hash omniscient
# @output: 1 => hash => LocusID->typeFeature =[ID,start,end]
sub gather_and_sort_l1_location_by_seq_id_and_strand{
	my ($omniscient) = @_;

	my %hash_sortBySeq;

  	foreach my $tag_level1 (keys %{$omniscient->{'level1'}}){
    	foreach my $level1_id (keys %{$omniscient->{'level1'}{$tag_level1}}){
	    	my $level1_feature = $omniscient->{'level1'}{$tag_level1}{$level1_id};
	    	my $ID = $level1_feature->_tag_value('ID');
	    	my $strand="+";
	    	if($level1_feature->strand != 1){$strand = "-";}
	    	my $position_l1=$level1_feature->seq_id."".$strand;
	    	push ( @{$hash_sortBySeq{$position_l1}{$tag_level1}}, [$ID, int($level1_feature->start), int($level1_feature->end)] );
        }

        foreach my $position_l1 (keys %hash_sortBySeq){
        	@{$hash_sortBySeq{$position_l1}{$tag_level1}} = sort { ncmp ( $a->[1], $b->[1] ) } @{$hash_sortBySeq{$position_l1}{$tag_level1}};
        }
	}
	return \%hash_sortBySeq;
}

# @Purpose: Create a hash of level1 location (location = [level1ID,start,end]) sorted by feature type and localisation. A localisation is the sequence_id
# @input: 1 => hash omniscient
# @output: 1 => hash => LocusID->typeFeature =[ID,start,end]
sub gather_and_sort_l1_location_by_seq_id{
	my ($omniscient) = @_;

	my %hash_sortBySeq;

  	foreach my $tag_level1 (keys %{$omniscient->{'level1'}}){
    	foreach my $level1_id (keys %{$omniscient->{'level1'}{$tag_level1}}){
	    	my $level1_feature = $omniscient->{'level1'}{$tag_level1}{$level1_id};
	    	my $ID = $level1_feature->_tag_value('ID');
	    	my $position_l1=$level1_feature->seq_id;
	    	push ( @{$hash_sortBySeq{$position_l1}{$tag_level1}}, [$ID, int($level1_feature->start), int($level1_feature->end)] );
        }

        foreach my $position_l1 (keys %hash_sortBySeq){
        	@{$hash_sortBySeq{$position_l1}{$tag_level1}} = sort { ncmp ( $a->[1], $b->[1] ) } @{$hash_sortBySeq{$position_l1}{$tag_level1}};
        }
	}
	return \%hash_sortBySeq;
}

# @Purpose: get position of the most left and right cds positions
# @input: 2 => hash(omniscient hash), [l1 feature /or/ l1 id]
# @output: (integer, integer)
sub get_most_right_left_cds_positions {
	my ($omniscient, $l1_feature) = @_;

	my $cds_start = undef;
	my $cds_end = undef;
 	my $gene_id=undef;
	if (ref($l1_feature) =~ "::"){
	    $gene_id = lc($l1_feature->_tag_value('ID'));
    }
    else{
		$gene_id=lc($l1_feature);
    }

    foreach my $tag_l1 (keys %{$omniscient->{'level1'}}){
		if (exists_keys ($omniscient, ('level1', $tag_l1, $gene_id) ) ){

			# == LEVEL 2 == #
			foreach my $tag_l2 (keys %{$omniscient->{'level2'}}){
				if (exists_keys ($omniscient, ('level2', $tag_l2, $gene_id) ) ){
					foreach my $feature_l2 ( @{$omniscient->{'level2'}{$tag_l2}{$gene_id}}) {
						# == LEVEL 3 == #
						my $l2_id = lc($feature_l2->_tag_value('ID') );
						if (exists_keys ($omniscient, ('level3', 'cds', $l2_id ) ) ){

							my @sorted_cds = sort {$a->start <=> $b->start} @{$omniscient->{'level3'}{'cds'}{$l2_id}};
    						my $local_cds_start  = $omniscient->{'level3'}{'cds'}{$l2_id}[0]->start; #first element of the array
							my $local_cds_end = $omniscient->{'level3'}{'cds'}{$l2_id}[$#{$omniscient->{'level3'}{'cds'}{$l2_id}}]->end; #last element of the array

							if ( ! $cds_start){
								$cds_start = $local_cds_start;
							}
							elsif( $local_cds_start < $cds_start){
								$cds_start = $local_cds_start;
							}
							if ( ! $cds_end){
								$cds_end = $local_cds_end;
							}
							elsif( $local_cds_end > $cds_end){
								$cds_end = $local_cds_end;
							}
						}
					}
				}
			}
		}
	}
	return $cds_start, $cds_end;
}

#				   +------------------------------------------------------+
#				   |+----------------------------------------------------+|
#				   || 			MANIPULATION from record 				 ||
#				   |+----------------------------------------------------+|
#				   +------------------------------------------------------+

# @Purpose: Check the ID of a record and fix it if duplicated.
# @input: 3 => id from level1 hash1, hash1 omnicient, hash1 of whole_IDs
# @output: l1_id
# /!\ Not tested yest
sub check_record_ids {
	my ($l1_id1, $hash_omniscient, $hash_whole_IDs, $hash_whole_IDs2, $miscCount)=@_;

	my $l1_id_final = $l1_id1;
	#################
	# == LEVEL 1 == #
	#################
	my $id_l1 = lc($l1_id1);
	foreach my $tag_l1 (keys %{$hash_omniscient->{'level1'}}){ # tag_l1 = gene or repeat etc...
		if (exists_keys ($hash_omniscient, ('level1', $tag_l1, $id_l1) ) ){

			my $new_parent=undef;
			my $uID = $hash_omniscient->{'level1'}{$tag_l1}{$id_l1}->_tag_value('ID');

			if ( exists ( $hash_whole_IDs->{$id_l1} ) ){
				my $feature = $hash_omniscient->{'level1'}{$tag_l1}{$id_l1};
				$uID = replace_by_uniq_ID( $feature, $hash_whole_IDs2, $hash_whole_IDs, $miscCount);
				$hash_omniscient->{'level1'}{$tag_l1}{lc($uID)} = delete $hash_omniscient->{'level1'}{$tag_l1}{$id_l1}; # save feature level1
				$new_parent=1;
				$l1_id_final = $uID ;
			}
			else{
				$hash_whole_IDs2->{lc($uID)}=$uID;
			}
			#################
			# == LEVEL 2 == #
			#################
			foreach my $tag_l2 (keys %{$hash_omniscient->{'level2'}}){ # tag_l2 = mrna or mirna or ncrna or trna etc...

				if (exists_keys ($hash_omniscient, ('level2', $tag_l2, $id_l1) ) ){ #Non present in hash2, we create a list with one element

					foreach my $feature_l2 ( @{$hash_omniscient->{'level2'}{$tag_l2}{$id_l1}}) {

						my $new_parent_l2=undef;
						if($new_parent){
							create_or_replace_tag($feature_l2, 'Parent', $uID);
						}

						my $id_l2 = $feature_l2->_tag_value('ID');
						my $uID_l2 = undef;
						if ( exists ( $hash_whole_IDs->{lc($id_l2)} ) ){
							$uID_l2 = replace_by_uniq_ID($feature_l2, $hash_whole_IDs2, $hash_whole_IDs, $miscCount);
							$new_parent_l2=1;
						}
						else{$hash_whole_IDs2->{lc($id_l2)} = $id_l2;}

						#################
						# == LEVEL 3 == #
						#################
						foreach my $tag_l3 (keys %{$hash_omniscient->{'level3'}}){

							if (exists_keys ($hash_omniscient, ('level3', $tag_l3, lc($id_l2) ) ) ){

								foreach my $feature_l3 ( @{$hash_omniscient->{'level3'}{$tag_l3}{lc($id_l2)}}) {

									if($new_parent_l2){
										create_or_replace_tag($feature_l3, 'Parent', $uID_l2);
									}

									my $id_l3 = $feature_l3->_tag_value('ID');
									my $uID_l3 = undef;
									if ( exists ( $hash_whole_IDs->{lc($id_l3)} ) ){
										$uID_l3 = replace_by_uniq_ID($feature_l3, $hash_whole_IDs2, $hash_whole_IDs, $miscCount);
									}
									else{$hash_whole_IDs2->{lc($id_l3)} = $id_l3;}
								}
								if($new_parent_l2){
									$hash_omniscient->{'level3'}{$tag_l3}{lc($uID_l2)} = delete $hash_omniscient->{'level3'}{$tag_l3}{lc($id_l2)}; # save feature level1
								}
							}
						}
					}
					if($new_parent){
							$hash_omniscient->{'level2'}{$tag_l2}{lc($uID)} = delete $hash_omniscient->{'level2'}{$tag_l2}{$id_l1}; # save feature level1
					}
				}
			}
		}
	}
	return $l1_id_final;
}

#				   +------------------------------------------------------+
#				   |+----------------------------------------------------+|
#				   || 			MANIPULATION on feature	 				 ||
#				   |+----------------------------------------------------+|
#				   +------------------------------------------------------+

# @Purpose: Check if has a l3 subfeature type.
# @input: 2 => hash(omniscient hash), [l1 feature /or/ l1 id]
# @output: bolean
sub l1_has_l3_type {
	my ($omniscient, $l1_feature, $type, $part_match) = @_;

	my $full_match=1;
	if($part_match){
		$full_match=undef;
	}


	my $gene_id=undef;
	if (ref($l1_feature) =~ "::"){
	    $gene_id = lc($l1_feature->_tag_value('ID'));
    }
    else{
		$gene_id=lc($l1_feature);
    }

    foreach my $tag_l1 (keys %{$omniscient->{'level1'}}){
		if (exists_keys ($omniscient, ('level1', $tag_l1, $gene_id) ) ){

			# == LEVEL 2 == #
			foreach my $tag_l2 (keys %{$omniscient->{'level2'}}){

				if (exists_keys ($omniscient, ('level2', $tag_l2, $gene_id) ) ){
					foreach my $feature_l2 ( @{$omniscient->{'level2'}{$tag_l2}{$gene_id}}) {

						# == LEVEL 3 == #
						if($full_match){
							if (exists_keys ($omniscient, ('level3', $type, lc($feature_l2->_tag_value('ID') ) ) ) ){
								return 1
							}
						}
						else{
							my $level2_ID = lc($feature_l2->_tag_value('ID'));
							foreach my $ptag_l3 (keys %{$omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...
								if (lc($ptag_l3) =~ lc($type)){
									if( exists_keys($omniscient, ('level3', $ptag_l3, $level2_ID)) ){
										return 1;
									}
								}
							}
						}
					}
				}
			}
		}
	}
	return 0;
}

# @Purpose: Check if has a cds l3 subfeature.
# @input: 2 => hash(omniscient hash), l2 feature
# @output: bolean
sub l2_has_cds {
	my ($omniscient, $l2_feature) = @_;

	my $gene_id = lc($l2_feature->_tag_value('Parent'));

	# == LEVEL 2 == #
	foreach my $tag_l2 (keys %{$omniscient->{'level2'}}){
		if (exists_keys ($omniscient, ('level2', $tag_l2, $gene_id) ) ){
			foreach my $feature_l2 ( @{$omniscient->{'level2'}{$tag_l2}{$gene_id}}) {
				# == LEVEL 3 == #
				if (exists_keys ($omniscient, ('level3', 'cds', lc($feature_l2->_tag_value('ID') ) ) ) ){
					return 1
				}
			}
		}
	}
	return 0;
}

# @Purpose: Check if has a cds l3 subfeature.
# @input: 2 => hash(omniscient hash), l2 feature
# @output: 1 => undef or Array ref of CDS features
sub get_cds_from_l2 {
	my ($omniscient, $l2_feature) = @_;

	my $gene_id = lc($l2_feature->_tag_value('Parent'));

	# == LEVEL 2 == #
	foreach my $tag_l2 (keys %{$omniscient->{'level2'}}){
		if (exists_keys ($omniscient, ('level2', $tag_l2, $gene_id) ) ){
			foreach my $feature_l2 ( @{$omniscient->{'level2'}{$tag_l2}{$gene_id}}) {
				# == LEVEL 3 == #
				my $l2_id = lc($feature_l2->_tag_value('ID') );
				if (exists_keys ($omniscient, ('level3', 'cds',  $l2_id) ) ){
					my @sorted_cds = sort {$a->start <=> $b->start} @{$omniscient->{'level3'}{'cds'}{$l2_id}};
					return \@sorted_cds;
				}
			}
		}
	}
	return undef;
}

# @Purpose: check if two features overlap.
# @input: 2 => l1 feature, l2 feature
# @output: bolean
sub check_if_feature_overlap{
	my($feature1, $feature2)=@_;
	my $result=undef;
	if (($feature1->start <= $feature2->end) and ($feature1->end >= $feature2->start)){
		$result="true";
	}

return $result
}

#				   +------------------------------------------------------+
#				   |+----------------------------------------------------+|
#				   || 			MANIPULATION on feature	List 			 ||
#				   |+----------------------------------------------------+|
#				   +------------------------------------------------------+


#				   +------------------------------------------------------+
#				   |+----------------------------------------------------+|
#				   || 			MANIPULATION from id 					 ||
#				   |+----------------------------------------------------+|
#				   +------------------------------------------------------+

#				   +------------------------------------------------------+
#				   |+----------------------------------------------------+|
#				   || 			Info from id/feature 					 ||
#				   |+----------------------------------------------------+|
#				   +------------------------------------------------------+

# @Purpose: check if it is a single exon gene
# @input: 2 => hash(omniscient hash), [l1 feature /or/ l1 id]
# @output: bolean
sub is_single_exon_gene {
	my ($omniscient, $l1_feature) = @_;

	my $gene_id = undef;
	if (ref($l1_feature) =~ "::"){
	    $gene_id = lc($l1_feature->_tag_value('ID'));
    }
    else{
		$gene_id=lc($l1_feature);
    }

    foreach my $tag_l1 (keys %{$omniscient->{'level1'}}){
		if (exists_keys ($omniscient, ('level1', $tag_l1, $gene_id) ) ){

			# == LEVEL 2 == #
			foreach my $tag_l2 (keys %{$omniscient->{'level2'}}){
				if (exists_keys ($omniscient, ('level2', $tag_l2, $gene_id) ) ){
					foreach my $feature_l2 ( @{$omniscient->{'level2'}{$tag_l2}{$gene_id}}) {
						# == LEVEL 3 == #
						my $l2_id = lc($feature_l2->_tag_value('ID') );
						if (exists_keys ($omniscient, ('level3', 'exon', $l2_id ) ) ){
  							if (scalar @{$omniscient->{'level3'}{'exon'}{$l2_id}} == 1){
  								return 1;
  							}
  							else{return 0;}
						}
						else{
							warn "WARNING No exon available to check if it is a single exon gene\n";
						}
					}
				}
			}
		}
	}
}

#				   +-------------------------  END -----------------------------+
1;
