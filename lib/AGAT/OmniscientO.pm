#!/usr/bin/perl -w

package AGAT::OmniscientO;

use strict;
use warnings;
use Text::Wrap qw(wrap);
use Sort::Naturally;
use Bio::Tools::GFF;
use URI::Escape;
use AGAT::OmniscientTool;
use AGAT::OmniscientI;
use AGAT::Utilities;
use Exporter;

our @ISA = qw(Exporter);
our @EXPORT = qw(print_ref_list_feature print_omniscient print_omniscient_as_match
print_omniscient_from_level1_id_list webapollo_compliant embl_compliant
convert_omniscient_to_ensembl_style write_top_features);

sub import {
  AGAT::OmniscientO->export_to_level(1, @_); # to be able to load the EXPORT functions when direct call; (normal case)
  AGAT::OmniscientO->export_to_level(2, @_); # to be able to load the EXPORT functions when called from one level up;
}

=head1 SYNOPSIS

  This is the code to output data store in Omniscient.

=head1 DESCRIPTION



=head1 AUTHOR

   Jacques Dainat - jacques.dainat@nbis.se

=cut

#				   +------------------------------------------------------+
#				   |+----------------------------------------------------+|
#				   || 					          Print Methods 					       ||
#				   |+----------------------------------------------------+|
#				   +------------------------------------------------------+

# omniscient is a hash containing a whole gXf file in memory sorted in a specific way (3 levels)
sub print_omniscient{

	my ($hash_omniscient, $gffout) = @_  ;

	my @l3_out_priority = ("tss", "exon", "cds", "tts");
	#uri_decode_omniscient($hash_omniscient);

  # --------- deal with header --------------
  write_headers($hash_omniscient, $gffout);

### OLD FASHION GOING TRHOUGH LEVEL1
	#foreach my $primary_tag_l1 ( sort {$a <=> $b or $a cmp $b} keys %{$hash_omniscient->{'level1'}}){ # primary_tag_l1 = gene or repeat etc...
	#	foreach my $id_tag_key_level1 ( sort { $hash_omniscient->{'level1'}{$primary_tag_l1}{$a}->start <=> $hash_omniscient->{'level1'}{$primary_tag_l1}{$b}->start } keys %{$hash_omniscient->{'level1'}{$primary_tag_l1}} ) { #sort by position

### NEW FASHION GOING TRHOUGH LEVEL1 - Have to first create a hash of seq_id -> level1_feature , then we can go through in alphanumerical order.

	# sort by seq id
	my ( $hash_sortBySeq, $hash_sortBySeq_std, $hash_sortBySeq_topf ) = collect_l1_info_sorted_by_seqid_and_location($hash_omniscient);

	# Read by seqId to sort properly the output by seq ID
	# sort { (($a =~ /(\d+)$/)[0] || 0) <=> (($b =~ /(\d+)$/)[0] || 0) will provide sorting like that: contig contig1 contig2 contig3 contig10 contig11 contig22 contig100 contig101
	foreach my $seqid (sort { (($a =~ /(\d+)$/)[0] || 0) <=> (($b =~ /(\d+)$/)[0] || 0) } keys %{$hash_sortBySeq}){ # loop over all the feature level1

  #################
	# == LEVEL 1 == # IF not in omniscient do that, otherwise we us within. Make a method for it.
	#################
		write_top_features($gffout, $seqid, $hash_sortBySeq_topf, $hash_omniscient);

		foreach my $locationid ( sort { ncmp ($a, $b) } keys %{$hash_sortBySeq->{$seqid} } ){

				my $primary_tag_l1 = $hash_sortBySeq->{$seqid}{$locationid}{'tag'};
				my $id_tag_key_level1 = $hash_sortBySeq->{$seqid}{$locationid}{'id'};
			  $gffout->write_feature($hash_omniscient->{'level1'}{$primary_tag_l1}{$id_tag_key_level1}); # print feature

				#################
				# == LEVEL 2 == #
				#################
				foreach my $primary_tag_l2 (sort {$a cmp $b} keys %{$hash_omniscient->{'level2'}}){ # primary_tag_l2 = mrna or mirna or ncrna or trna etc...

					if ( exists_keys( $hash_omniscient, ('level2', $primary_tag_l2, $id_tag_key_level1) ) ){
						foreach my $feature_level2 ( sort { ncmp ($a->start.$a->end.$a->_tag_value('ID'), $b->start.$b->end.$b->_tag_value('ID') ) } @{$hash_omniscient->{'level2'}{$primary_tag_l2}{$id_tag_key_level1}}) {
							$gffout->write_feature($feature_level2);

							#################
							# == LEVEL 3 == #
							#################
							my $level2_ID = lc($feature_level2->_tag_value('ID'));

              ###########
              # Before tss
              if ( exists_keys($hash_omniscient,('level3','tss',$level2_ID)) ){
                foreach my $feature_level3 ( @{$hash_omniscient->{'level3'}{'tss'}{$level2_ID}}) {
                  #_uri_encode_one_feature($feature_level3);
                  $gffout->write_feature($feature_level3);
                }
              }

							######
							# FIRST EXON
							if ( exists_keys( $hash_omniscient, ('level3', 'exon', $level2_ID) ) ){
								foreach my $feature_level3 ( sort {$a->start <=> $b->start} @{$hash_omniscient->{'level3'}{'exon'}{$level2_ID}}) {
									$gffout->write_feature($feature_level3);
								}
							}
							###########
							# SECOND CDS
							if ( exists_keys( $hash_omniscient, ('level3', 'cds', $level2_ID) ) ){
								foreach my $feature_level3 ( sort {$a->start <=> $b->start} @{$hash_omniscient->{'level3'}{'cds'}{$level2_ID}}) {
									$gffout->write_feature($feature_level3);
								}
							}

              ###########
              # Last tts
              if ( exists_keys($hash_omniscient,('level3','tts',$level2_ID)) ){
                foreach my $feature_level3 ( @{$hash_omniscient->{'level3'}{'tts'}{$level2_ID}}) {
                  #_uri_encode_one_feature($feature_level3);
                  $gffout->write_feature($feature_level3);
                }
              }

							############
							# THEN ALL THE REST
							foreach my $primary_tag_l3 (sort {$a cmp $b} keys %{$hash_omniscient->{'level3'}}){ # primary_tag_l3 = cds or exon or start_codon or utr etc...
								if (! grep { $_ eq $primary_tag_l3 } @l3_out_priority){
									if ( exists_keys( $hash_omniscient, ('level3', $primary_tag_l3, $level2_ID) ) ){
										foreach my $feature_level3 ( sort {$a->start <=> $b->start} @{$hash_omniscient->{'level3'}{$primary_tag_l3}{$level2_ID}}) {
											$gffout->write_feature($feature_level3);
										}
									}
								}
							}
						}
					}
				}
			}

	}

	# --------- deal with fasta seq --------------
	write_fasta($gffout, $hash_omniscient);
}

# omniscient is a hash containing a whole gXf file in memory sorted in a specific way (3 levels)
sub print_omniscient_as_match{

	my ($hash_omniscient, $gffout) = @_  ;

	#uri_decode_omniscient($hash_omniscient);

  # --------- deal with header --------------
  write_headers($hash_omniscient, $gffout);

	# sort by seq id
	my ( $hash_sortBySeq, $hash_sortBySeq_std, $hash_sortBySeq_topf ) = collect_l1_info_sorted_by_seqid_and_location($hash_omniscient);

	#Read by seqId to sort properly the output by seq ID
	foreach my $seqid ( sort { (($a =~ /(\d+)$/)[0] || 0) <=> (($b =~ /(\d+)$/)[0] || 0) } keys %{$hash_sortBySeq}){ # loop over all the feature level1

  #################
  # == LEVEL 1 == #
  #################

    write_top_features($gffout, $seqid, $hash_sortBySeq_topf, $hash_omniscient);

		foreach my $locationid ( sort { ncmp ($a, $b) } keys %{$hash_sortBySeq->{$seqid} } ){

				my $primary_tag_l1 = $hash_sortBySeq->{$seqid}{$locationid}{'tag'};
				my $id_tag_key_level1 = $hash_sortBySeq->{$seqid}{$locationid}{'id'};

				if($primary_tag_l1 =~ "match"){
					$gffout->write_feature($hash_omniscient->{'level1'}{$primary_tag_l1}{$id_tag_key_level1}); # print feature
				}
				#################
				# == LEVEL 2 == #
				#################
				foreach my $primary_tag_l2 (sort {$a cmp $b} keys %{$hash_omniscient->{'level2'}}){ # primary_tag_l2 = mrna or mirna or ncrna or trna etc...

					if ( exists_keys( $hash_omniscient, ('level2', $primary_tag_l2, $id_tag_key_level1) ) ){
						foreach my $feature_level2 ( sort {$a->start <=> $b->start} @{$hash_omniscient->{'level2'}{$primary_tag_l2}{$id_tag_key_level1}}) {

							if($primary_tag_l2 =~ "match"){
								$gffout->write_feature($feature_level2);
							}
							else{
								$feature_level2->primary_tag('match');
								if( $feature_level2->has_tag('Parent')){
									$feature_level2->remove_tag('Parent');
								}

								$gffout->write_feature($feature_level2);

								#################
								# == LEVEL 3 == #
								#################
								my $level2_ID = $feature_level2->_tag_value('ID');

								######
								# EXON
								if ( exists_keys( $hash_omniscient, ('level3', 'exon', lc($level2_ID)) ) ){
									my $current_start=0;
									foreach my $feature_level3 ( sort {$a->start <=> $b->start} @{$hash_omniscient->{'level3'}{'exon'}{lc($level2_ID)}}) {

										$current_start++;
										my $end=($feature_level3->end - $feature_level3->start)+$current_start;
										$feature_level3->primary_tag('match_part');

										if(! $feature_level3->has_tag('Target')){
											my @target=();
											create_or_replace_tag($feature_level3, "Target", "$level2_ID $current_start $end +"); # Target has value has to be a list correctly formated
										}
										$current_start=$end;

										$gffout->write_feature($feature_level3);
									}
								}
							}
						}
					}
				}
			}
		}

	# --------- deal with fasta seq --------------
	write_fasta($gffout, $hash_omniscient);
}

# omniscient is a hash containing a whole gXf file in memory sorted in a specific way (3 levels)
sub print_omniscient_from_level1_id_list {

	my ($hash_omniscient, $level_id_list, $gffout) = @_  ;
  my %topfeature_printed;
  my @l3_out_priority = ("tss", "exon", "cds", "tts");

  #uri_decode_omniscient($hash_omniscient);

  # --------- deal with header --------------
  write_headers($hash_omniscient, $gffout);

  # sort by seq id
	my ( $hash_sortBySeq, $hash_sortBySeq_std, $hash_sortBySeq_topf ) = collect_l1_info_sorted_by_seqid_and_location($hash_omniscient, $level_id_list);

  #Read by seqId to sort properly the output by seq ID
  foreach my $seqid ( sort { (($a =~ /(\d+)$/)[0] || 0) <=> (($b =~ /(\d+)$/)[0] || 0) } keys %{$hash_sortBySeq}){ # loop over all the feature level1

    #################
    # == LEVEL 1 == #
    #################
		write_top_features($gffout, $seqid, $hash_sortBySeq_topf, $hash_omniscient);

		foreach my $locationid ( sort { ncmp ($a, $b) } keys %{$hash_sortBySeq->{$seqid} } ){

			my $primary_tag_l1 = $hash_sortBySeq->{$seqid}{$locationid}{'tag'};
			my $id_tag_key_level1 = $hash_sortBySeq->{$seqid}{$locationid}{'id'};

    	#_uri_encode_one_feature($hash_omniscient->{'level1'}{$primary_tag_l1}{$id_tag_key_level1});

    	$gffout->write_feature($hash_omniscient->{'level1'}{$primary_tag_l1}{$id_tag_key_level1}); # print feature

    	#################
    	# == LEVEL 2 == #
    	#################
    	foreach my $primary_tag_key_level2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...

    		if ( exists ($hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1} ) ){
    			foreach my $feature_level2 ( @{$hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1}}) {

    				#_uri_encode_one_feature($feature_level2);

    				$gffout->write_feature($feature_level2);

    				#################
    				# == LEVEL 3 == #
    				#################
    				my $level2_ID ;
    				if($feature_level2->has_tag('ID')){
    					$level2_ID = lc($feature_level2->_tag_value('ID'));
    				}
    				elsif($feature_level2->has_tag('transcript_id')){
    					$level2_ID = lc( $feature_level2->_tag_value('transcript_id'));
    				}
    				else{
    					warn "Cannot retrieve the parent feature of the following feature: ".gff_string($feature_level2);
    				}

    				###########
    				# Before tss
    				if ( exists_keys($hash_omniscient,('level3','tss',$level2_ID)) ){
    					foreach my $feature_level3 ( @{$hash_omniscient->{'level3'}{'tss'}{$level2_ID}}) {
    	  				#_uri_encode_one_feature($feature_level3);
    		  			$gffout->write_feature($feature_level3);
    					}
    			  }

						######
						# FIRST EXON
						if ( exists_keys($hash_omniscient,('level3','exon',$level2_ID)) ){
							foreach my $feature_level3 ( @{$hash_omniscient->{'level3'}{'exon'}{$level2_ID}}) {
								#_uri_encode_one_feature($feature_level3);
								$gffout->write_feature($feature_level3);
							}
						}
						###########
						# SECOND CDS
						if ( exists_keys($hash_omniscient,('level3','cds',$level2_ID)) ){
							foreach my $feature_level3 ( @{$hash_omniscient->{'level3'}{'cds'}{$level2_ID}}) {
								#_uri_encode_one_feature($feature_level3);
								$gffout->write_feature($feature_level3);
							}
						}
						###########
						# Last tts
						if ( exists_keys($hash_omniscient,('level3','tts',$level2_ID)) ){
							foreach my $feature_level3 ( @{$hash_omniscient->{'level3'}{'tts'}{$level2_ID}}) {
								#_uri_encode_one_feature($feature_level3);
								$gffout->write_feature($feature_level3);
							}
						}
						###########
    				# The rest
    				foreach my $primary_tag_key_level3 ( sort {$a cmp $b} keys %{$hash_omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...
              if (! grep { $_ eq $primary_tag_key_level3 } @l3_out_priority){
                if ( exists ($hash_omniscient->{'level3'}{$primary_tag_key_level3}{$level2_ID} ) ){
    							foreach my $feature_level3 ( @{$hash_omniscient->{'level3'}{$primary_tag_key_level3}{$level2_ID}}) {
    								#_uri_encode_one_feature($feature_level3);
    								$gffout->write_feature($feature_level3);
    							}
    						}
    					}
    				}
    			}
        }
      }
		}
	}

	# --------- deal with fasta seq --------------
	write_fasta($gffout, $hash_omniscient);
}

# @Purpose: Print the fasta seq at the end of the gff
# @input: 2 =>  gff/gtf file output, omniscient
# @output none => none
sub write_fasta {
	my ($gffout, $hash_omniscient) = @_;

	if ( exists ($hash_omniscient->{'other'}{'fasta'} ) ){
		$gffout->_print("##FASTA\n");
		
		my $gffin = $hash_omniscient->{'other'}{'fasta'};
		my @Bio_Seq_objs =  $gffin->get_seqs();

		for my $Bio_Seq_obj (sort { ncmp ($a->display_id, $b->display_id) } @Bio_Seq_objs){
		
			if( $Bio_Seq_obj->desc ){
				$gffout->_print(">".$Bio_Seq_obj->display_id." ".$Bio_Seq_obj->desc."\n");
			}
			else{
				$gffout->_print(">".$Bio_Seq_obj->display_id."\n");
			}
			
			local $Text::Wrap::columns = 80;
			$gffout->_print(wrap('', '', $Bio_Seq_obj->seq));
			$gffout->_print("\n");
		}

		# Close the gff input FH opened by OmniscientI
		$gffin->close();
	}

}

# @Purpose: Print a list of feature simple apporach (not sorting)
# @input: 2 =>  reference list, gff fh
# @output none => none
sub print_ref_list_feature {

	my ($list, $gffout) = @_  ;

	foreach my $feature (@$list) {
		$gffout->write_feature($feature);
	}
}

# @Purpose: Print the headers when first time we access the fh
# @input: 2 =>  ref omniscient, gff fh
# @output none => none
sub write_headers{
  	my ($hash_omniscient, $gffout) = @_  ;

  # If first time we write in this fh
  if ($gffout->{'_first'} ){

    # we write the very top header describing gff version
    $gffout->_print("##gff-version ".$gffout->gff_version()."\n");

    # Now we inject the header catched when parsing input file
    if (exists_keys( $hash_omniscient, ('other', 'header') ) ){
      my $gffXtra=$gffout->{"_filehandle"}; #to add extra lines to gff!!
      foreach my $header_line ( @{$hash_omniscient->{'other'}{'header'} } ) {
        print $gffXtra $header_line;
      }
    }
    # to avoid to write again the header later in the file
    $gffout->{'_first'} = 0;
  }
}

sub write_top_features{

  my ($gffout, $seqid, $hash_sortBySeq_topf, $hash_omniscient ) = @_;

    if ( exists_keys( $hash_sortBySeq_topf, ($seqid) ) ){

      foreach my $locationid ( sort { ncmp ($a, $b) } keys %{$hash_sortBySeq_topf->{$seqid}} ){
				my $tag_l1 = $hash_sortBySeq_topf->{$seqid}{$locationid}{'tag'};
				my $id_l1 = $hash_sortBySeq_topf->{$seqid}{$locationid}{'id'};
				my $feature_l1 = $hash_omniscient->{'level1'}{$tag_l1}{$id_l1};
        $gffout->write_feature($feature_l1); # print feature
      }
    }
}

#				   +------------------------------------------------------+
#				   |+----------------------------------------------------+|
#				   || 					     webapollo compliant 					       ||
#				   |+----------------------------------------------------+|
#				   +------------------------------------------------------+

###############################
# METHOD RELATED TO WEBAPOLLO #
###############################
use constant CWA_skip_feature => { "non_canonical_three_prime_splice_site" => 1 , "non_canonical_five_prime_splice_site" => 2};
#Transform omniscient data to be Webapollo compliant
sub webapollo_compliant {
		my ($hash_omniscient) = @_  ;

	#################
	# == LEVEL 1 == #
	#################
	foreach my $primary_tag_l1 (keys %{$hash_omniscient->{'level1'}}){ # primary_tag_l1 = gene or repeat etc...
		if(exists (CWA_skip_feature->{$primary_tag_l1})){delete $hash_omniscient->{'level1'}{$primary_tag_l1}; next;}
		foreach my $id_l1 (keys %{$hash_omniscient->{'level1'}{$primary_tag_l1}}){
			webapollo_rendering_l1($hash_omniscient->{'level1'}{$primary_tag_l1}{$id_l1});
			#################
			# == LEVEL 2 == #
			#################
			foreach my $primary_tag_l2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_l2 = mrna or mirna or ncrna or trna etc...
				if(exists (CWA_skip_feature->{$primary_tag_l2})){delete $hash_omniscient->{'level2'}{$primary_tag_l2}; next;}
				if ( exists ($hash_omniscient->{'level2'}{$primary_tag_l2}{$id_l1} ) ){
					foreach my $feature_level2 ( @{$hash_omniscient->{'level2'}{$primary_tag_l2}{$id_l1}}) {
						webapollo_rendering_l2($feature_level2);

						#################
						# == LEVEL 3 == #
						#################
						my $level2_ID  = lc($feature_level2->_tag_value('ID'));

						foreach my $primary_tag_l3 (keys %{$hash_omniscient->{'level3'}}){ # primary_tag_l3 = cds or exon or start_codon or utr etc...
							if(exists (CWA_skip_feature->{$primary_tag_l3})){delete $hash_omniscient->{'level3'}{$primary_tag_l3}; next;}
							if ( exists ($hash_omniscient->{'level3'}{$primary_tag_l3}{$level2_ID} ) ){
								foreach my $feature_level3 ( @{$hash_omniscient->{'level3'}{$primary_tag_l3}{$level2_ID}}) {
									webapollo_rendering_l3($feature_level3);
								}
							}
						}
					}
				}
			}
		}
	}
}

#follow webapollo description for a correct visualisation of data
sub webapollo_rendering_l1 {
	my ($feature)=@_;

	my @tags = $feature->get_all_tags();

	## check Name attribute
	my @f = grep /^\Qname\E$/i, @tags; #look at tag that match the whole string to NAME (case insensitive)
	my $name_tag = $f[0];
	if(! $name_tag){
		my $value = $feature->_tag_value('ID');
		$feature->add_tag_value('Name', $value);
	}
}

	my ($feature)=@_;

#follow webapollo description for a correct visualisation of data
sub webapollo_rendering_l2 {
	my ($feature)=@_;

	## check primary tag
	my $primary_tag = lc($feature->primary_tag);

	my %corrections = (
			mrna => 'mRNA',
		);
	if ( exists $corrections{$primary_tag}) {
		$feature->primary_tag( $corrections{$primary_tag});
	}

	my @tags = $feature->get_all_tags();

	## check product/description attribute
	#	if($feature->has_tag('product')){
	my @f = grep /\Qproduct\E/i, @tags;
	my $product_tag = $f[0];
	if($product_tag){
		my @values = $feature->get_tag_values($product_tag);
		$feature->add_tag_value('description', @values);
		$feature->remove_tag($product_tag);
	}
}

#follow webapollo description for a correct visualisation of data
sub webapollo_rendering_l3 {
	my ($feature)=@_;

	## check primary tag
	my $primary_tag = lc($feature->primary_tag);

	my %corrections = (
			cds => 'CDS',
			exon => 'exon',
			three_prime_utr => 'three_prime_UTR',
			five_prime_utr => 'five_prime_UTR',
			utr => 'UTR',
		);
	if ( exists $corrections{$primary_tag}) {
		$feature->primary_tag( $corrections{$primary_tag});
	}

	my @tags = $feature->get_all_tags();

	foreach my $tag (@tags){
		if(lc($tag) ne 'id' and lc($tag) ne 'parent'){
			$feature->remove_tag($tag);
		}
	}
}

#				   +------------------------------------------------------+
#				   |+----------------------------------------------------+|
#				   || 					         EMBL compliant 					       ||
#				   |+----------------------------------------------------+|
#				   +------------------------------------------------------+

#Transform omniscient data to be embl compliant
sub embl_compliant {
		my ($hash_omniscient) = @_  ;

	#################
	# == LEVEL 1 == #
	#################
	foreach my $primary_tag_key_level1 (keys %{$hash_omniscient->{'level1'}}){ # primary_tag_key_level1 = gene or repeat etc...
		foreach my $id_tag_key_level1 (keys %{$hash_omniscient->{'level1'}{$primary_tag_key_level1}}){
			_embl_rendering($hash_omniscient->{'level1'}{$primary_tag_key_level1}{$id_tag_key_level1});

			#################
			# == LEVEL 2 == #
			#################
			foreach my $primary_tag_key_level2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...

				if ( exists ($hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1} ) ){
					foreach my $feature_level2 ( @{$hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1}}) {
						_embl_rendering($feature_level2);

						#################
						# == LEVEL 3 == #
						#################
						my $level2_ID  = lc($feature_level2->_tag_value('ID'));

						foreach my $primary_tag_key_level3 (keys %{$hash_omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...

							if ( exists ($hash_omniscient->{'level3'}{$primary_tag_key_level3}{$level2_ID} ) ){
								foreach my $feature_level3 ( @{$hash_omniscient->{'level3'}{$primary_tag_key_level3}{$level2_ID}}) {
									_embl_rendering($feature_level3);
								}
							}
						}
					}
				}
			}
		}
	}
}

sub _embl_rendering {

	my ($feature)=@_;

	## check primary tag
	my $primary_tag = lc($feature->primary_tag);

	my @feature_list=["assembly_gap","C_region","CDS","centromere","D-loop","D_segment","exon","gap","gene","iDNA","intron","J_segment","LTR","mat_peptide","misc_binding","misc_difference","misc_feature","misc_recomb","misc_RNA","misc_structure","mobile_element","modified_base","mRNA","ncRNA","N_region","old_sequence","operon","oriT","polyA_site","precursor_RNA","prim_transcript","primer_bind","protein_bind","regulatory","repeat_region","rep_origin","rRNA","S_region","sig_peptide","source","stem_loop","STS","telomere","tmRNA","transit_peptide","tRNA","unsure","V_region","V_segment","variation","3'UTR","5'UTR"];

	foreach my $element (@feature_list){
		if(lc($element) =~ /$primary_tag/){
			$feature->$primary_tag = $element;
		}
		else{
			#repeat exception rule
			if( $primary_tag =~ /repeat/ ){
				$feature->$primary_tag = "repeat_region";
			}
			#utr5 exception rule
			elsif($primary_tag =~ /utr/ and ($primary_tag =~ /3/ or $primary_tag =~ /three/) ){
				$feature->$primary_tag = "3'UTR";
			}
			#utr5 exception rule
			elsif($primary_tag =~ /utr/ and ($primary_tag =~ /5/ or $primary_tag =~ /five/) ){
				$feature->$primary_tag = "5'UTR";
			}
			print "WARNING: this primary tag ".$primary_tag." is not recognized among those expected to be EMBL compliant. Please check it or create an exception rule.\n";
		}
	}
}

#				   +------------------------------------------------------+
#				   |+----------------------------------------------------+|
#				   || 					         ENSEMBL compliant 					     ||
#				   |+----------------------------------------------------+|
#				   +------------------------------------------------------+

# @Purpose: Convert the omniscient to be constistent with the ensembl gff format
# @input: 1 =>  omniscient Hash reference
# @output 1 =>  omniscient Hash reference
sub convert_omniscient_to_ensembl_style{
	my ($omniscient) = @_;

	_convert_3th_column($omniscient, "mRNA", "transcript");
	_add_exon_id($omniscient, "ID");
	_add_transcript_id($omniscient, "ID");
	_add_gene_id($omniscient, "ID");
}

# @Purpose: Convert the 3th column (feature type) from old value to new value
# @input: 3 =>  omniscient Hash reference, String = featurey type original, String = new feature type
# @output none =>  The hash itself is modified
sub _convert_3th_column{
	my ($omniscient, $original, $new) = @_;

	foreach my $primary_tag_l1 (keys %{$omniscient->{'level1'}}){ # primary_tag_l1 = gene or repeat etc...
		foreach my $id_l1 (keys %{$omniscient->{'level1'}{$primary_tag_l1}}){

			if( lc($omniscient->{'level1'}{$primary_tag_l1}{$id_l1}->primary_tag) eq lc($original) ){
				$omniscient->{'level1'}{$primary_tag_l1}{$id_l1}->primary_tag($new);
			}

			#################
			# == LEVEL 2 == #
			#################
			foreach my $primary_tag_l2 (keys %{$omniscient->{'level2'}}){ # primary_tag_l2 = mrna or mirna or ncrna or trna etc...
				if ( exists ($omniscient->{'level2'}{$primary_tag_l2}{$id_l1} ) ){
					foreach my $feature_level2 ( @{$omniscient->{'level2'}{$primary_tag_l2}{$id_l1}}) {

						if(lc($feature_level2->primary_tag) eq lc($original) ){
							$feature_level2->primary_tag($new);
						}
						#################
						# == LEVEL 3 == #
						#################
						my $level2_ID  = lc($feature_level2->_tag_value('ID'));

						foreach my $primary_tag_l3 (keys %{$omniscient->{'level3'}}){ # primary_tag_l3 = cds or exon or start_codon or utr etc...
							if ( exists ($omniscient->{'level3'}{$primary_tag_l3}{$level2_ID} ) ){
								foreach my $feature_level3 ( @{$omniscient->{'level3'}{$primary_tag_l3}{$level2_ID}}) {
									if(lc($feature_level3->primary_tag) eq lc($original) ){
										$feature_level3->primary_tag($new);
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

# @Purpose: add exon_id to exon feature if does not exist.
# @input: 3 =>  omniscient Hash reference, String = attribute to use to create exon_id
# @output none =>  The hash itself is modified
sub _add_exon_id{
	my ($omniscient, $original_attribute) = @_;

	foreach my $primary_tag_l1 (keys %{$omniscient->{'level1'}}){ # primary_tag_l1 = gene or repeat etc...
		foreach my $id_l1 (keys %{$omniscient->{'level1'}{$primary_tag_l1}}){
			#################
			# == LEVEL 2 == #
			#################
			foreach my $primary_tag_l2 (keys %{$omniscient->{'level2'}}){ # primary_tag_l2 = mrna or mirna or ncrna or trna etc...
				if ( exists_keys ($omniscient, ('level2', $primary_tag_l2, $id_l1) ) ){
					foreach my $feature_l2 ( @{$omniscient->{'level2'}{$primary_tag_l2}{$id_l1}}) {
						#################
						# == LEVEL 3 == #
						#################
						my $level2_ID  = lc($feature_l2->_tag_value('ID'));

						if ( exists_keys ($omniscient, ('level3', 'exon', $level2_ID) ) ){
							foreach my $feature_l3 ( @{$omniscient->{'level3'}{'exon'}{$level2_ID}}) {

								if( $feature_l3->has_tag($original_attribute) and ! $feature_l3->has_tag("exon_id") ){
									create_or_replace_tag($feature_l3, "exon_id", $feature_l3->_tag_value($original_attribute) );
								}
							}
						}
					}
				}
			}
		}
	}
}

# @Purpose: add transcript_id to all feature l2 if does not exist.
# @input: 3 =>  omniscient Hash reference, String = attribute to use to create transcript_id
# @output none =>  The hash itself is modified
sub _add_transcript_id{
	my ($omniscient, $original_attribute) = @_;

	foreach my $primary_tag_l1 (keys %{$omniscient->{'level1'}}){ # primary_tag_l1 = gene or repeat etc...
		foreach my $id_l1 (keys %{$omniscient->{'level1'}{$primary_tag_l1}}){
			#################
			# == LEVEL 2 == #
			#################
			foreach my $primary_tag_l2 (keys %{$omniscient->{'level2'}}){ # primary_tag_l2 = mrna or mirna or ncrna or trna etc...
				if ( exists ($omniscient->{'level2'}{$primary_tag_l2}{$id_l1} ) ){
					foreach my $feature_l2 ( @{$omniscient->{'level2'}{$primary_tag_l2}{$id_l1}}) {
						my $level2_ID  = lc($feature_l2->_tag_value('ID'));

 						if( $feature_l2->has_tag($original_attribute) and ! $feature_l2->has_tag('transcript_id') ){
							create_or_replace_tag($feature_l2, "transcript_id", $feature_l2->_tag_value($original_attribute) );
						}
					}
				}
			}
		}
	}
}

# @Purpose: add gene_id to all feature l1 if does not exist.
# @input: 3 =>  omniscient Hash reference, String = attribute to use to create gene_id
# @output none =>  The hash itself is modified
sub _add_gene_id{
	my ($omniscient, $original_attribute) = @_;

	foreach my $primary_tag_l1 (keys %{$omniscient->{'level1'}}){ # primary_tag_l1 = gene or repeat etc...
		foreach my $id_l1 (keys %{$omniscient->{'level1'}{$primary_tag_l1}}){
			my $feature_l1  = $omniscient->{'level1'}{$primary_tag_l1}{$id_l1};

 			if( $feature_l1->has_tag($original_attribute) and ! $feature_l1->has_tag('gene_id') ){
				create_or_replace_tag($feature_l1, "gene_id", $feature_l1->_tag_value($original_attribute) );
			}
		}
	}
}

#				   +------------------------------------------------------+
#				   |+----------------------------------------------------+|
#				   || 					         Deal with URI 					         ||
#				   |+----------------------------------------------------+|
#				   +------------------------------------------------------+

# check all the attribute to URI encode the values
sub uri_encode_omniscient {

	my ($hash_omniscient) = @_  ;

	#################
	# == LEVEL 1 == #
	#################
	foreach my $primary_tag_key_level1 (keys %{$hash_omniscient->{'level1'}}){ # primary_tag_key_level1 = gene or repeat etc...
		foreach my $id_tag_key_level1 (keys %{$hash_omniscient->{'level1'}{$primary_tag_key_level1}}){
			_uri_encode_one_feature($hash_omniscient->{'level1'}{$primary_tag_key_level1}{$id_tag_key_level1});

			#################
			# == LEVEL 2 == #
			#################
			foreach my $primary_tag_key_level2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...

				if ( exists ($hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1} ) ){
					foreach my $feature_level2 ( @{$hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1}}) {
						_uri_encode_one_feature($feature_level2);

						#################
						# == LEVEL 3 == #
						#################
						my $level2_ID  = lc($feature_level2->_tag_value('ID'));

						foreach my $primary_tag_key_level3 (keys %{$hash_omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...

							if ( exists ($hash_omniscient->{'level3'}{$primary_tag_key_level3}{$level2_ID} ) ){
								foreach my $feature_level3 ( @{$hash_omniscient->{'level3'}{$primary_tag_key_level3}{$level2_ID}}) {
									_uri_encode_one_feature($feature_level3);
								}
							}
						}
					}
				}
			}
		}
	}
}

sub uri_decode_omniscient {

	my ($hash_omniscient) = @_  ;

	#################
	# == LEVEL 1 == #
	#################
	foreach my $primary_tag_key_level1 (keys %{$hash_omniscient->{'level1'}}){ # primary_tag_key_level1 = gene or repeat etc...
		foreach my $id_tag_key_level1 (keys %{$hash_omniscient->{'level1'}{$primary_tag_key_level1}}){
			_uri_decode_one_feature($hash_omniscient->{'level1'}{$primary_tag_key_level1}{$id_tag_key_level1});

			#################
			# == LEVEL 2 == #
			#################
			foreach my $primary_tag_key_level2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...

				if ( exists ($hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1} ) ){
					foreach my $feature_level2 ( @{$hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1}}) {
						_uri_decode_one_feature($feature_level2);

						#################
						# == LEVEL 3 == #
						#################
						my $level2_ID  = lc($feature_level2->_tag_value('ID'));

						foreach my $primary_tag_key_level3 (keys %{$hash_omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...

							if ( exists ($hash_omniscient->{'level3'}{$primary_tag_key_level3}{$level2_ID} ) ){
								foreach my $feature_level3 ( @{$hash_omniscient->{'level3'}{$primary_tag_key_level3}{$level2_ID}}) {
									_uri_decode_one_feature($feature_level3);
								}
							}
						}
					}
				}
			}
		}
	}
}

# reencode in uri all value of all attributes of a feature
sub _uri_encode_one_feature {

	my ($feature)=@_;

	my @list_tag = $feature->get_all_tags;

	foreach my $tag (@list_tag){
		my @values = $feature->get_tag_values($tag);
		$feature->remove_tag($tag);
		foreach my $val (@values){
			my $val_checked = uri_unescape($val);
			my $new_val = uri_escape($val_checked );
			$feature->add_tag_value($tag, $new_val);
		}
	}
}

sub _uri_decode_one_feature {

	my ($feature)=@_;

	my @list_tag = $feature->get_all_tags;

	foreach my $tag (@list_tag){
		my @values = $feature->get_tag_values($tag);
		$feature->remove_tag($tag);
		foreach my $val (@values){
			my $new_val = uri_unescape($val);
			$feature->add_tag_value($tag, $new_val);
		}
	}
}

1;
