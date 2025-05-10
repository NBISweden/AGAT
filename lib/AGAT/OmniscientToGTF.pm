#!/usr/bin/perl -w

package AGAT::OmniscientToGTF;

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use Sort::Naturally;
use AGAT::AGAT;
use AGAT::OmniscientTool;
use AGAT::Utilities;
use AGAT::Levels;
use Exporter;

our @ISA = qw(Exporter);
our @EXPORT = qw(print_omniscient_as_gtf);

# --------------------------- INFO ---------------------------
# Set GTF version definitions
#my @GTF3 = ("gene", "transcript", "exon", "CDS", "Selenocysteine", "start_codon", "stop_codon", "three_prime_utr", "five_prime_utr");
#my @GTF2_5 = ("gene", "transcript", "exon", "CDS", "UTR", "start_codon", "stop_codon", "Selenocysteine");
#my @GTF2_2 = ("CDS", "start_codon", "stop_codon", "5UTR", "3UTR", "inter", "inter_CNS", "intron_CNS", "exon");
#my @GTF2_1 = ("CDS", "start_codon", "stop_codon", "exon", "5UTR", "3UTR");
#my @GTF2 = ("CDS", "start_codon", "stop_codon", "exon");
#my @GTF1 = ("CDS", "start_codon", "stop_codon", "exon", "intron");


# --------------------------- SHARED FUNCTIONS ---------------------------

sub print_omniscient_as_gtf{

	my ($args) = @_;
	# Check we receive a hash as ref
	if(ref($args) ne 'HASH'){ warn "Hash Arguments expected for print_omniscient_as_gtf. Please check the call.\n";exit;	}
	
	# Fill the parameters
	my ($hash_omniscient, $gtf_out);
	if( defined($args->{omniscient})) { $hash_omniscient = $args->{omniscient}; } else{ print "Omniscient parameter mandatory to use print_omniscient_as_gtf!"; exit; }
	if( defined($args->{output})) { $gtf_out = $args->{output}; } else{ print "Output parameter mandatory to use print_omniscient_as_gtf!"; exit; }

	my $gtf_version = $CONFIG->{"gtf_output_version"};
	my $verbose = $CONFIG->{"verbose"};
	my $debug = $CONFIG->{"debug"};
	my $relax = undef;
	if ($gtf_version eq "relax"){
		$relax = 1;
	}
	
	# rebuild gene_id and transcript_id feature;
	my %keep_track_gene_id;
	my %keep_track_transcript_id;
	my $previous_tag_l1 = "previous_gene_id";
	my $previous_tag_l2 = "previous_transcript_id";
	# sort by seq id
	my $hash_sortBySeq = gather_and_sort_l1_by_seq_id($hash_omniscient);

	#################
	# == LEVEL 1 == #
	#################
	foreach my $seqid (sort { (($a =~ /(\d+)$/)[0] || 0) <=> (($b =~ /(\d+)$/)[0] || 0) } keys %{$hash_sortBySeq}){ # loop over all the feature level1

		print "seqid $seqid\n" if $debug;
		foreach my $primary_tag_key_level1 (sort {$a cmp $b} keys %{$hash_omniscient->{'level1'}}){
			foreach my $feature_level1 ( @{$hash_sortBySeq->{$seqid}{$primary_tag_key_level1}} ){
				print "\n\n\nlevel1: ".$feature_level1->gff_string."\n" if $debug;
				my $id_tag_key_level1 = lc($feature_level1->_tag_value('ID'));

	      # Gene ID level1
				my $gene_id=undef;
	      my $gene_id_att=undef;
	      if($feature_level1->has_tag('gene_id')){
	        $gene_id_att=$feature_level1->_tag_value('gene_id');
	      }

	      my $transcript_id=undef;
	      my $level3_gene_id=undef;
	      #################
	      # == LEVEL 2 == #
	      #################
	      foreach my $primary_tag_key_level2 ( sort {$a cmp $b} keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...

	        if ( exists_keys ($hash_omniscient, ('level2', $primary_tag_key_level2, $id_tag_key_level1) ) ){
	          foreach my $feature_level2 ( sort {$a->start <=> $b->start} @{$hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1}}) {
							print "\nlevel2: ".$feature_level2->gff_string."\n" if $debug;


	            # Gene ID level2
	            my $gene_id_mrna_att=undef;
	            if($feature_level2->has_tag('gene_id')){
	              $gene_id_mrna_att=$feature_level2->_tag_value('gene_id');
	            }

	            my $transcript_id_mrna_att=undef;
	            if($feature_level2->has_tag('transcript_id')){
	              $transcript_id_mrna_att=$feature_level2->_tag_value('transcript_id');
	            }

	            # get gff3 feature (ID)
	            my $level2_ID = lc($feature_level2->_tag_value('ID'));

	            my $level3_transcript_id=undef;
	            #################
	            # == LEVEL 3 == #
	            #################

	            ############
	            # Go through one time to check if gene_id and transcript_id are present and save them
	            foreach my $primary_tag_key_level3 ( sort {$a cmp $b} keys %{$hash_omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...
	              if ( exists_keys ($hash_omniscient, ('level3', $primary_tag_key_level3, $level2_ID) ) ){
	                foreach my $feature_level3 ( sort { $a->start <=> $b->start } @{$hash_omniscient->{'level3'}{$primary_tag_key_level3}{$level2_ID}}) {
										print "level3: ".$feature_level3->gff_string."\n" if $debug;
	                  #Get level3 gene_id
	                  if(! $level3_gene_id){
	                    if($feature_level3->has_tag('gene_id')){
	                      $level3_gene_id=$feature_level3->_tag_value('gene_id');
	                    }
	                  }

	                  #Get level3 transcript_id
	                  if(! $level3_transcript_id){
	                    if($feature_level3->has_tag('transcript_id')){
	                      $level3_transcript_id=$feature_level3->_tag_value('transcript_id');
	                    }
	                  }
	                  if($level3_gene_id and $level3_transcript_id){last;}
	                }
	              }
	              if($level3_gene_id and $level3_transcript_id){last;}
	            }

	            #################
	            # CHOOSE the gene_id. We take the first from level1 to level3.
	            if($gene_id_att and ! _does_id_exist("gene_id", $gene_id_att, \%keep_track_gene_id) ){
	              $gene_id=$gene_id_att;
	            }
	            elsif($gene_id_mrna_att and ! _does_id_exist("gene_id", $gene_id_mrna_att, \%keep_track_gene_id) ){
	              $gene_id=$gene_id_mrna_att
	            }
	            elsif($level3_gene_id and ! _does_id_exist("gene_id", $level3_gene_id, \%keep_track_gene_id) ){
	              $gene_id=$level3_gene_id;
	            }
	            else{ # We didn't find any gene_id we will the ID of level1 as gene_id.
	              $gene_id=$feature_level1->_tag_value('ID');
	            }

	            #################
	            # CHOOSE the transcript_id. We take the first from level2 to level3.
	            if($transcript_id_mrna_att and ! _does_id_exist("transcript_id", $transcript_id_mrna_att, \%keep_track_transcript_id) ){
	              $transcript_id=$transcript_id_mrna_att;
	            }
	            elsif($level3_transcript_id and ! _does_id_exist("transcript_id", $level3_transcript_id, \%keep_track_transcript_id) ){
	              $transcript_id=$level3_transcript_id;
	            }
	            else{ # We didn't find any gene_id we will the ID of level2 as transcript_id.
	              $transcript_id=$feature_level2->_tag_value('ID');
	            }

	            ##############
	            # Second pass of level3 features
	            # Add gene_id and transcript_id to level3 feature that don't have this information
	            foreach my $primary_tag_key_level3 (keys %{$hash_omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...
	              if ( exists_keys ($hash_omniscient, ('level3', $primary_tag_key_level3, $level2_ID) ) ){
	                foreach my $feature_level3 ( @{$hash_omniscient->{'level3'}{$primary_tag_key_level3}{$level2_ID}}) {

	                  #Check add gene_id
	                  if(! $feature_level3->has_tag('gene_id')) {
	                    $feature_level3->add_tag_value('gene_id', $gene_id);
	                  }
	                  elsif($feature_level3->_tag_value('gene_id') ne $gene_id) { #gene_id different, we replace it.
	                    warn("Level3 ".$feature_level3->_tag_value('ID').": We replace the gene_id ".$feature_level3->_tag_value('gene_id')." by ".$gene_id.". We save original gene_id into $previous_tag_l1 attribute.\n");
											create_or_replace_tag($feature_level3, $previous_tag_l1, $feature_level3->_tag_value('gene_id'));
											create_or_replace_tag($feature_level3, 'gene_id', $gene_id);
	                  }
	                  #Check add transcript_id
	                  if(! $feature_level3->has_tag('transcript_id')){
	                    $feature_level3->add_tag_value('transcript_id', $transcript_id);
	                  }
	                  elsif($feature_level3->_tag_value('transcript_id') ne $transcript_id){ #transcript_id different, we replace it.
	                    warn("Level3 ".$feature_level3->_tag_value('ID').": We replace the transcript_id ".$feature_level3->_tag_value('transcript_id')." by ".$transcript_id.". We save original transcript_id into $previous_tag_l2 attribute.\n");
										  create_or_replace_tag($feature_level3, $previous_tag_l2, $feature_level3->_tag_value('transcript_id'));
											create_or_replace_tag($feature_level3, 'transcript_id', $transcript_id);
	                  }
	                }
	              }
	            }

	            ## add level2 missing information gene_id
	            if(! $feature_level2->has_tag('gene_id')) {
	               $feature_level2->add_tag_value('gene_id', $gene_id);
	            }
	            elsif($feature_level2->_tag_value('gene_id') ne $gene_id) { #gene_id different, we replace it.
	              warn("Level2 ".$feature_level2->_tag_value('ID').": We replace the gene_id ".$feature_level2->_tag_value('gene_id')." by ".$gene_id.". We save original gene_id into $previous_tag_l1 attribute.\n");
								create_or_replace_tag($feature_level2, $previous_tag_l1, $feature_level2->_tag_value('gene_id'));
								create_or_replace_tag($feature_level2, 'gene_id', $gene_id);
	            }
	            # add level2 missing information transcript_id
	            if(! $feature_level2->has_tag('transcript_id')){
	              $feature_level2->add_tag_value('transcript_id', $transcript_id);
	            }
	            elsif($feature_level2->_tag_value('transcript_id') ne $transcript_id){ #gene_id transcript_id, we replace it.
	              warn("Level2 ".$feature_level2->_tag_value('ID').": We replace the transcript_id ".$feature_level2->_tag_value('transcript_id')." by ".$transcript_id.". We save original transcript_id into $previous_tag_l2 attribute.\n");
								create_or_replace_tag($feature_level2, $previous_tag_l2, $feature_level2->_tag_value('transcript_id'));
	              create_or_replace_tag($feature_level2, 'transcript_id', $transcript_id);
	            }
	          }
	        }
	      }

	      ## add level1 missing information gene_id
	      if(! $feature_level1->has_tag('gene_id')) {
	        $gene_id = $feature_level1->_tag_value('ID') if (! $gene_id);
	        $feature_level1->add_tag_value('gene_id', $gene_id);
	      }
	      elsif($feature_level1->_tag_value('gene_id') ne $gene_id) { #gene_id different, we replace it.
	        warn("Level1 ".$feature_level1->_tag_value('ID').": We replace the gene_id ".$feature_level1->_tag_value('gene_id')." by ".$gene_id.". We save original gene_id into $previous_tag_l1 attribute.\n");
					create_or_replace_tag($feature_level1, $previous_tag_l1, $feature_level1->_tag_value('gene_id'));
					create_or_replace_tag($feature_level1,'gene_id', $gene_id);
	      }

				# Save used ID
				$keep_track_gene_id{$gene_id}++;
				$keep_track_transcript_id{$transcript_id}++ if ($transcript_id);
	    }
	  }
	}

	if (! $relax){
		# convert correct
		_convert_feature_type($hash_omniscient, $gtf_version, $verbose);

	}

	# print results
	_print_omniscient_filter($hash_omniscient, $gtf_version, $relax, $gtf_out);

}

# --------------------------- LOCAL FUNCTIONS ---------------------------

# ---------------------------------
# check if the ID has already been used for another record. In that case we will not use it but we need to inform the user.
sub _does_id_exist{
	my ($id_type, $id, $hash_to_check) = @_;
	my $result=0;

	if( exists_keys($hash_to_check,(lc($id))) ){
		warn "$id_type $id has already been used earlier, I will use another uniq id value.\n";
		$result = 1;
	}

	return $result;
}

# convert feature type to correct one expected.
# All l1 will become gene type excepted for topfeature and standalone features
# that will be discarded.
sub _convert_feature_type{
	my ($hash_omniscient, $gtf_version, $verbose)=@_;

	my $topfeatures = get_feature_type_by_agat_value($hash_omniscient, 'level1', 'topfeature');
	my $standalones = get_feature_type_by_agat_value($hash_omniscient, 'level1', 'standalone');

	# all l1 are gene now
	foreach my $tag_l1 ( keys %{$hash_omniscient->{'level1'}}){
		if(exists_keys($topfeatures,($tag_l1))){ print "throw $tag_l1\n" if $verbose; next; }
		if(exists_keys($standalones,($tag_l1))){ print "throw $tag_l1\n" if $verbose; next; }

		foreach my $id_l1 ( keys %{$hash_omniscient->{'level1'}{$tag_l1}}){
			if (lc($tag_l1) ne "gene"){
				print "convert $tag_l1 to gene feature\n" if $verbose;
				my $feature = $hash_omniscient->{'level1'}{$tag_l1}{$id_l1};
				create_or_replace_tag( $feature, "original_biotype", $tag_l1 );
				$feature->primary_tag("gene");
			}

			# all l2 are transcript now
			foreach my $tag_l2 ( keys %{$hash_omniscient->{'level2'}}){
				if ( exists_keys ($hash_omniscient, ('level2',$tag_l2,$id_l1) ) ){
					foreach my $feature_level2 ( @{$hash_omniscient->{'level2'}{$tag_l2}{$id_l1}}) {
						if (lc($tag_l2) ne "transcript"){

							create_or_replace_tag($feature_level2, "original_biotype", $tag_l2 );
							$feature_level2->primary_tag("transcript");
						}

						# We do not check for GTF1 and GTF2 UTR specifically because we will remove them later
						if($gtf_version == "1" or $gtf_version == "2"){next;}

						# get gff3 feature (ID)
						my $level2_ID = lc($feature_level2->_tag_value('ID'));


						foreach my $tag_l3_lc (keys %{$hash_omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...

							#shortened CDS to remove stop codon if any
							if ( $tag_l3_lc =~ "cds"){
								if (exists_keys($hash_omniscient, ('level3', 'stop_codon', $level2_ID) ) ){
									my @list_cds = sort {$a->start <=> $b->start} @{$hash_omniscient->{'level3'}{'cds'}{$level2_ID}};
									my $first_cds =  $list_cds[0];
									my $first_cds_size = $first_cds->end - $first_cds->start + 1;
									my $last_cds = $list_cds[ $#list_cds ];
									my $last_cds_size = $last_cds->end - $last_cds->start + 1;
									my $strand = $first_cds->strand;

									my @list_stop = sort {$a->start <=> $b->start} @{$hash_omniscient->{'level3'}{'stop_codon'}{$level2_ID}};

									if ( ($strand eq "+" ) or ($strand eq "1" ) ) {
										my $stop_codon = $list_stop[$#list_stop];
										my $stop_codon_size = ($stop_codon->end - $stop_codon->start + 1);

										if (check_features_overlap($stop_codon, $last_cds) ){
											if( $stop_codon_size == 3){
												if($last_cds_size > 3){
													$last_cds->end($stop_codon->start-1);
												}
												else{_remove_cds($hash_omniscient, $level2_ID, $last_cds);}
											}
											elsif( $stop_codon_size < 3){
												if($last_cds_size > $stop_codon_size){
													$last_cds->end($stop_codon->start);
												}
												else{_remove_cds($hash_omniscient, $level2_ID, $last_cds);}
												# there is another stop_codon
												if($#list_stop > 0){
													$stop_codon = $list_stop[$#list_stop-1];
													$last_cds = $list_cds[ $#list_cds - 1];
													$last_cds->end($stop_codon->start - 1);
												}
											}
											else{
												print "Warning: stop codon longer than 3 nucleotides for $level2_ID\n";
											}
										}
									}
									else{ # minus strand
										my $stop_codon = $list_stop[0];
										my $stop_codon_size = ($stop_codon->end - $stop_codon->start + 1);

										if (check_features_overlap($stop_codon, $first_cds) ){
											if( $stop_codon_size == 3){
												if($first_cds_size > 3){
													$first_cds->start($stop_codon->end+1);
												}
												else{_remove_cds($hash_omniscient, $level2_ID, $first_cds);}
											}
											elsif( $stop_codon_size < 3){
												if($first_cds_size > $stop_codon_size){
													$last_cds->start($stop_codon->end);
												}
												else{_remove_cds($hash_omniscient, $level2_ID, $first_cds);}
												# there is another stop_codon
												if($#list_stop > 0){
													$stop_codon = $list_stop[1];
													$first_cds =  $list_cds[1];
													$first_cds->start($stop_codon->end+1);
												}
											}
											else{
												print "Warning: stop codon longer than 3 nucleotides for $level2_ID\n";
											}
										}
									}
								}
							}

							if ( $tag_l3_lc =~ "utr"){
								if ( exists_keys ($hash_omniscient, ('level3', $tag_l3_lc, $level2_ID) ) ){
									foreach my $feature_level3 ( @{$hash_omniscient->{'level3'}{$tag_l3_lc}{$level2_ID}}) {
										my $tag_l3 = $feature_level3->primary_tag;
										if ($gtf_version eq "2.5"){
											if ($tag_l3 ne "UTR"){
												create_or_replace_tag($feature_level3, "original_biotype", $tag_l3 );
												$feature_level3->primary_tag("UTR");
											}
										}
										else{
											if ( ($tag_l3 =~ "5") or ( $tag_l3_lc =~ "five") ){
												if ($gtf_version eq "3"){
													if ($tag_l3 ne "five_prime_utr"){
														create_or_replace_tag($feature_level3, "original_biotype", $tag_l3 );
														$feature_level3->primary_tag("five_prime_utr");
													}
												}
												else{
													if ($tag_l3 ne "5UTR"){
														create_or_replace_tag($feature_level3, "original_biotype", $tag_l3 );
														$feature_level3->primary_tag("5UTR");
													}
												}
											}
											if ( ($tag_l3 =~ "3") or ( $tag_l3_lc =~ "three") ){
												if ($gtf_version eq "3"){
													if ($tag_l3 ne "three_prime_utr"){
														create_or_replace_tag($feature_level3, "original_biotype", $tag_l3 );
														$feature_level3->primary_tag("three_prime_utr");
													}
												}
												else{
													if ($tag_l3 ne "3UTR"){
														create_or_replace_tag($feature_level3, "original_biotype", $tag_l3 );
														$feature_level3->primary_tag("3UTR");
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
	}
}

sub _remove_cds{
	my ($hash_omniscient, $id_l2, $cds_feature)=@_;

	my @new_cds_list=();
	my $cds_ft_string_comparison = _uniq_comparison($cds_feature);

	foreach my $feature ( @{$hash_omniscient->{'level3'}{'cds'}{$id_l2}} ) {
		my $ft_string_comparison = _uniq_comparison($feature);

		if( $ft_string_comparison eq  $cds_ft_string_comparison){
			next;
		}
		push @new_cds_list, $feature;
	}
	@{$hash_omniscient->{'level3'}{'cds'}{$id_l2}} = @new_cds_list;
}

# Make a uniq string id for comprison.
# Needed because ID is not enough because CDS can share identifiers
sub _uniq_comparison{
	my ($feature)=@_;
	return lc($feature->_tag_value('ID')).$feature->start().$feature->end()
}

# filter feature type to remove not expected ones
sub _print_omniscient_filter{
	my ($hash_omniscient, $gtf_version, $relax, $gffout)=@_;

	# --------- deal with header --------------
  	_write_headers_gtf($hash_omniscient, $gffout, $gtf_version, $relax);

	# sort by seq id
	my ( $hash_sortBySeq, $hash_sortBySeq_stdf,  $hash_sortBySeq_topf) = collect_l1_info_sorted_by_seqid_and_location($hash_omniscient);

	# Read by seqId to sort properly the output by seq ID
	# sort { (($a =~ /(\d+)$/)[0] || 0) <=> (($b =~ /(\d+)$/)[0] || 0) will provide sorting like that: contig contig1 contig2 contig3 contig10 contig11 contig22 contig100 contig101
	foreach my $seqid (sort { ncmp ($a, $b) } keys %{$hash_sortBySeq}){ # loop over all the feature level1

		# ----- LEVEL 1 -----
		_write_top_features_gtf($gffout, $seqid, $hash_sortBySeq_topf, $hash_omniscient);

		foreach my $locationid ( sort { ncmp ($a, $b) } keys %{$hash_sortBySeq->{$seqid} } ){
			my $primary_tag_l1 = $hash_sortBySeq->{$seqid}{$locationid}{'tag'};
			my $id_tag_key_level1 = $hash_sortBySeq->{$seqid}{$locationid}{'id'};

			my $feature_l1 = $hash_omniscient->{'level1'}{$primary_tag_l1}{$id_tag_key_level1};
			my $primary_tag_l1_gtf = lc($feature_l1->primary_tag());
			$gffout->write_feature($feature_l1); # print feature

			# ----- LEVEL 2 -----
			foreach my $primary_tag_l2 (sort {$a cmp $b} keys %{$hash_omniscient->{'level2'}}){ # primary_tag_l2 = mrna or mirna or ncrna or trna etc...
				if ( exists_keys( $hash_omniscient, ('level2', $primary_tag_l2, $id_tag_key_level1) ) ){
					foreach my $feature_level2 ( sort { ncmp ($a->start.$a->end.$a->_tag_value('ID'), $b->start.$b->end.$b->_tag_value('ID') ) } @{$hash_omniscient->{'level2'}{$primary_tag_l2}{$id_tag_key_level1}}) {

						my $primary_tag_l2_gtf = lc($feature_level2->primary_tag());
						$gffout->write_feature($feature_level2); # print feature

						# ----- LEVEL 3 -----
						my $level2_ID = lc($feature_level2->_tag_value('ID'));
						my @l3_done;

						######
						# FIRST EXON
						if ( exists_keys( $hash_omniscient, ('level3', 'exon', $level2_ID) ) ){
							foreach my $feature_level3 ( sort {$a->start <=> $b->start} @{$hash_omniscient->{'level3'}{'exon'}{$level2_ID}}) {
								$gffout->write_feature($feature_level3);
								push @l3_done, 'exon';
							}
						}
						###########
						# SECOND CDS
						if ( exists_keys( $hash_omniscient, ('level3', 'cds', $level2_ID) ) ){
							foreach my $feature_level3 ( sort {$a->start <=> $b->start} @{$hash_omniscient->{'level3'}{'cds'}{$level2_ID}}) {
								$gffout->write_feature($feature_level3);
								push @l3_done, 'cds';
							}
						}

						############
						# THEN ALL THE REST
						foreach my $primary_tag_l3 (sort {$a cmp $b} keys %{$hash_omniscient->{'level3'}}){ # primary_tag_l3 = cds or exon or start_codon or utr etc...
							if (! grep { $_ eq $primary_tag_l3 } @l3_done){
								if ( exists_keys( $hash_omniscient, ('level3', $primary_tag_l3, $level2_ID) ) ){
									my $primary_tag_l3_gtf = lc($hash_omniscient->{'level3'}{$primary_tag_l3}{$level2_ID}->[0]->primary_tag() );

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
}

# @Purpose: Print the headers when first time we access the fh
# @input: 2 =>  ref omniscient, gff fh
# @output none => none
sub _write_headers_gtf{
  	my ($hash_omniscient, $gffout, $gtf_version, $relax ) = @_  ;

  # If first time we write in this fh
  if ($gffout->{'_first'} ){

    # we write the very top header describing gff version
		if ($relax){
			$gffout->_print("##gtf-version X\n");
			my $gffXtra=$gffout->{"_filehandle"}; #to add extra lines to gff!!
			print $gffXtra "# GFF-like GTF i.e. not checked against any GTF specification. Conversion based on GFF input, standardised by AGAT.\n";
		}else{
    	$gffout->_print("##gtf-version ".$gtf_version."\n");
		}

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

sub _write_top_features_gtf{

  my ( $gffout, $seqid, $hash_sortBySeq_topf, $hash_omniscient ) = @_;

  if ( exists_keys( $hash_sortBySeq_topf, ($seqid) ) ){

    foreach my $locationid ( sort { ncmp ($a, $b) } keys %{$hash_sortBySeq_topf->{$seqid}} ){
			my $tag_l1 = $hash_sortBySeq_topf->{$seqid}{$locationid}{'tag'};
			my $id_l1 = $hash_sortBySeq_topf->{$seqid}{$locationid}{'id'};
			my $feature_l1 = $hash_omniscient->{'level1'}{$tag_l1}{$id_l1};

			$gffout->write_feature($feature_l1); # print feature
    }
  }
}

1;