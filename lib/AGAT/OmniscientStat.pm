#!/usr/bin/perl -w

package AGAT::OmniscientStat;

use strict;
use warnings;
use Bio::Tools::GFF;
use Bio::SeqIO;
use AGAT::OmniscientTool;
use AGAT::OmniscientJson;
use AGAT::Utilities;
use AGAT::PlotR;
use Try::Tiny;
use Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw( print_omniscient_statistics );

sub import {
  AGAT::OmniscientStat->export_to_level(1, @_); # to be able to load the EXPORT functions when direct call; (normal case)
  AGAT::OmniscientStat->export_to_level(2, @_); # to be able to load the EXPORT functions when called from one level up;
}

=head1 SYNOPSIS

This is the code to perform statisctis of data store in Omniscient.

=head1 DESCRIPTION

	A library to get statistics from an omniscient hash of a gff3 file


	We create a complex hash of hash containing all information needeed.
	The data are scaned from level 2 to level 1 and 3. We do that because different type of feature from level 2 can have same type of feature of level1. (e.g: mRNA => gene and tRNA => gene).
	So the structure of the hash created is the following:
	{type_feature_level2}{'level'}{type_feature_level}{'flag'}='value';
	'level' can be level1, level2 or level accordingly, allow to go all over the data for printing by driving the data form level1 to level3.
	'flag' correspond to the type of information that has been saved in 'value'

=head1 AUTHOR

   Jacques Dainat - jacques.dainat@nbis.se


=cut

sub print_omniscient_statistics{

#	---	HANDLE ARGUMENTS ---
	my ($args) = @_	;

	# Check we receive a hash as ref
	if(ref($args) ne 'HASH'){ print "Hash Arguments expected for print_omniscient_statistics. Please check the call.\n";exit;	}

	# Declare all variables and fill them
	my ($omniscient, $genome_size, $output, $verbose, $distri, $isoform);

	# omniscient
	if( defined($args->{input})) {$omniscient = $args->{input};}
		else{ print "Input omniscient mandatory to use print_omniscient_statistics!"; exit;}

	#genome size
	if( ! defined($args->{genome}) ) {
		$genome_size = undef;
	}
	else{
		my $genome = $args->{genome};
		# --- check genome size ---
		if($genome){
			if( $genome =~ /^[0-9]+$/){ #check if it's a number
				$genome_size = $genome;
			}
			elsif($genome){
				my $seqio = Bio::SeqIO->new(-file => $genome, '-format' => 'Fasta');
				while(my $seq = $seqio->next_seq) {
		  		my $string = $seq->seq;
		  		$genome_size += length($string);
		  	}
			}
			printf("%-45s%d%s", "Total sequence length", $genome_size,"\n");
		}
	}

	# statistics output
	if( ! defined($args->{output}) ) {
		$output = IO::File->new();
		$output->fdopen( fileno(STDOUT), 'w' );
	}
		else{	$output = $args->{output};}
	# add verbosity
	if( ! defined($args->{verbose}) ) {$verbose = 0;}
		else{ $verbose = $args->{verbose}; }
	# Path to the folder where to put distribution plot
	if( ! defined($args->{distri}) ) {$distri = 0;}
		else{ $distri = $args->{distri}; }
	# Should we deal with isoform (remove them and re-compute the statistics)
	if( ! defined($args->{isoform}) ) {$isoform = 0;}
		else{ $isoform = $args->{isoform}; }
	print "get_omniscient_statistics\n" if $verbose;
	my $result_by_type = get_omniscient_statistics($omniscient, $genome_size, $verbose);
	my $omniscientNew = undef ; #if isoform has to be removed
	my $result_by_type2 = undef; #if isoform will be a computed without isoforms

	print $output ("-"x80)."\n\n";

	# --- print statistics ---
	# by_main_type = 1(topfeatures), 2(standalone features), or 3 (L1 features with children)
	foreach my $by_main_type  ( sort {$a <=> $b } keys %{$result_by_type} ){
		my $isoform_type = ($by_main_type eq 3) ? $isoform : undef;
		foreach my $by_type ( sort keys %{ $result_by_type->{$by_main_type} } ){

			my $stat = $result_by_type->{$by_main_type}{$by_type}{'info'};
			my $distri_hash = $result_by_type->{$by_main_type}{$by_type}{'distri'};
			my $l1l2 = $result_by_type->{$by_main_type}{$by_type}{'iso'};

			# print sentences/info
			foreach my $infoList (@$stat){
				if($isoform_type){
					print $output "Compute $by_type with isoforms if any\n\n";
				}
				else{
					print $output "Compute $by_type\n\n";
				}

				foreach my $info (@$infoList){
			    print $output "$info";
			  }
			  print $output "\n";
			}

			if($distri){
				_print_distribution($distri, "with_isoforms", $distri_hash);
			}

			#------- DEAL WITH ISOFORMS -----
			if($isoform_type){
				print $output "Re-compute $by_type without isoforms asked. We remove shortest isoforms if any\n\n";

				if(! $omniscientNew){ # re-compute wihtout isoforms only once!
					# create a new omniscient with only one mRNA isoform per gene
					my ($nb_iso_removed_cds,  $nb_iso_removed_exon) = remove_shortest_isoforms($omniscient);
					$omniscientNew = 1;
				}

				#get stat without isoform
				print "get_omniscient_statistics_from_l2 wihtout iso for $by_type\n" if $verbose;
				# get nb of each feature in omniscient;
				my ($info_l2, $extra_l2) = get_omniscient_statistics_from_l2($omniscient, $by_type, $verbose);
				my $stat2 = get_info_sentences($info_l2, $extra_l2, $genome_size);
				my $distri2 = get_distributions($info_l2, $extra_l2);

				# print sentences/info
				foreach my $infoList2 (@$stat2){
					foreach my $info2 (@$infoList2){
						print $output "$info2";
					}
					print $output "\n";
				}

				if($distri){;
					_print_distribution($distri, "without_isoforms", $distri_hash);
				}
			}
			print $output ("-"x80)."\n\n";
		}
	}
}

# @Purpose: Purpose print distribution from feature statistics
# @input: 2 =>	String (folder),	string (sub folder with or without iso) and hash (distribution)
# @output: 0
sub _print_distribution{
  my ($folder, $subfolder, $distri)=@_;

  foreach my $type (keys %{$distri} ) {

    foreach my $level (keys %{$distri->{$type}} ) {
      foreach my $tag ( keys %{$distri->{$type}{$level}} ) {
        if( exists_keys ($distri,($type, $level, $tag, 'whole') ) ){

          if(! -d $folder){
            mkdir $folder;
          }

          if(! -d $folder."/".$subfolder){
            mkdir $folder."/".$subfolder;
          }

          my $outputPDF = $folder."/".$subfolder."/".$type."Class_".$tag.".pdf";

          #CREATE THE R COMMAND
          my $nbValues = @{$distri->{$type}{$level}{$tag}{'whole'}};
          my $R_command = rcc_plot_from_list($distri->{$type}{$level}{$tag}{'whole'}, "", "histogram", "$tag"." size (nt)", "Number of $tag", "Distribution of $tag sizes\nMade with $nbValues $tag", $outputPDF);
          #EXECUTE THE R COMMAND
          execute_R_command($R_command);
        }

        if( exists_keys ($distri,($type, $level, $tag, 'piece') ) ){
        }
      }
    }
  }
}

# Calculate information necessary going through the omniscient only once
# return a lisf of sub_list - Sub list contain all inforamtion level1,2,3 of all feature linked to a type of feature of level 2.
# (eg: Gene(l1),mRNA(l2),cds(l3),exon(l3), where the type of level1 and level3 feature are only those linked to mRNA.)
sub get_omniscient_statistics {

	my ($hash_omniscient, $genomeSize, $verbose) = @_  ;

	my %result_by_type;

	#my $out = Bio::Tools::GFF->new(-fh => \*STDOUT, -gff_version => 3);

	# --- get features by type ---
	my ( $hash_sortBySeq, $hash_sortBySeq_stdf, $hash_sortBySeq_topf) = collect_l1_info_sorted_by_seqid_and_location($hash_omniscient);

	# --- get statistics from topfeatures -------------------------
	my $topfeatures = get_feature_type_by_agat_value($hash_omniscient, 'level1', 'topfeature');
	foreach my $tag_l1 ( sort keys %{ $topfeatures }){
		if ( exists_keys ($hash_omniscient, ('level1', $tag_l1) ) ){
			print "get_omniscient_statistics_for_topfeature for $tag_l1\n" if $verbose;
			my ($info_l1, $extra_l1) = get_omniscient_statistics_for_topfeature($hash_omniscient, $tag_l1);
			my $info_l1_sentence = get_info_sentences($info_l1, $extra_l1, $genomeSize);
			my $info_l1_distri = get_distributions($info_l1, $extra_l1);
			$result_by_type{1}{$tag_l1} = { info => $info_l1_sentence, distri => $info_l1_distri, iso =>  undef};
		}
	}

	# --- get statistics from standalone features -------------------------
	my $stdfeatures = get_feature_type_by_agat_value($hash_omniscient, 'level1', 'standalone');
	foreach my $tag_l1 ( sort keys %{ $stdfeatures }){
		if ( exists_keys ($hash_omniscient, ('level1', $tag_l1) ) ){
				print "get_omniscient_statistics_for_standalone\n" if $verbose;
				my ($info_l1, $extra_l1) = get_omniscient_statistics_for_topfeature($hash_omniscient, $tag_l1); #normal title is topfeature
				my $info_l1_sentence = get_info_sentences($info_l1, $extra_l1, $genomeSize);
				my $info_l1_distri = get_distributions($info_l1, $extra_l1);
				$result_by_type{2}{$tag_l1} = { info => $info_l1_sentence, distri => $info_l1_distri, iso =>  undef};
		}
	}

	# ------------------------- get statistic from l2 -------------------------
	print "get_omniscient_statistics_from_l2\n" if $verbose;
	# get nb of each feature in omniscient;
	foreach my $tag_l2 ( sort keys %{$hash_omniscient->{'level2'} }){
		print "tag_l2 $tag_l2\n" if $verbose;
		my ($info_l2, $extra_l2) = get_omniscient_statistics_from_l2($hash_omniscient, $tag_l2, $verbose);

		my $info_l2_sentence = get_info_sentences($info_l2, $extra_l2, $genomeSize);
		my $info_l2_distri = get_distributions($info_l2, $extra_l2);


		#chech if isoforms
		my $nbLevel1 = 0;
		my $nbLevel2 = 0;
		foreach my $idl1 ( keys %{ $hash_omniscient->{'level2'}{$tag_l2} } ){
			$nbLevel2 += scalar @{$hash_omniscient->{'level2'}{$tag_l2}{$idl1} };
			$nbLevel1++;
		}

		my $l1l2 = [$nbLevel1,$nbLevel2];

		$result_by_type{3}{$tag_l2} = { info => $info_l2_sentence, distri => $info_l2_distri, iso =>  $l1l2};
	}

	return \%result_by_type;
}

# Get statistics for top features
sub get_omniscient_statistics_for_topfeature{
	my ($omniscient, $tag_l1) = @_;

	my %all_info;
	my %extra_info; #For info not sorted by Level.
	my $sortBySeqL1 = gather_and_sort_l1_by_seq_id_for_l1type($omniscient, $tag_l1);

	foreach my $id_l1 ( sort keys %{$omniscient->{'level1'}{$tag_l1}}){

		my $feature_l1=$omniscient->{'level1'}{$tag_l1}{$id_l1};

		#count number of feature
		$all_info{$tag_l1}{'level1'}{$tag_l1}{'nb_feat'}++;

		#compute feature size
		my $sizeFeature=($feature_l1->end-$feature_l1->start)+1;
		$all_info{$tag_l1}{'level1'}{$tag_l1}{'size_feat'} += $sizeFeature;

		#create distribution list
		push @{$all_info{$tag_l1}{'level1'}{$tag_l1}{'distribution'}}, $sizeFeature;

		# grab longest
		if ((! $all_info{$tag_l1}{'level1'}{$tag_l1}{'longest'}) or ($all_info{$tag_l1}{'level1'}{$tag_l1}{'longest'} < $sizeFeature)){
			$all_info{$tag_l1}{'level1'}{$tag_l1}{'longest'}=$sizeFeature;
		}

		# grab shorter
		if ((! $all_info{$tag_l1}{'level1'}{$tag_l1}{'shortest'}) or ($all_info{$tag_l1}{'level1'}{$tag_l1}{'shortest'} > $sizeFeature)){
			$all_info{$tag_l1}{'level1'}{$tag_l1}{'shortest'}=$sizeFeature;
		}

		# count how many overlaping genes
		my $nb_overlap_gene = _detect_overlap_features($omniscient, $sortBySeqL1, 'level1');
		$extra_info{"overlap"}{$tag_l1}{"level1"}{$tag_l1} = $nb_overlap_gene;
	}
	return \%all_info, \%extra_info;
}

# Parse omiscient by L2 to seprate statistics e.g not mixing exon from mRNA of
# those tRNA
sub get_omniscient_statistics_from_l2{
	my ($hash_omniscient, $tag_l2, $verbose) = @_;

	my %all_info;
	my %extra_info; #For info not sorted by Level.
	my $sortBySeqL2 = gather_and_sort_l1_by_seq_id_for_l2type($hash_omniscient, $tag_l2);
	my %tags_l1; # to hold all tag L1 related to the sutdied L2

	foreach my $id_l1 ( sort keys %{$hash_omniscient->{'level2'}{$tag_l2}}){
	#               +----------------------------------------------------+
	#               |                     FEATURE LEVEL1                 |
	#               +----------------------------------------------------+

		my $feature_l1=undef;
		my $tag_l1;
		# retrieve the l1 tag	on the fly because a same L2 type can be attached to several l1 types
		foreach my $tag_level1 (keys %{$hash_omniscient->{'level1'}}){
			if ( exists_keys ($hash_omniscient, ('level1', $tag_level1, $id_l1) ) ){
				$feature_l1=$hash_omniscient->{'level1'}{$tag_level1}{$id_l1};
				$tag_l1=$tag_level1;
				$tags_l1{$tag_l1}++;
				last;
			}
		}
		if(! $feature_l1){print "Problem ! We didnt retrieve the level1 feature with id $id_l1\n";exit;}

		#count number of feature
		$all_info{$tag_l2}{'level1'}{$tag_l1}{'nb_feat'}++;

		#compute feature size
		my $sizeFeature=($feature_l1->end-$feature_l1->start)+1;
		$all_info{$tag_l2}{'level1'}{$tag_l1}{'size_feat'} += $sizeFeature;

		#create distribution list
		push @{$all_info{$tag_l2}{'level1'}{$tag_l1}{'distribution'}}, $sizeFeature;

		# grab longest
		if ((! $all_info{$tag_l2}{'level1'}{$tag_l1}{'longest'}) or ($all_info{$tag_l2}{'level1'}{$tag_l1}{'longest'} < $sizeFeature)){
			$all_info{$tag_l2}{'level1'}{$tag_l1}{'longest'}=$sizeFeature;
		}

		# grab shorter
		if ((! $all_info{$tag_l2}{'level1'}{$tag_l1}{'shortest'}) or ($all_info{$tag_l2}{'level1'}{$tag_l1}{'shortest'} > $sizeFeature)){
			$all_info{$tag_l2}{'level1'}{$tag_l1}{'shortest'}=$sizeFeature;
		}

	#               +----------------------------------------------------+
	#               |                     FEATURE LEVEL2                 |
	#               +----------------------------------------------------+
		my $counterL2_match=-1;
		my $All_l2_single=1;
		foreach my $feature_l2 ( @{$hash_omniscient->{'level2'}{$tag_l2}{$id_l1}} ){
			#print $feature_l2->gff_string()."\n";
			#count number of feature
			$all_info{$tag_l2}{'level2'}{$tag_l2}{'nb_feat'}++;

			#compute feature size
			my $sizeFeature=($feature_l2->end-$feature_l2->start)+1;
				$all_info{$tag_l2}{'level2'}{$tag_l2}{'size_feat'} += $sizeFeature;

			#create distribution list
			push @{$all_info{$tag_l2}{'level2'}{$tag_l2}{'distribution'}}, $sizeFeature;

				# grab longest
				if ((! $all_info{$tag_l2}{'level2'}{$tag_l2}{'longest'}) or ($all_info{$tag_l2}{'level2'}{$tag_l2}{'longest'} < $sizeFeature)){
				$all_info{$tag_l2}{'level2'}{$tag_l2}{'longest'}=$sizeFeature;
			}
			# grab shorter
			if ((! $all_info{$tag_l2}{'level2'}{$tag_l2}{'shortest'}) or ($all_info{$tag_l2}{'level2'}{$tag_l2}{'shortest'} > $sizeFeature)){
				$all_info{$tag_l2}{'level2'}{$tag_l2}{'shortest'}=$sizeFeature;
			}

			########################################################
			# Special case match match_part => calcul the introns
			########################################################
			if($tag_l2 =~ "match"){
				my @sortedList = sort {$a->start <=> $b->start} @{$hash_omniscient->{'level2'}{$tag_l2}{$id_l1}};
				#if(! exists ($all_info{$tag_l2}{'level2'}{'intron'}{'nb_feat'}))	{$all_info{$tag_l2}{'level2'}{'intron'}{'nb_feat'}=0;}
				#if(! exists ($all_info{$tag_l2}{'level2'}{'intron'}{'size_feat'}))	{$all_info{$tag_l2}{'level2'}{'intron'}{'size_feat'}=0;}
				my $indexLastL2 = $#{$hash_omniscient->{'level2'}{$tag_l2}{$id_l1}};
				$counterL2_match++;

				if($counterL2_match > 0 and $counterL2_match <= $indexLastL2){
	  			my $intronSize = $sortedList[$counterL2_match]->start - $sortedList[$counterL2_match-1]->end - 1;

	  			#compute feature size
	  			$all_info{$tag_l2}{'level2'}{'intron'}{'size_feat'} += $intronSize;

	  			#create distribution list
					push @{$all_info{$tag_l2}{'level2'}{'intron'}{'distribution'}}, $sizeFeature;

					# grab longest
	    		if ((! $all_info{$tag_l2}{'level2'}{'intron'}{'longest'}) or ($all_info{$tag_l2}{'level2'}{'intron'}{'longest'} < $intronSize)){
						$all_info{$tag_l2}{'level2'}{'intron'}{'longest'}=$intronSize;
					}
					# grab shorter
	    		if ((! $all_info{$tag_l2}{'level2'}{'intron'}{'shortest'}) or ($all_info{$tag_l2}{'level2'}{'intron'}{'shortest'} > $intronSize)){
	    			$all_info{$tag_l2}{'level2'}{'intron'}{'shortest'}=$intronSize;
	    		}
	  			#Count number
	    		$all_info{$tag_l2}{'level2'}{'intron'}{'nb_feat'} += 1;
				}
			}

	#               +----------------------------------------------------+
	#               |                     FEATURE LEVEL3                 |
	#               f+----------------------------------------------------+
			my $utr3 = undef;
			my $utr5 = undef;
			my $id_l2=lc($feature_l2->_tag_value('ID'));
	  	foreach my $tag_l3 (keys %{$hash_omniscient->{'level3'}}){

	  		if( exists_keys ($hash_omniscient, ('level3', $tag_l3, $id_l2) ) ){
					my $sizeMultiFeat=0;
					my $counterL3=-1;
					my $indexLast = $#{$hash_omniscient->{'level3'}{$tag_l3}{$id_l2}};

					my @sortedList = sort {$a->start <=> $b->start} @{$hash_omniscient->{'level3'}{$tag_l3}{$id_l2}};
	  			foreach my $feature_l3 ( @sortedList ){

	  				#count number feature of tag_l3 type
	  				$counterL3++;

	  				#-------------------------------------------------
	  				#				Manage Introns
	  				#-------------------------------------------------
	  				# from the second intron to the last (from index 1 to last index of the table sortedList)
	  				# We go inside this loop only if we have more than 1 feature.
	  				if($counterL3 > 0 and $counterL3 <= $indexLast){
	  					my $intronSize = $sortedList[$counterL3]->start - $sortedList[$counterL3-1]->end - 1;

	  					#compute feature size
	  					$all_info{$tag_l2}{'level3'}{$tag_l3}{'intron'}{'size_feat'} += $intronSize;

	  					#create distribution list
							push @{$all_info{$tag_l2}{'level3'}{$tag_l3}{'intron'}{'distribution'}}, $sizeFeature;

	  					# grab longest
	    	  		if ((! $all_info{$tag_l2}{'level3'}{$tag_l3}{'intron'}{'longest'}) or ($all_info{$tag_l2}{'level3'}{$tag_l3}{'intron'}{'longest'} < $intronSize)){
								$all_info{$tag_l2}{'level3'}{$tag_l3}{'intron'}{'longest'}=$intronSize;
							}

							# grab shorter
	    				if ((! $all_info{$tag_l2}{'level3'}{$tag_l3}{'intron'}{'shortest'}) or ($all_info{$tag_l2}{'level3'}{$tag_l3}{'intron'}{'shortest'} > $intronSize)){
	    					$all_info{$tag_l2}{'level3'}{$tag_l3}{'intron'}{'shortest'}=$intronSize;
	    				}

	  					#Count number
	  					$all_info{$tag_l2}{'level3'}{$tag_l3}{'intron'}{'nb_feat'} += 1;
	  				}

	  				#compute cumulative feature size
	  				my $sizeFeature=($feature_l3->end-$feature_l3->start)+1;
	  				$all_info{$tag_l2}{'level3'}{$tag_l3}{'size_feat'} += $sizeFeature;

	  				#-------------------------------------------------
	  				# MANAGE SPREAD FEATURES (multi exon features)
	  				#-------------------------------------------------
	  	  		if(($tag_l3 =~ /cds/) or ($tag_l3 =~ /utr/)){
	  	  			$sizeMultiFeat += $sizeFeature;
	  	  			$all_info{$tag_l2}{'level3'}{$tag_l3}{'exon'}{'nb_feat'}++;

	  	  			#### MANAGE piece of multi exon features (spread features)

	  	  			#create distribution list of multifeature piece
							push @{$all_info{$tag_l2}{'level3'}{$tag_l3}{'piece'}{'distribution'}}, $sizeFeature;

							# grab longest
	    	  		if ((! $all_info{$tag_l2}{'level3'}{$tag_l3}{'piece'}{'longest'}) or ($all_info{$tag_l2}{'level3'}{$tag_l3}{'piece'}{'longest'} < $sizeFeature)){
								$all_info{$tag_l2}{'level3'}{$tag_l3}{'piece'}{'longest'}=$sizeFeature;
							}
							# grab shorter
	    				if ((! $all_info{$tag_l2}{'level3'}{$tag_l3}{'piece'}{'shortest'}) or ($all_info{$tag_l2}{'level3'}{$tag_l3}{'piece'}{'shortest'} > $sizeFeature)){
	    					$all_info{$tag_l2}{'level3'}{$tag_l3}{'piece'}{'shortest'}=$sizeFeature;
	    				}
	  	  		}
	  	  		#-------------------------------------------------
	  				# MANAGE single FEATURES (multi exon features)
	  				#-------------------------------------------------
	  	  		else{
	  	  			#count number of feature
	  					$all_info{$tag_l2}{'level3'}{$tag_l3}{'nb_feat'}++;

	  	  			#create distribution list of multifeature piece
							push @{$all_info{$tag_l2}{'level3'}{$tag_l3}{'distribution'}}, $sizeFeature;

	    	  		# grab longest
	    	  		if ((! $all_info{$tag_l2}{'level3'}{$tag_l3}{'longest'}) or ($all_info{$tag_l2}{'level3'}{$tag_l3}{'longest'} < $sizeFeature)){
								$all_info{$tag_l2}{'level3'}{$tag_l3}{'longest'}=$sizeFeature;
							}
							# grab shorter
	    				if ((! $all_info{$tag_l2}{'level3'}{$tag_l3}{'shortest'}) or ($all_info{$tag_l2}{'level3'}{$tag_l3}{'shortest'} > $sizeFeature)){
	    					$all_info{$tag_l2}{'level3'}{$tag_l3}{'shortest'}=$sizeFeature;
	    				}
	    			}
	    			####################
	  				#mange utr per mRNA
	  				if ($tag_l3 =~ /three_prime_utr/){
							$utr3=1;
	  				}
	  				if ($tag_l3 =~ /five_prime_utr/){
	  					$utr5=1;
	  				}
	  			}# END FOREACH L3

	  			#----------------------------------------
	  			# NOW TAKE CARE OF MULTIFEATURE AND L2
	  			#in that case the feature was split in several peaces that have been glue together
	  			if (($tag_l3 =~ /utr/) or ($tag_l3 =~ /cds/)){
	  				#count number of feature
	  				$all_info{$tag_l2}{'level3'}{$tag_l3}{'nb_feat'}++;

	  				#create distribution list
						push @{$all_info{$tag_l2}{'level3'}{$tag_l3}{'distribution'}}, $sizeMultiFeat;

	  				# grab longest
	    	  	if ((! $all_info{$tag_l2}{'level3'}{$tag_l3}{'longest'}) or ($all_info{$tag_l2}{'level3'}{$tag_l3}{'longest'} < $sizeMultiFeat)){
							$all_info{$tag_l2}{'level3'}{$tag_l3}{'longest'}=$sizeMultiFeat;
						}

						# grab shorter
	    			if ((! $all_info{$tag_l2}{'level3'}{$tag_l3}{'shortest'}) or ($all_info{$tag_l2}{'level3'}{$tag_l3}{'shortest'} > $sizeMultiFeat)){
	    				$all_info{$tag_l2}{'level3'}{$tag_l3}{'shortest'}=$sizeMultiFeat;
	    			}
	  			}

	  			if ($tag_l3 =~ /exon/){
	  				if ($indexLast == 0) {
	    				# body...
	    				$extra_info{'single'}{$tag_l2}{'level2'}{$tag_l2}++;
	    			}
	    			else{
	    				$All_l2_single=undef;
	    			}
	  			}
				}
			}# END all feature level 3

	  	# 1) Manage UTR both side
	  	if ($utr3  and $utr5){
	    	$all_info{$tag_l2}{'level2'}{$tag_l2}{'utr_both_side'}++;
	    	$all_info{$tag_l2}{'level2'}{$tag_l2}{'utr_at_least_one_side'}++;
	  	} # 2) Manage UTR at least one side
	  	elsif ($utr3  or $utr5){
	 			$all_info{$tag_l2}{'level2'}{$tag_l2}{'utr_at_least_one_side'}++;
			}
		}# END all feature level 2

		if($All_l2_single){
			$extra_info{'single'}{$tag_l2}{'level1'}{$tag_l1}++;
		}
	}
	# count how many overlaping genes
	my $nb_overlap_gene = _detect_overlap_features($hash_omniscient, $sortBySeqL2);
	#One l2 type can be linked to several l1 type. Create a dedicated name
	my $name;
	foreach my $oneTag (keys %tags_l1){
		if (!$name){$name=$oneTag;}
		else{$name.="..".$oneTag;}
	}
	if (!$name){$name="gene";}
	$extra_info{"overlap"}{$tag_l2}{"level1"}{$name} = $nb_overlap_gene;

	return \%all_info, \%extra_info;
}

sub get_info_sentences{
	my ($all_info, $extra_info, $genomeSize) = @_;

	my @result_list;

	# create the list of sentences that resume the results
	foreach my $type (sort keys %{$all_info}){

		my $hashType = $all_info->{$type};
		my @result;

		my $info_number = _info_number($hashType);
		push @result, @$info_number;

		foreach my $extra_key (sort keys %{$extra_info} ){

			if( exists_keys ( $extra_info, ($extra_key, $type) ) ){

				if ($extra_key eq "single" ){
					my $info_single = _info_single( $extra_info->{'single'}{$type});
					push @result, @$info_single;
				}
				if ($extra_key eq "overlap" ){
					my $info_single = _info_overlap( $extra_info->{'overlap'}{$type});
					push @result, @$info_single;
				}
			}
		}

		my $info_mean_per = _info_mean_per($hashType);
		push @result, @$info_mean_per;

		my $info_length = _info_length($hashType);
		push @result, @$info_length;

		my $info_mean_length = _info_mean_length($hashType);
		push @result, @$info_mean_length;

		if($genomeSize){
			my $info_coverage = _info_coverage($hashType, $genomeSize);
	 		push @result, @$info_coverage;
		}

		my $info_longest = _info_longest($hashType);
		push @result, @$info_longest;

		my $info_shortest = _info_shortest($hashType);
		push @result, @$info_shortest;

		push @result_list, \@result;
	}
	return \@result_list;
}


sub get_distributions{
	my ($all_info, $extra_info) = @_;

	my %distribution;

	#extract distribution values
	foreach my $type (sort keys %{$all_info} ) {
		foreach my $level (keys %{$all_info->{$type}} ) {

			foreach my $tag ( keys %{$all_info->{$type}{$level}} ) {

				if( exists_keys ( $all_info, ($type, $level, $tag, 'distribution') ) ){
					$distribution{$type}{$level}{$tag}{'whole'} =  delete $all_info->{$type}{$level}{$tag}{'distribution'};
				}
				if( exists_keys ( $all_info, ($type, $level, $tag, 'piece', 'distribution') ) ){
					$distribution{$type}{$level}{$tag}{'piece'} = delete $all_info->{$type}{$level}{$tag}{'piece'}{'distribution'};
				}
			}
		}
	}

return \%distribution;
}

# ----------------------------------  -----------------------------------------

#####
# Give info about single exon gene and mRNA
sub _info_single{

	my ($all_info) = @_  ;
	my @resu;

	#print level1
	foreach my $tag_l1 (sort keys %{$all_info->{'level1'}}){
		push @resu, sprintf("%-45s%d%s", "Number of single exon $tag_l1", $all_info->{'level1'}{$tag_l1},"\n");
	}

	#print level2
	foreach my $tag_l2 (sort keys %{$all_info->{'level2'}}){
	    push @resu, sprintf("%-45s%d%s", "Number of single exon $tag_l2", $all_info->{'level2'}{$tag_l2},"\n");
	}

	return \@resu;
}

#####
# Give info about snumber of overlaping genes
sub _info_overlap{

	my ($all_info) = @_  ;
	my @resu;

	#print level1
	foreach my $tag_l1 (sort keys %{$all_info->{'level1'}}){
		push @resu, sprintf("%-45s%d%s", "Number $tag_l1 overlapping", $all_info->{'level1'}{$tag_l1},"\n");
	}

	return \@resu;
}

#####
# Give info about number of feature of each type
sub _info_number {

	my ($all_info) = @_  ;
	my @resu;
	my $there_is_utr=undef;

	#print level1
	foreach my $tag_l1 (sort keys %{$all_info->{'level1'}}){
		push @resu, sprintf("%-45s%d%s", "Number of $tag_l1", $all_info->{'level1'}{$tag_l1}{'nb_feat'},"\n");
	}

	#print level2
	foreach my $tag_l2 (sort keys %{$all_info->{'level2'}}){
	    push @resu, sprintf("%-45s%d%s", "Number of $tag_l2", $all_info->{'level2'}{$tag_l2}{'nb_feat'},"\n");
	    #manage utr both side
	    if(exists ($all_info->{'level2'}{$tag_l2}{'utr_both_side'})){
			push @resu, sprintf("%-45s%d%s", "Number of mrnas with utr both sides", $all_info->{'level2'}{$tag_l2}{'utr_both_side'},"\n");
		}
		#manage utr both side
		if(exists ($all_info->{'level2'}{$tag_l2}{'utr_at_least_one_side'})){
			push @resu, sprintf("%-45s%d%s", "Number of mrnas with at least one utr", $all_info->{'level2'}{$tag_l2}{'utr_at_least_one_side'},"\n");
		}
	}

	#print level3 -
	foreach my $tag_l3 (sort keys %{$all_info->{'level3'}}){
	    push @resu, sprintf("%-45s%d%s", "Number of $tag_l3", $all_info->{'level3'}{$tag_l3}{'nb_feat'},"\n");
	}

 	#print level3 - exon case
	foreach my $tag_l3 (sort keys %{$all_info->{'level3'}}){
		if( exists_keys ($all_info, ('level3', $tag_l3, 'exon') ) ) {
	    	push @resu, sprintf("%-45s%d%s", "Number of exon in $tag_l3", $all_info->{'level3'}{$tag_l3}{'exon'}{'nb_feat'},"\n");
	    }
	}
	#print level3 - intron case
	foreach my $tag_l3 (sort keys %{$all_info->{'level3'}}){
		if( exists_keys ($all_info, ('level3', $tag_l3, 'intron') ) ) {
	    	push @resu, sprintf("%-45s%d%s", "Number of intron in $tag_l3", $all_info->{'level3'}{$tag_l3}{'intron'}{'nb_feat'},"\n");
	    }
	}

	return \@resu;
}

#############
# Give info about shortest feature of each type
sub _info_shortest {

	my ($all_info) = @_  ;
	my @resu;

	#print level1
	foreach my $tag_l1 (sort keys %{$all_info->{'level1'}}){
		push @resu, sprintf("%-45s%d%s", "Shortest $tag_l1", $all_info->{'level1'}{$tag_l1}{'shortest'},"\n");
	}

	#print level2
	foreach my $tag_l2 (sort keys %{$all_info->{'level2'}}){
	    push @resu, sprintf("%-45s%d%s", "Shortest $tag_l2", $all_info->{'level2'}{$tag_l2}{'shortest'},"\n");
	 }

	#print level3
	foreach my $tag_l3 (sort keys %{$all_info->{'level3'}}){
		if( ! exists_keys ( $all_info, ('level3', $tag_l3, 'shortest') ) or $all_info->{'level3'}{$tag_l3}{'shortest'} == 0 ) {
			print "No shortest for $tag_l3\n";
		}
		else{
	    	push @resu, sprintf("%-45s%d%s", "Shortest $tag_l3", $all_info->{'level3'}{$tag_l3}{'shortest'},"\n");
		}
	}

	#print level3 - spread feature cases
	foreach my $tag_l3 (sort keys %{$all_info->{'level3'}}){
		if( exists_keys ($all_info, ('level3', $tag_l3, 'piece') ) ) {
	    	push @resu, sprintf("%-45s%d%s", "Shortest $tag_l3 piece", $all_info->{'level3'}{$tag_l3}{'piece'}{'shortest'},"\n");
	    }
	}

	#print level3 - intron
	foreach my $tag_l3 (sort keys %{$all_info->{'level3'}}){
		if(exists_keys($all_info, ('level3',$tag_l3,'intron'))){
	    	push @resu, sprintf("%-45s%d%s", "Shortest intron into $tag_l3 part", $all_info->{'level3'}{$tag_l3}{'intron'}{'shortest'},"\n");
	    }
	}

	return \@resu;
}

#############
# Give info about longest feature of each type
sub _info_longest {

	my ($all_info) = @_  ;
	my @resu;

	#print level1
	foreach my $tag_l1 (sort keys %{$all_info->{'level1'}}){
		push @resu, sprintf("%-45s%d%s", "Longest $tag_l1", $all_info->{'level1'}{$tag_l1}{'longest'},"\n");
	}

	#print level2
	foreach my $tag_l2 (sort keys %{$all_info->{'level2'}}){
	    push @resu, sprintf("%-45s%d%s", "Longest $tag_l2", $all_info->{'level2'}{$tag_l2}{'longest'},"\n");
	 }

	#print level3
	foreach my $tag_l3 (sort keys %{$all_info->{'level3'}}){
		if( ! exists_keys($all_info, ('level3', $tag_l3, 'longest') ) or $all_info->{'level3'}{$tag_l3}{'longest'} == 0 ) {
			print "No longest for $tag_l3\n";
		}
		else{
	    	push @resu, sprintf("%-45s%d%s", "Longest $tag_l3", $all_info->{'level3'}{$tag_l3}{'longest'},"\n");
		}
	}

	#print level3 - spread feature cases
	foreach my $tag_l3 (sort keys %{$all_info->{'level3'}}){
		if( exists_keys ($all_info, ('level3', $tag_l3, 'piece') ) ) {
	    	push @resu, sprintf("%-45s%d%s", "Longest $tag_l3 piece", $all_info->{'level3'}{$tag_l3}{'piece'}{'longest'},"\n");
	    }
	}

	#print level3 - intron
	foreach my $tag_l3 (sort keys %{$all_info->{'level3'}}){
		if( exists_keys($all_info, ('level3',$tag_l3,'intron'))){
	    	push @resu, sprintf("%-45s%d%s", "Longest intron into $tag_l3 part", $all_info->{'level3'}{$tag_l3}{'intron'}{'longest'},"\n");
	    }
	}

	return \@resu;
}

#############
# Give info about number mean of feature of each type per Parent type (mRNA per gene / cds per mRNA / etc)
sub _info_mean_per {
	my ($all_info) = @_  ;
	my @resu;

	#print level2
	foreach my $tag_l2 (sort keys %{$all_info->{'level2'}}){
		foreach my $tag_l1 (sort keys %{$all_info->{'level1'}}){
			my $mean=  $all_info->{'level2'}{$tag_l2}{'nb_feat'}/$all_info->{'level1'}{$tag_l1}{'nb_feat'};
		    push @resu, sprintf("%-45s%.1f%s", "mean $tag_l2"."s per $tag_l1", $mean,"\n");
		}
	}
	#print level3
	foreach my $tag_l3 (sort keys %{$all_info->{'level3'}}){
	    foreach my $tag_l2 (sort keys %{$all_info->{'level2'}}){
			my $mean=  $all_info->{'level3'}{$tag_l3}{'nb_feat'}/$all_info->{'level2'}{$tag_l2}{'nb_feat'};
		    push @resu, sprintf("%-45s%.1f%s", "mean $tag_l3"."s per $tag_l2", $mean,"\n");
		}
	}
	#print level3 - spread feature cases
	foreach my $tag_l3 (sort keys %{$all_info->{'level3'}}){
		if( exists_keys ($all_info, ('level3', $tag_l3, 'exon') ) ) {
			my $mean=  $all_info->{'level3'}{$tag_l3}{'exon'}{'nb_feat'}/$all_info->{'level3'}{$tag_l3}{'nb_feat'};
		    push @resu, sprintf("%-45s%.1f%s", "mean exons per $tag_l3", $mean,"\n");
	    }
	}
	#print introns
	foreach my $tag_l3 (sort keys %{$all_info->{'level3'}}){
		foreach my $tag_l2 (sort keys %{$all_info->{'level2'}}){
			if(exists_keys($all_info, ('level3',$tag_l3,'intron'))){
			    my $mean=  $all_info->{'level3'}{$tag_l3}{'intron'}{'nb_feat'}/$all_info->{'level2'}{$tag_l2}{'nb_feat'};
			    push @resu, sprintf("%-45s%.1f%s", "mean introns in $tag_l3"."s per $tag_l2", $mean,"\n");
			}
		}
	}
	return \@resu;
}

#############
# Give info about lenght of the total of features by type
sub _info_length {
	my ($all_info) = @_  ;
	my @resu;

	#print level1
	foreach my $tag_l1 (sort keys %{$all_info->{'level1'}}){
		push @resu, sprintf("%-45s%d%s", "Total $tag_l1 length", $all_info->{'level1'}{$tag_l1}{'size_feat'},"\n");
	}

	#print level2
	foreach my $tag_l2 (sort keys %{$all_info->{'level2'}}){
	    push @resu, sprintf("%-45s%d%s", "Total $tag_l2 length", $all_info->{'level2'}{$tag_l2}{'size_feat'},"\n");
	 }

	#print level3
	foreach my $tag_l3 (sort keys %{$all_info->{'level3'}}){
	    push @resu, sprintf("%-45s%d%s", "Total $tag_l3 length", $all_info->{'level3'}{$tag_l3}{'size_feat'},"\n");
	}

	#print introns
	foreach my $tag_l3 (sort keys %{$all_info->{'level3'}}){
		if(exists_keys($all_info, ('level3',$tag_l3,'intron'))){
	    	push @resu, sprintf("%-45s%d%s", "Total intron length per $tag_l3", $all_info->{'level3'}{$tag_l3}{'intron'}{'size_feat'},"\n");
	    }
	}

	return \@resu;
}

#############
# Give info about mean lenght of features by type
sub _info_mean_length {
	my ($all_info) = @_  ;
	my @resu;


	#print level1
	foreach my $tag_l1 (sort keys %{$all_info->{'level1'}}){
		my $meanl= $all_info->{'level1'}{$tag_l1}{'size_feat'}/$all_info->{'level1'}{$tag_l1}{'nb_feat'};
		push @resu, sprintf("%-45s%d%s", "mean $tag_l1 length", $meanl,"\n");
	}

	#print level2
	foreach my $tag_l2 (sort keys %{$all_info->{'level2'}}){
		my $size_feat = $all_info->{'level2'}{$tag_l2}{'size_feat'};
		my $nb_feat = $all_info->{'level2'}{$tag_l2}{'nb_feat'};

		if($size_feat !=0 and $nb_feat != 0){
			my $meanl= $all_info->{'level2'}{$tag_l2}{'size_feat'}/$all_info->{'level2'}{$tag_l2}{'nb_feat'};
		    push @resu, sprintf("%-45s%d%s", "mean $tag_l2 length", $meanl,"\n");
		}
		else{warn "Problem in the calcul of level2 - $tag_l2 - size_feat";}
	 }

	#print level3
	foreach my $tag_l3 (sort keys %{$all_info->{'level3'}}){
		if( ($all_info->{'level3'}{$tag_l3}{'size_feat'} == 0) or ($all_info->{'level3'}{$tag_l3}{'nb_feat'} == 0) ) {
			print "No size_feat for $tag_l3\n";
		}
		else{
			my $meanl= $all_info->{'level3'}{$tag_l3}{'size_feat'}/$all_info->{'level3'}{$tag_l3}{'nb_feat'};
		    push @resu, sprintf("%-45s%d%s", "mean $tag_l3 length", $meanl,"\n");
		}
	}

	#print level3 - multifeature cases
	foreach my $tag_l3 (sort keys %{$all_info->{'level3'}}){
		if( exists_keys ($all_info, ('level3', $tag_l3, 'exon') ) ) {
			my $meanl= $all_info->{'level3'}{$tag_l3}{'size_feat'}/$all_info->{'level3'}{$tag_l3}{'exon'}{'nb_feat'};
	    	push @resu, sprintf("%-45s%d%s", "mean $tag_l3 piece length", $meanl,"\n");
	    }
	}

	#print introns
	foreach my $tag_l3 (sort keys %{$all_info->{'level3'}}){
		if(exists_keys($all_info, ('level3',$tag_l3,'intron'))){
	    	my $meanl= $all_info->{'level3'}{$tag_l3}{'intron'}{'size_feat'}/$all_info->{'level3'}{$tag_l3}{'intron'}{'nb_feat'};
	    	push @resu, sprintf("%-45s%d%s", "mean intron in $tag_l3 length", $meanl,"\n");
	    }
	}

	return \@resu;
}

#############
# Give info about the features' coverage (by types) within/among the genome
sub _info_coverage {
	my ($all_info, $genomeSize) = @_  ;
	my @resu;

	#print level1
	foreach my $tag_l1 (sort keys %{$all_info->{'level1'}}){
		my $perc= ($all_info->{'level1'}{$tag_l1}{'size_feat'}*100)/$genomeSize;
		push @resu, sprintf("%-45s%.1f%s", "% of genome covered by $tag_l1", $perc,"\n");
	}

	#print level2
	foreach my $tag_l2 (sort keys %{$all_info->{'level2'}}){
		my $perc= ($all_info->{'level2'}{$tag_l2}{'size_feat'}*100)/$genomeSize;
	    push @resu, sprintf("%-45s%.1f%s", "% of genome covered by $tag_l2", $perc,"\n");
	 }

	#print level3
	foreach my $tag_l3 (sort keys %{$all_info->{'level3'}}){
	    my $perc= ($all_info->{'level3'}{$tag_l3}{'size_feat'}*100)/$genomeSize;
	    push @resu, sprintf("%-45s%.1f%s", "% of genome covered by $tag_l3", $perc,"\n");
	}

	#print level3
	foreach my $tag_l3 (sort keys %{$all_info->{'level3'}}){
		if(exists_keys($all_info, ('level3',$tag_l3,'intron'))){
		    my $perc= ($all_info->{'level3'}{$tag_l3}{'intron'}{'size_feat'}*100)/$genomeSize;
		    push @resu, sprintf("%-45s%.1f%s", "% of genome covered by intron from $tag_l3", $perc,"\n");
		}
	}

	return \@resu;
}

#############
# @Purpose: Give info about the overlaping genes
# Could check only at L1 if level1 only feature (standalone, topfeature)
# check overlap at L3 if not a level1 only feature (standalone, topfeature)
# return a string summerizing the number of overlap
# @input: 3 => hash(omniscient hash), featureL2,  primary tag, if it is for L1 feature only (topfeature and standalone feature)
# @output: 1 => int (nb feature removed)
sub _detect_overlap_features{
	my ($omniscient, $sortBySeq, $toplevel) = @_;
	my $resume_case = 0;

	foreach my $locusID ( keys %{$sortBySeq}){ # tag_l1 = gene or repeat etc...
		foreach my $tag_l1 ( keys %{$sortBySeq->{$locusID}} ) {

			#create list to keep track of l1
			my %to_check;
			foreach my $feature_l1 ( @{$sortBySeq->{$locusID}{$tag_l1}} ) {
				my $id_l1 = lc($feature_l1->_tag_value('ID'));
				$to_check{$id_l1}++;
			}

			# Go through location from left to right ###
			while ( @{$sortBySeq->{$locusID}{$tag_l1}} ){

				my $feature_l1 = shift @{$sortBySeq->{$locusID}{$tag_l1}};
				my $id_l1 = lc($feature_l1->_tag_value('ID'));
				my @location = ($id_l1, int($feature_l1->start()), int($feature_l1->end())); # This location will be updated on the fly

				# Go through location from left to right ### !!
				foreach my $l1_feature2 ( @{$sortBySeq->{$locusID}{$tag_l1}} ) {
					my $id2_l1 = lc($l1_feature2->_tag_value('ID'));
					my @location_to_check = ($id2_l1, int($l1_feature2->start()), int($l1_feature2->end()));

					#If location_to_check start if over the end of the reference location, we stop
					if($location_to_check[1] > $location[2]) {last;}

					# Let's check at Gene LEVEL
					if(location_overlap(\@location, \@location_to_check)){

						# Need to check only at level1
						if($toplevel){
							$resume_case++;
						}
						# Need to check at level3
						else{
							#let's check at CDS level
							if(check_gene_overlap_at_level3($omniscient, $omniscient , $id_l1, $id2_l1, "exon")){ #If contains CDS it has to overlap at CDS level to be merged, otherwise any type of feature level3 overlaping is sufficient to decide to merge the level1 together
								#they overlap in the CDS we should give them the same name
								$resume_case++;
							}
						}
					}
				}
			}
		}
	}
	return $resume_case;
}

1;
