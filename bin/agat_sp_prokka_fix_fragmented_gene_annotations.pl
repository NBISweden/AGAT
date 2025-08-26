#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Clone 'clone';
use File::Basename;
use Pod::Usage;
use LWP::UserAgent;
use List::MoreUtils qw(uniq);
use Sort::Naturally;
use Bio::DB::Fasta;
use Bio::SeqIO;
use AGAT::AGAT;

BEGIN {
   package case_info;
   use Moose;

   has name => (is => 'rw', isa => 'Str', required => 1);
	 has expected_length => (is => 'rw', isa => 'Int');
	 has expected_length_10 => (is => 'rw', isa => 'Int');
	 has current_aa_length_together => (is => 'rw', isa => 'Int');
	 has list_aa_sizes => (is => 'rw', isa => 'Int');
	 has list_gene => ('is' => 'rw', isa => 'ArrayRef');
	 has list_gene_feature => ('is' => 'rw', isa => 'ArrayRef');
	 has hash_sub_gene_obj => ('is' => 'rw', isa => 'ArrayRef');
	 has insert_size => ('is' => 'rw', isa => 'ArrayRef');
	 has insert_position => ('is' => 'rw', isa => 'ArrayRef');
	 has frame => ('is' => 'rw', isa => 'ArrayRef');
	 has total_insert_before => ('is' => 'rw', isa => 'ArrayRef');

   $INC{"case_info.pm"} = 1;

	 package sub_gene;
   use Moose;

   has id => (is => 'rw', isa => 'Str', required => 1);
	 has inference_db => (is => 'rw', isa => 'Str');
	 has inference_value => (is => 'rw', isa => 'Str');
	 has inference_aa_length => (is => 'rw', isa => 'Int');
	 has inference_aa_seq => (is => 'rw', isa => 'Str');
	 has inference_dna_length => (is => 'rw', isa => 'Int');
	 has inference_dna_seq => (is => 'rw', isa => 'Str');
	 has overlap_dna => (is => 'rw', isa => 'Int');
	 has current_aa_length => (is => 'rw', isa => 'Int');
	 has current_aa_seq => (is => 'rw', isa => 'Str');
	 has current_dna_length => (is => 'rw', isa => 'Int');
	 has current_dna_seq => (is => 'rw', isa => 'Str');

   $INC{"sub_gene.pm"} = 1;
}

use case_info;

# maximum percent size over original protein size to merge chunks
my $SIZE_OPT=21;

my $header = get_agat_header();
my ( $opt, $usage, $config ) = AGAT::AGAT::describe_script_options(
    $header,
    [ 'gff=s',          'Input reference gff file',   { required => 1 } ],
    [ 'fasta|fa|f=s',   'Input reference fasta file', { required => 1 } ],
    [ 'db=s',           'Input db file',              { required => 1 } ],
    [ 'frags!',         'Output fragments' ],
    [ 'pseudo!',        'Output pseudogenes' ],
    [ 'hamap_size=s',   'HAMAP size category', { default => 'high',
        callbacks => { allowed => sub { $_[0] =~ /^(?:high|low|middle)$/ or die 'Wrong value provided for the option --hamap_size: accepted value: high, low or middle'; } } } ],
    [ 'table|codon|ct=i', 'Codon translation table', { default => 1,
        callbacks => { positive => sub { $_[0] > 0 or die 'Codon translation table must be positive'; } } } ],
    [ 'skip_hamap!',    'Skip hamap processing' ],
);

my $gff        = $opt->gff;
my $file_fasta = $opt->fasta;
my $file_db    = $opt->db;
my $frags      = $opt->frags;
my $pseudo     = $opt->pseudo;
my $hamap_size = $opt->hamap_size;
my $codonTable = $opt->table;
my $skip_hamap = $opt->skip_hamap;
my $outfolder  = $opt->out;
my $verbose    = $config->{verbose};

my @copyARGV = @ARGV;

my $log;
if ( my $log_name = $config->{log_path} ) {
    open( $log, '>', $log_name ) or die "Can not open $log_name for printing: $!";
}
dual_print($log, $header, 3);

# Check codon table
$codonTable = get_proper_codon_table($codonTable, $log, $verbose);

######################
# Manage output file #

my $gff_out; my $fasta_out; my $report; #my $gffout4;
my ($file_fasta_in,$path_fasta,$ext_fasta) = fileparse($file_fasta,qr/\.[^.]*/);
my ($file_gff_in,$path_gff,$ext_gff) = fileparse($gff,qr/\.[^.]*/);

if ($outfolder) {
        if (-d $outfolder){
                  my $msg = "Provided output folder exists. Exit!\n";
                  dual_warn($log, $msg, 3);
                  exit;
        }
        mkdir $outfolder;

	if($frags or $pseudo){
		# gff out
		my $gff_out_path = "$outfolder/$file_gff_in$ext_gff";

		$gff_out = prepare_gffout($config, $gff_out_path);
	}
	if($frags){
		# fasta out
		my $fasta_out_path = "$outfolder/$file_fasta_in$ext_fasta";
		open(my $fh_fasta, '>', $fasta_out_path) or die "Could not open file '$fasta_out_path' $!";
		$fasta_out=  Bio::SeqIO->new(-fh => $fh_fasta , -format => 'Fasta');
	}

	# report out
	my $report_out_path = "$outfolder/report.txt";
	open($report, '>', $report_out_path) or die "Could not open file '$report_out_path' $!";
}
else{
  my $msg = "No output folder provided. Exit!\n";
  dual_warn($log, $msg, 3);
  exit;
}
$hamap_size = lc($hamap_size);

                #####################
                #     MAIN          #
                #####################


######################
### Parse GFF input #
my ($hash_omniscient, $hash_mRNAGeneLink) =slurp_gff3_file_JD({ input => $gff,
                                                                config => $config
                                                              });
dual_print($log, "GFF3 file parsed\n");


####################
# index the genome #
my $db_fasta = Bio::DB::Fasta->new($file_fasta);
# save ID in lower case to avoid cast problems
my %all_db_fasta_IDs;
my @ids_db_fasta     = $db_fasta->get_all_primary_ids;
foreach my $id (@ids_db_fasta ){$all_db_fasta_IDs{lc($id)}=$id;}
dual_print($log, "Fasta file parsed\n");


my $db_db = Bio::DB::Fasta->new($file_db);
# save ID in lower case to avoid cast problems
my %all_db_db_IDs;
my @ids_db_db     = $db_db->get_all_primary_ids;
foreach my $id (@ids_db_db ){$all_db_db_IDs{lc($id)}=$id;}
dual_print($log, "db file parsed\n");

####################

my $total_gene = 0;
my $total_gene_with_name = 0;
my $total_putative_gene_fixed_before = 0;
my $total_putative_gene_fixed_after = 0;
my $total_gene_fixed_before = 0;
my $total_gene_fixed_after = 0;
my $fixed_gene = 0;
my $case_frame1 = 0; # case where to contiguous gene cannot be merged because in hte same frame. It means the stop codon of the first gene cannot be skipped.
my $case_frame2 = 0;
my $case_frame3 = 0;
my $incongruent_size = 0;
my $congruent_size = 0;
my $total_insert_before = 0;
my %gff_shift;

# sort by seq id
my ( $hash_sortBySeq, $hash_sortBySeq_std, $hash_sortBySeq_topf ) = collect_l1_info_sorted_by_seqid_and_location($hash_omniscient);

# Read by seqId to sort properly the output by seq ID
foreach my $seqid (sort { (($a =~ /(\d+)$/)[0] || 0) <=> (($b =~ /(\d+)$/)[0] || 0) } keys %{$hash_sortBySeq}){ # loop over all the feature level1

	my $previousName=undef;
	my $previous_gene_feature=undef;
	my @listGeneToMerge;
	my @listGeneNameToMerge;
	my $Name_ok = undef;
	$total_insert_before = 0;

	my $seq_id_correct = $all_db_fasta_IDs{lc($seqid)};
	# get description
	my $header = $db_fasta->header($seq_id_correct);
	my @headers =  split(' ',$header);
	shift @headers;
	my $description = "";
	if(@headers){
		$description .= join(' ',@headers);
	}
	#print "seq description: $description\n";
	my $sequence = $db_fasta->seq($seq_id_correct);
	my $seqObj  = Bio::Seq->new( '-format' => 'fasta' , -id => $seq_id_correct, -seq => $sequence, -description => $description);

	foreach my $locationid ( sort { ncmp ($a, $b) } keys %{$hash_sortBySeq->{$seqid} } ){

		$total_gene++;
		my $primary_tag_l1 = $hash_sortBySeq->{$seqid}{$locationid}{'tag'};
		my $id_tag_key_level1 = $hash_sortBySeq->{$seqid}{$locationid}{'id'};
		my $gene_feature = $hash_omniscient->{'level1'}{$primary_tag_l1}{$id_tag_key_level1};


		#get Name value's attribute
		$Name_ok = undef;
		if ($gene_feature->has_tag("Name")){
			$total_gene_with_name ++;
			my $Name_att = $gene_feature->_tag_value("Name");
			if ( $Name_att =~ /(.*)_[0-9]$/){
				$Name_ok = $1;
			}
		}

		#check name and strand is the same
		if ($Name_ok and $previousName and ($Name_ok eq $previousName) and ($gene_feature->strand eq $previous_gene_feature->strand) ){
			#add the gene to list of gene to merge together
			if (! @listGeneToMerge){ # only first time
				push @listGeneToMerge,$previous_gene_feature;
				push @listGeneNameToMerge, $previous_gene_feature->_tag_value("Name");
			}
			push @listGeneToMerge,$gene_feature;
			push @listGeneNameToMerge, $gene_feature->_tag_value("Name");
		}
		else{
		  #We have a case to process
			if(@listGeneToMerge){
				my @listGeneToMerge_sorted = sort {$a->start <=> $b->start } @listGeneToMerge;
				my $obj_case = case_info->new(name => $previousName, list_gene =>\@listGeneNameToMerge, list_gene_feature => \@listGeneToMerge_sorted);

				if ( check_protein_size_congruency($obj_case, $seqObj) ){
					if (check_long_orf_if_can_be_merged($obj_case)){
						$total_gene_fixed_before += @listGeneToMerge;
						$total_gene_fixed_after ++;
						if($frags or $pseudo){
							merge_case($obj_case, $seqObj);
						}
					}
				}
			}
			@listGeneToMerge=();
			@listGeneNameToMerge=();
		}
		#
		$previousName = $Name_ok;
		$previous_gene_feature = $gene_feature;
	}
	# deal case where it was the last element
	if(@listGeneToMerge){
		my @listGeneToMerge_sorted = sort {$a->start <=> $b->start } @listGeneToMerge;
		my $obj_case = case_info->new(name => $previousName, list_gene =>\@listGeneNameToMerge, list_gene_feature => \@listGeneToMerge_sorted);

		if( check_protein_size_congruency($obj_case, $seqObj) ){
			if ( check_long_orf_if_can_be_merged($obj_case) ) {
				$total_gene_fixed_before += @listGeneToMerge;
				$total_gene_fixed_after ++;
				if($frags or $pseudo){
					merge_case($obj_case, $seqObj);
				}
			}
		}
	}

	if($frags){
		# write the tested seq modified or intact
		$fasta_out->write_seq($seqObj);
	}
}


if($frags){
	# add non modified sequences
	foreach my $id_seq (keys %all_db_fasta_IDs){
		my $seq_id_correct = $all_db_fasta_IDs{lc($id_seq)};

		if (! exists_keys( $hash_sortBySeq, ($seq_id_correct) )){
				my $sequence = $db_fasta->seq($seq_id_correct);

				# get description
				my $header = $db_fasta->header($seq_id_correct);
				my @headers =  split(' ',$header);
				shift @headers;
				my $description = "";
				if(@headers){
					$description .= join(' ',@headers);
				}

				my $seqObj  = Bio::Seq->new( '-format' => 'fasta' , -id => $seq_id_correct, -seq => $sequence, -description => $description);
				$fasta_out->write_seq($seqObj);
		}
	}

	#shift annotation
	# need to be parsed again because we might have removed some features
	( $hash_sortBySeq, $hash_sortBySeq_std, $hash_sortBySeq_topf ) = collect_l1_info_sorted_by_seqid_and_location($hash_omniscient);
	foreach my $seqid (sort { (($a =~ /(\d+)$/)[0] || 0) <=> (($b =~ /(\d+)$/)[0] || 0) } keys %{$hash_sortBySeq}){

		my $total_shift = 0;
		my $shift_location;
		my @shift_locations;

		# make a list of location
		foreach my $val (keys %{$gff_shift{$seqid}}) {
			push @shift_locations, $val;
		}
		@shift_locations = sort { $a <=> $b} @shift_locations; # sort
		$shift_location = shift @shift_locations; # get first value
		if (exists_keys (\%gff_shift, ($seqid) ) ){
			# loop over feature in order
			foreach my $locationid ( sort { ncmp ($a, $b) } keys %{$hash_sortBySeq->{$seqid} } ){

				my $primary_tag_l1 = $hash_sortBySeq->{$seqid}{$locationid}{'tag'};
				my $id_tag_key_level1 = $hash_sortBySeq->{$seqid}{$locationid}{'id'};
				my $gene_feature = $hash_omniscient->{'level1'}{$primary_tag_l1}{$id_tag_key_level1};
				my $id_l1 = lc($gene_feature->_tag_value('ID'));

				# get new start
				while (@shift_locations and $shift_location < $gene_feature->start) {
					$total_shift += $gff_shift{$seqid}{$shift_location}; # append value
					$shift_location = shift @shift_locations; # get first value
				}
				if($shift_location and $shift_location < $gene_feature->start and !@shift_locations){
					$total_shift += $gff_shift{$seqid}{$shift_location};
					$shift_location = undef;
				}
				$gene_feature->start($gene_feature->start + $total_shift);

				# get new end
				while (@shift_locations and $shift_location <= $gene_feature->end) {
					$total_shift += $gff_shift{$seqid}{$shift_location}; # append value
					$shift_location = shift @shift_locations; # get first value
				}
				if($shift_location and $shift_location <= $gene_feature->end and !@shift_locations){
					$total_shift += $gff_shift{$seqid}{$shift_location};
					$shift_location = undef;
				}
				$gene_feature->end($gene_feature->end + $total_shift);

				foreach my $tag_l2 (sort keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
					if ( exists_keys( $hash_omniscient, ('level2', $tag_l2, $id_l1) ) ){
						foreach my $feature_l2 ( @{$hash_omniscient->{'level2'}{$tag_l2}{$id_l1}}){
							$feature_l2->start( $gene_feature->start );
							$feature_l2->end( $gene_feature->end );
							my $id_l2 = lc($feature_l2->_tag_value('ID'));

							foreach my $tag_l3 (sort keys %{$hash_omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...
								if ( exists_keys( $hash_omniscient, ('level3', $tag_l3, $id_l2) ) ){
									foreach my $feature_l3 ( @{$hash_omniscient->{'level3'}{$tag_l3}{$id_l2}} ){
										$feature_l3->start( $gene_feature->start );
										$feature_l3->end( $gene_feature->end );
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

if($frags or $pseudo){
	#print gff
	print_omniscient( {omniscient => $hash_omniscient, output => $gff_out} );
}

# print results
my $stringprint = "total gene: $total_gene\n";
$stringprint .= "total gene with name: $total_gene_with_name\n";
$stringprint .=  "total gene not checked:".($total_gene-$total_gene_with_name)."\n";
$stringprint .=  "Among the $total_gene_with_name gene checked, $total_putative_gene_fixed_before putative gene may be part of $total_putative_gene_fixed_after FRAGS.\n";
$stringprint .=  "$incongruent_size case where contiguous genes that do not pass the length test. When appending them we got something too long compared to the reference protein.\n";
$stringprint .=  "$congruent_size detected FRAGS:\n";
$stringprint .=  "* We got $case_frame1 cases where the FRAGS is in the same frame.\n";
$stringprint .=  "* We got $case_frame2 cases where the FRAGS is due to a frameshift of 1 nucleotide.\n";
$stringprint .=  "* We got $case_frame3 cases where the FRAGS is due to a frameshift of 2 nucleotides.\n";

if ($pseudo and $frags) {
	$stringprint .=  "$case_frame1 FRAGS merged and annotated as pseudo.\n";
	my $cases = $case_frame2 + $case_frame3;
	$stringprint .=  "$cases FRAGS merged and frameshift fixed.\n";
	$stringprint .=  "In total $total_gene_fixed_before gene have been merged into $total_gene_fixed_after\n";
}
elsif ($pseudo and ! $frags){
	my $cases = $case_frame1 + $case_frame2 + $case_frame3;
	$stringprint .=  "$cases FRAGS merged and annotated as pseudo.\n";
}
elsif (! $pseudo and $frags){
	my $cases = $case_frame2 + $case_frame3;
	$stringprint .=  "$cases FRAGS merged and frameshift fixed.\n";
	$stringprint .=  "In total $total_gene_fixed_before gene have been merged into $total_gene_fixed_after\n";
}

$stringprint .= "\nWhen a FRAGS is detected you shoud think about few things...\n".
								"Does the stop codon really codes for a stop codon?\n".
								"Is the stop codon from the first gene real?\n".
								"    => If the two fragments (genes) are in the same frame, it can be due to a substitution...\n".
								"    => If the two fragments (genes) are in a different frame, it can be due to an indel...\n".
								"Indel or substitution, is it a pseudogene or a sequencing error?\n".
								"In a case of pseudogene, if the pseudogene is not recent, other mutations should be accumulated in the sequence after the first premature stop codon...\n";
$stringprint .=  "\nBye Bye.\n";

dual_print($log, $stringprint);
print $report $stringprint;


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

sub check_protein_size_congruency{
	my ($obj_case, $seqObj) = @_;

	my $size_congruency;
	my $Name_ok = $obj_case->{name};
	my $listGeneToMerge = $obj_case->{list_gene_feature};
	my $listGeneNameToMerge =  $obj_case->{list_gene};

	$total_putative_gene_fixed_before += @$listGeneToMerge;
	$total_putative_gene_fixed_after ++;

        dual_print($log, "\nBased on their names and the fact they are conitguous, those ".@$listGeneToMerge." genes might be merged: @$listGeneNameToMerge \n");
  retrieve_expected_protein_length($obj_case);

	my $expected_protein_length = $obj_case->{expected_length};
	if ($expected_protein_length){
                dual_print($log, "Average of the expected length = $expected_protein_length \n");

	 	#add 10 percent to expected protein $seq_lengt
	 	my $expected_protein_length_plus10 = int (($expected_protein_length * 110) / 100);
		$obj_case->{expected_length_10} = $expected_protein_length_plus10;
                dual_print($log, "Average of the expected length + 20 % = $expected_protein_length_plus10 \n");

		# need to get the oervlap to remove from the calculated length
		my $aa_overlap = get_overlap($obj_case);

		#Compute current size of predicted proteins put together
		my $total_current_size;
	 	foreach my $gene_feature (@{$listGeneToMerge}) {
	 		$total_current_size += (($gene_feature->end - $gene_feature->start + 1) / 3);
	 	}
		$total_current_size -= $aa_overlap; # remove all overlap parts
		$obj_case->{current_aa_length_together} = $total_current_size;
                dual_print($log, "current length adding all genes together (and removing overlaping part): $total_current_size\n");

	 	if ($total_current_size < $expected_protein_length_plus10){
                        dual_print($log, "$total_current_size < $expected_protein_length_plus10 => Let's merge them. (The length of the appended proteins is shorter than the size of the protein use for the inference)\n\n");
			$congruent_size++;
			$size_congruency = 1;
		}
		else{
			$incongruent_size++;
		}
	}
	else{
                dual_print($log, "No expected size found - skip the case\n");
	}
 return $size_congruency;
}


# We have a case that passed the tests. We can update the fasta and the gff
sub merge_case{
	my ($obj_case, $seqObj)=@_;

	my $listGeneToMerge = $obj_case->{list_gene_feature};
	my @listGene =  @$listGeneToMerge;

	while (@listGene > 1){
		my $gene_feature1 = shift @listGene;
		my $gene_feature2 = $listGene[0];

		# deal with the fasta sequence
		my $seq_id = lc($gene_feature1->seq_id);
		my $value_insert_before = shift @{$obj_case->{total_insert_before}};
		my $value_insert_size = shift  @{$obj_case->{insert_size}};
		my $value_insert_position = shift  @{$obj_case->{insert_position}};
		$gff_shift{$gene_feature1->seq_id}{$value_insert_position}=$value_insert_size; #keep track of shifts
		my $sequence = $seqObj->seq();

                dual_print($log, "add N x $value_insert_size at position $value_insert_position\n", 2);
                dual_print($log, "piece of sequence before: ".substr($sequence, $value_insert_before + $value_insert_position - 10, 20)."\n", 2);
                substr($sequence, $value_insert_before + $value_insert_position-1, 1) = "N" x ($value_insert_size+1);
                dual_print($log, "piece of sequence after: ".substr($sequence, $value_insert_before + $value_insert_position - 10, 20)."\n", 2);
		$seqObj->seq($sequence);

		# Should I add Pseudo attribure?
		my $add_pseudo = undef;
		if($pseudo) {
			$add_pseudo = 1;
			if($frags){
				$add_pseudo = undef;
				# if frags and case frame 1, we have to add the pseudo attribute
				if( grep( /1/, @{$obj_case->{frame}} ) ){
					$add_pseudo = 1;
				}
			}
		}

		# create first stop codon position
		my $first_stop = ($gene_feature1->strand == -1 or $gene_feature1->strand eq "-") ? $gene_feature2->start + 2 : $gene_feature1->end - 2;

		# modify end position
		$gene_feature1->add_tag_value('agat_pseudo', $first_stop) if ($add_pseudo);
		$gene_feature1->end( $gene_feature2->end );
		my $id_l1 = lc($gene_feature1->_tag_value('ID'));
		foreach my $tag_l2 (sort keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
			if ( exists_keys( $hash_omniscient, ('level2', $tag_l2, $id_l1) ) ){
				foreach my $feature_l2 ( @{$hash_omniscient->{'level2'}{$tag_l2}{$id_l1}}){
					$feature_l2->add_tag_value('agat_pseudo', $first_stop ) if ($add_pseudo);
					$feature_l2->end( $gene_feature2->end );
					my $id_l2 = lc($feature_l2->_tag_value('ID'));
					foreach my $tag_l3 (sort keys %{$hash_omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...
						if ( exists_keys( $hash_omniscient, ('level3', $tag_l3, $id_l2) ) ){
							foreach my $feature_l3 ( @{$hash_omniscient->{'level3'}{$tag_l3}{$id_l2}} ){
								$feature_l3->add_tag_value('agat_pseudo', $first_stop ) if ($add_pseudo);
								$feature_l3->end( $gene_feature2->end );
							}
						}
					}
				}
			}
		}
		# remove gene2
		remove_l1_and_relatives($hash_omniscient, $gene_feature2);
	}
}

# test if we can merge genes
# test will look if we can have one only long ORF by shifting the frame by adding 1 or 2 nucleotides.
# To avoid the stop codon we also replace the last nucleotide from the codon stop of the first gene by a N
sub check_long_orf_if_can_be_merged{
	my ($obj_case)=@_;
	my $merged = 0;

	my $listGeneToMerge = $obj_case->{list_gene_feature}; #list is sorted
	my @listGene =  @$listGeneToMerge;

	while (@listGene > 1){
		my $gene_feature1 = shift @listGene;
		my $gene_feature2 = $listGene[0];
		my $intergenic = 1;
                dual_print($log, $gene_feature1->gff_string()."\n", 2);
                dual_print($log, $gene_feature2->gff_string()."\n", 2);
                if($gene_feature1->end > $gene_feature2->start){
                        dual_print($log, "No intergenic region!\n");
                        $intergenic=undef;

                }
                else{
                        my $intergenic_region = [$gene_feature1->end,$gene_feature2->start];
                        dual_print($log, "intergenic_region = @$intergenic_region\n");
                }

		my $ID_correct = $all_db_fasta_IDs{lc($gene_feature1->seq_id)};

		# ---- get full prot1 ---
		my $subseq1_cds_obj_full;
		my $subseq1_full = $db_fasta->seq($ID_correct, $gene_feature1->start, $gene_feature1->end); #allow to remove the start codon M that can be another AA when
		$subseq1_cds_obj_full = Bio::Seq->new(-seq => $subseq1_full, -alphabet => 'dna' );
		# revcom if minus sstrand
		if($gene_feature1->strand == -1 or $gene_feature1->strand eq "-"){
			$subseq1_cds_obj_full = $subseq1_cds_obj_full->revcom();
		}
		# now translate
		my $subseq1_prot_obj_full = $subseq1_cds_obj_full->translate(-codontable_id => $codonTable);

		# --- get protein part to check ----
		my $subseq1_cds_obj;
		if( $intergenic ){
			my $subseq1 = $db_fasta->seq($ID_correct, $gene_feature1->start, $gene_feature1->end); #allow to remove the start codon M that can be another AA when
			$subseq1_cds_obj = Bio::Seq->new(-seq => $subseq1, -alphabet => 'dna' );

			# revcom if minus strand
			if($gene_feature1->strand == -1 or $gene_feature1->strand eq "-"){
				$subseq1_cds_obj = $subseq1_cds_obj->revcom();
			}
		}
		# ! $intergenic
		else{
			if ($gene_feature1->strand == -1 or $gene_feature1->strand eq "-") {

				# --- need to remove overlap part to test if protein contained in the resulting new protein
				# --- this overlaping part will be used as intergenec region to shift frame

				my $overlap_part = $gene_feature1->end - $gene_feature2->start + 1;
                                dual_print($log, "overlap_part $overlap_part\n", 2);
                                my $to_shrink = 3 + ( 3 - ($overlap_part % 3) );
                                dual_print($log, "to_shrink (last codon plus offset needed to be in frame:) $to_shrink\n", 2);
				my $subseq1 = $db_fasta->seq($ID_correct, $gene_feature1->start, ($gene_feature2->start - $to_shrink - 1 ));
				$subseq1_cds_obj = Bio::Seq->new(-seq => $subseq1, -alphabet => 'dna' );
				$subseq1_cds_obj = $subseq1_cds_obj->revcom();
			}
			#plus strand
			else{
				my $overlap_part = $gene_feature1->end - $gene_feature2->start + 1;
				my $to_shrink = 3 + ( 3 - ($overlap_part % 3) );
				my $subseq1 = $db_fasta->seq($ID_correct, $gene_feature1->start, ($gene_feature2->start - $to_shrink - 1 ));
				$subseq1_cds_obj = Bio::Seq->new(-seq => $subseq1, -alphabet => 'dna' );
			}
		}
		# now translate
		my $subseq1_prot_obj = $subseq1_cds_obj->translate(-codontable_id => $codonTable) ;


		#get protein to second gene
		# ---- get full prot2 ---
		my $subseq2_cds_obj_full;
		my $subseq2_full = $db_fasta->seq($ID_correct, $gene_feature2->start, $gene_feature2->end); #allow to remove the start codon M that can be another AA when
		$subseq2_cds_obj_full = Bio::Seq->new(-seq => $subseq2_full, -alphabet => 'dna' );
		# revcom if minus sstrand
		if($gene_feature2->strand == -1 or $gene_feature2->strand eq "-"){
			$subseq2_cds_obj_full = $subseq2_cds_obj_full->revcom();
		}
		# now translate
		my $subseq2_prot_obj_full = $subseq2_cds_obj_full->translate(-codontable_id => $codonTable) ;

		# ---- get protein part to check ----
		my $subseq2_cds_obj;
		if( $intergenic ){
			my $subseq2 = $db_fasta->seq($ID_correct, $gene_feature2->start, $gene_feature2->end); #allow to remove the start codon M that can be another AA when
			$subseq2_cds_obj = Bio::Seq->new(-seq => $subseq2, -alphabet => 'dna' );

			# revcom if minus sstrand
			if($gene_feature2->strand == -1 or $gene_feature2->strand eq "-"){
				$subseq2_cds_obj = $subseq2_cds_obj->revcom();
			}
		}
		# ! $intergenic
		else{
			if ($gene_feature2->strand == -1 or $gene_feature2->strand eq "-") {
				my $overlap_part = $gene_feature1->end - $gene_feature2->start + 1;
				my $to_shrink = 3 + ( 3 - ( $overlap_part %3 ) );
				my $subseq2 = $db_fasta->seq($ID_correct, $gene_feature1->end + $to_shrink + 1, $gene_feature2->end );
				$subseq2_cds_obj = Bio::Seq->new(-seq => $subseq2, -alphabet => 'dna' );
				$subseq2_cds_obj = $subseq2_cds_obj->revcom();
			}
			else{
				my $overlap_part = $gene_feature1->end - $gene_feature2->start + 1;
                                dual_print($log, "overlap_part $overlap_part\n", 2);
                                my $to_shrink = 3 + (  3 - ( $overlap_part % 3 ) ) ;
                                dual_print($log, "to_shrink (last codon plus offset needed to be in frame:) $to_shrink\n", 2);
				my $subseq2 = $db_fasta->seq($ID_correct, $gene_feature1->end + $to_shrink + 1, $gene_feature2->end );
				$subseq2_cds_obj = Bio::Seq->new(-seq => $subseq2, -alphabet => 'dna' );
			}
		}
		# now translate
		my $subseq2_prot_obj = $subseq2_cds_obj->translate(-codontable_id => $codonTable) ;

		# -----------------------------------------------
		# ---- strand- ----
		if ($gene_feature1->strand == -1 or $gene_feature1->strand eq "-"){
			my $found = 0;
			my $prot_second=$subseq1_prot_obj;
                        dual_print($log, "Minus strand swithching seq1 and seq2\n", 2);
                        dual_print($log, "full seq1: ".$subseq2_prot_obj_full->seq."\n", 2);
                        dual_print($log, "subseq1: ".$subseq2_prot_obj->seq."\n", 2);
                        dual_print($log, "full seq2: ".$subseq1_prot_obj_full->seq."\n", 2);
                        dual_print($log, "subseq2: ".$subseq1_prot_obj->seq."\n", 2);

			#getting part1 for merging
			if( exists $all_db_fasta_IDs{lc($gene_feature1->seq_id)}){
				my $seq = $db_fasta->seq($ID_correct, $gene_feature2->start+1, $gene_feature2->end);

				#create seq test1 by merging with part2
				my $seq1 = $db_fasta->seq($ID_correct, $gene_feature1->start, $gene_feature2->start-1)."N".$seq;
				my $cds_obj1 = Bio::Seq->new(-seq => $seq1, -alphabet => 'dna' );
				$cds_obj1 = $cds_obj1->revcom();
				my $prot_obj1 = $cds_obj1->translate(-codontable_id => $codonTable) ;
                                dual_print($log, "Test frame 1 : ".$prot_obj1->seq."\n", 2);
                                if (index($prot_obj1->seq, $prot_second->seq) != -1) {
                                        dual_print($log, "frame 1 contains protein2\nIs the codon stop from the first gene real? Does the codon code for a stop codon? Is there a substition? Is it a pseudogene?\n");
					$found++;
					$case_frame1++;
					$merged++ if ($pseudo);
					push @{$obj_case->{frame}}, 1 ;
				}

				#create seq test2 by merging with part2
				my $seq2 = $db_fasta->seq($ID_correct, $gene_feature1->start, $gene_feature2->start-1)."NN".$seq;
				my $cds_obj2 = Bio::Seq->new(-seq => $seq2, -alphabet => 'dna' );
				$cds_obj2 = $cds_obj2->revcom();
				my $prot_obj2 = $cds_obj2->translate(-codontable_id => $codonTable) ;
                                dual_print($log, "Test frame 2 : ".$prot_obj2->seq."\n", 2);
                                if (index($prot_obj2->seq, $prot_second->seq) != -1) {
                                        dual_print($log, "frame 2 contains protein2\n");
					push @{$obj_case->{insert_size}}, 1;
					push @{$obj_case->{total_insert_before}}, $total_insert_before ;
					$total_insert_before += 1;
					push @{$obj_case->{insert_position}}, $gene_feature2->start;
					push @{$obj_case->{frame}}, 2 ;
					$found++;
					$case_frame2++;
					$merged++;
				}

				#create seq test3  by merging with part2
				my $seq3 = $db_fasta->seq($ID_correct, $gene_feature1->start, $gene_feature2->start-1)."NNN".$seq;
				my $cds_obj3 = Bio::Seq->new(-seq => $seq3, -alphabet => 'dna' );
				$cds_obj3 = $cds_obj3->revcom();
				my $prot_obj3 = $cds_obj3->translate(-codontable_id => $codonTable) ;
                                dual_print($log, "Test frame 3 : ".$prot_obj3->seq."\n", 2);
                                if (index($prot_obj3->seq, $prot_second->seq) != -1) {
                                        dual_print($log, "frame 3 contains protein2\n");
					push @{$obj_case->{insert_size}}, 2;
					push @{$obj_case->{total_insert_before}}, $total_insert_before;
					$total_insert_before += 2;
					push @{$obj_case->{insert_position}}, $gene_feature2->start;
					push @{$obj_case->{frame}}, 3 ;
					$found++;
					$case_frame3++;
					$merged++;
				}
			}
			else{
                                dual_print($log, "ERROR ".$gene_feature1->seq_id." not found among the db!\n");
			}
			if (! $found){
                                my $msg = "subseq2 not found in any frame. There is a problem when preparing subseq2\n";
                                dual_warn($log, $msg, 3);
			}
			elsif($found>1){
                                dual_print($log, "interesting, subseq2 found in $found frames\n");
			}
		}
		# ---- strand + ----
		else{
			my $found = 0;
			my $prot_second=$subseq2_prot_obj;
                        dual_print($log, "full seq1: ".$subseq1_prot_obj_full->seq."\n", 2);
                        dual_print($log, "subseq1: ".$subseq1_prot_obj->seq."\n", 2);
                        dual_print($log, "full seq2: ".$subseq2_prot_obj_full->seq."\n", 2);
                        dual_print($log, "subseq2: ".$subseq2_prot_obj->seq."\n", 2);

			#getting part1 for merging
			if( exists $all_db_fasta_IDs{lc($gene_feature1->seq_id)}){
				my $seq = $db_fasta->seq($ID_correct, $gene_feature1->start, $gene_feature1->end-1);

				#create seq test1 by merging with part2
				my $seq1 = $seq."N".$db_fasta->seq($ID_correct, $gene_feature1->end + 1, $gene_feature2->end);
				my $cds_obj1 = Bio::Seq->new(-seq => $seq1, -alphabet => 'dna' );
				my $prot_obj1 = $cds_obj1->translate(-codontable_id => $codonTable) ;
                                dual_print($log, "Test frame 1 : ".$prot_obj1->seq."\n", 2);
                                if (index($prot_obj1->seq, $prot_second->seq) != -1) {
                                        dual_print($log, "frame 1 contains protein2\nIs the codon stop from the first gene real? Does the codon code for a stop codon? Is there a substition? Is it a pseudogene?\n");
					$found++;
					$case_frame1++;
					$merged++ if ($pseudo);
					push @{$obj_case->{frame}}, 1 ;
				}

				#create seq test2 by merging with part2
				my $seq2 = $seq."NN".$db_fasta->seq($ID_correct, $gene_feature1->end + 1, $gene_feature2->end);
				my $cds_obj2 = Bio::Seq->new(-seq => $seq2, -alphabet => 'dna' );
				my $prot_obj2 = $cds_obj2->translate(-codontable_id => $codonTable) ;
                                dual_print($log, "frame 2 : ".$prot_obj2->seq."\n", 2);
                                if (index($prot_obj2->seq, $prot_second->seq) != -1) {
                                        dual_print($log, "frame 2 contains protein2\n");
					push @{$obj_case->{insert_size}}, 1;
					push @{$obj_case->{total_insert_before}}, $total_insert_before;
					$total_insert_before += 1;
					push @{$obj_case->{insert_position}}, $gene_feature1->end;
					push @{$obj_case->{frame}}, 2 ;
					$found++;
					$case_frame2++;
					$merged++;
				}

				#create seq test3  by merging with part2
				my $seq3 = $seq."NNN".$db_fasta->seq($ID_correct, $gene_feature1->end + 1, $gene_feature2->end);
				my $cds_obj3 = Bio::Seq->new(-seq => $seq3, -alphabet => 'dna' );
				my $prot_obj3 = $cds_obj3->translate(-codontable_id => $codonTable) ;
                                dual_print($log, "frame 3 : ".$prot_obj3->seq."\n", 2);
                                if (index($prot_obj3->seq, $prot_second->seq) != -1) {
                                        dual_print($log, "frame 3 contains protein2\n");
					push @{$obj_case->{insert_size}}, 2;
					push @{$obj_case->{total_insert_before}}, $total_insert_before;
					$total_insert_before += 2;
					push @{$obj_case->{insert_position}}, $gene_feature1->end;
					push @{$obj_case->{frame}}, 3 ;
					$found++;
					$case_frame3++;
					$merged++;
				}
				if (! $found){
                                my $msg = "subseq2 not found in any frame. There is a problem when preparing subseq2\n";
                                dual_warn($log, $msg, 3);
				}
				elsif($found>1){
                                        dual_print($log, "interesting, subseq2 found in $found frames\n");
				}
			}
			else{
                                dual_print($log, "ERROR ".$gene_feature1->seq_id." not found among the db!\n");
			}
		}
	}
	return $merged;
}

#return AA size of the overlaping region between the genes if any
sub get_overlap {
	my ($obj_case)=@_;

	my $overlap = 0;
	my $list_gene_feature = $obj_case->{list_gene_feature};
	my @listGene =  @$list_gene_feature;

	while (@listGene > 1){
		my $gene_feature1 = shift @listGene;
		my $gene_feature2 = $listGene[0];

		if($gene_feature1->end > $gene_feature2->start){
			my $overlap_region = $gene_feature1->end - $gene_feature2->start + 1;
                        dual_print($log, "overlap_region nt = $overlap_region\n", 2);
			$overlap += ($overlap_region * 2); # overlap touch both gene
		}
		else{
                        dual_print($log, "No overlap_region region!\n", 2);
			$overlap+=0;
		}
	}
	return int($overlap/3);
}

# check size of protein use for ithe inference and create an average value
sub retrieve_expected_protein_length{
	my ($obj_case)=@_;
	my $expected_protein_length;
	my $total_size = undef;
	my $nb_prot = undef;
	my $gene_list = $obj_case->{list_gene_feature};

	#fill object information
	foreach my $gene_feature (@{$gene_list}) {
		#print $gene_feature->gff_string()."\n";
		my $obj_sub_gene = sub_gene->new(id => $gene_feature->_tag_value('ID'));
		my $gene_size = ($gene_feature->end - $gene_feature->start + 1) / 3;
		$obj_sub_gene->{current_dna_length} = $gene_feature->end - $gene_feature->start + 1;
		$obj_sub_gene->{current_aa_length} = $gene_size;

                dual_print($log, $gene_feature->_tag_value('Name')." has a AA size of: $gene_size\n");

		if ($gene_feature->has_tag("inference")){
			my @inference_atts = $gene_feature->get_tag_values("inference");
			foreach my $inference_att (@inference_atts){
				#similar for UniProtKB and motif for HAMAP
				if ( $inference_att =~ /^similar/ or $inference_att =~ /^protein motif/){
					my @data = split /:/, $inference_att ;
					my $inference_value = $data[$#data];
					my $inference_db = $data[$#data-1];
					$obj_sub_gene->{inference_db} = $inference_db;
					$obj_sub_gene->{inference_value} = $inference_value;
				}
			}
		}
		else{
                        my $msg = "No inference attribute found\n";
                        dual_warn($log, $msg, 3);
		}
		$obj_case->{hash_sub_gene_obj}{$obj_sub_gene->{id}} = $obj_sub_gene;
	}

	#get original protein length
	foreach my $id_obj_sub_gene ( keys %{$obj_case->{hash_sub_gene_obj}}){
		my $obj_sub_gene = $obj_case->{hash_sub_gene_obj}{$id_obj_sub_gene};
		my $prot_name = $obj_sub_gene->{inference_value};

		if (lc($obj_sub_gene->{inference_db}) eq "uniprotkb"){
                        dual_print($log, "Inference made with Uniprot, looking for protein size: $prot_name\n");
			if( exists $all_db_db_IDs{lc($prot_name)}){
				my $protID_correct = $all_db_db_IDs{lc($prot_name)};
				$obj_sub_gene->{inference_aa_length} = $db_db->length( $protID_correct );
				$obj_sub_gene->{inference_aa_seq} = $db_db->seq($protID_correct);
			}
			else{
                                dual_print($log, "ERROR $prot_name not found among the db!\n");
			}
		}
		elsif (lc($obj_sub_gene->{inference_db}) eq "hamap"){
                        dual_print($log, "Inference made with HAMAP, looking for protein size using internet: $prot_name\n");
			fetcher_HAMAP($obj_sub_gene) if (! $skip_hamap);
		}
		else{
                        dual_print($log, "Inference made with ".$obj_sub_gene->{inference_db}.", not yet implemented\n");
		}

		#I have a length, let's check if we can merge the case
		if ( $obj_sub_gene->{inference_aa_length} ){
			push @{$obj_case->{list_aa_size}}, $obj_sub_gene->{inference_aa_length};
			$total_size += $obj_sub_gene->{inference_aa_length};
			$nb_prot++;
                        dual_print($log, "AA length found ".$obj_sub_gene->{inference_aa_length}." \n");
		}
	}
	$obj_case->{expected_length} = int($total_size/$nb_prot) if $nb_prot;
}

sub fetcher_HAMAP {
		my ($obj_sub_gene) = @_;

		my $id = $obj_sub_gene->{inference_value};
		my $ua = LWP::UserAgent->new;
		$ua->timeout(10);
		$ua->env_proxy;

		my $size = undef;
		my $url="https://hamap.expasy.org/rule/".$id;
		#	print $url."\n";
		my $response = $ua->get($url);
		if ($response->is_success) {
				my $string = $response->decoded_content;
				#print $string."\n";
				if ( $string =~ /.*<td>(.*)amino acids<\/td>.*/){
                                        dual_print($log, "size range: $1\n");
					my @data = split /-/, $1 ;

					if($hamap_size eq "low"){
						$size = $data[0];
					}
					elsif($hamap_size eq "middle"){
						my $nb_values = scalar @data;
						foreach my $value (@data){
							$size+=$value
						}
						$size = int($size/$nb_values);
					}
					elsif($hamap_size eq "high"){
						$size = $data[$#data];
					}

				}
				else{
                                        dual_print($log, "No size found for $id\n");
				}
			$obj_sub_gene->{inference_aa_length} = $size;
		}
		else {
                        dual_print($log, $response->status_line."\n");
                        die;
		}
}


__END__

=head1 NAME

agat_sp_prokka_fragmented_gene_annotations.pl

=head1 DESCRIPTION

The script aims to look at fragmented gene annotations (FRAGS) within prokka annotations.
The FRAGS represent two (or more) ORFs that are in close proximity and are annotated
with homology to the same gene. In such cases, Prokka ads an _n suffix to the gene ID.
For example, a splitted genX can then be found as genX_1 and genX_2 in the GFF.
See here for a case: https://github.com/tseemann/prokka/issues/502

* The script will inform you how many case there is in your annotation.
* If you think the FRAGS is due to a sequencing error (frameshift due to short indel),
using the --frags parameter will fix the FRAGS if genX_1 and genX_2 are not in the same frame.
The gff and the fasta file will be modified. The gene are merged, an insertion of
one or two N will be added in between the genes to fix the frameshift.
* If you think the FRAGS is not due to a sequencing error, use the --pseudo parameter,
the gff will be fix (gene merged) and the agat_pseudo attribute (the value is the position of the codon stop)
will be added to the related features.
* using --frags and --pseudo is similar to use only --frags, except when no frameshift
is found for a detected FRAGS (both gene are in the same frame), the agat_pseudo
attribute is also added to the related features.

How the tool detecte the FRAGS?
* Search for cases where contiguous genes have the same name (e.g. lpxA_1 lpxA_2).
* If so we look at the size of the protein of each of those genes (lpxA_1 AA=175 ; lpxA_2 AA=116),
and compute the size when merged togeter (devoided of the overlap if any) => here 270 AA
* Then we look at the size of the protein used to infer the name (lpxA_1 inferred from Q9PIM1 = 263 AA ; lpxA_2 inferred from P0A722 = 262 AA )
and compute the average length of the reference protein: here 262AA. We add 20% to the length to be sure to include border cases => 282AA.
* Compare the length of the merged proteins (262 AA) against the reference protein length (282).
If the the expected protein length (282 AA) is longer we have a FRAGS.


=head1 SYNOPSIS

    agat_sp_prokka_fragmented_gene_annotations.pl -gff infile.gff --fasta genome.fa --db prokka/prokka_bacteria_sprot.fa  -o outfolder
    agat_sp_prokka_fragmented_gene_annotations.pl --help

=head1 OPTIONS

=over 8

=item B<--gff>

Input genome GTF/GFF file. Mandatory.

=item B<-f>, B<--fa> or B<--fasta>

Input genome fasta file. Mandatory.

=item B<--db>

Input Uniprot fasta file used by prokka. Mandatory.

=item B<--frags>

Merge and fix detected FRAGS if not in the same frame

=item B<--pseudo>

Merge detected FRAGS and add the agat_pseudo attribute (value will be the location of the first stop codon met).

=item B<--hamap_size>

Some protein function are not infered by Uniprot but by Hamap. In such case the information
is retrieved from the web. As hamap provide a family profile, the protein size if a range.
"low" option will use the low value of the range,
"middle" option will use the average of the range,
"high" option will the the high value of the range.
Default "high".

=item B<--ct>, B<--codon> or B<--table>

Codon table to use. [default 1]

=item B<--skip_hamap>

For test purpose it could be useful to skip hamap, because it requires fetching information through internet.

=item B<-o> , B<--output> or B<--out>

Output folder. Mandatory.

=item B<-v>

verbose mode. Default off.

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
