#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use File::Path qw(remove_tree); # to remove directory easily (tmp directory)

=head1 DESCRIPTION

Test to verify the output of scripts using gff_syntax cases

=cut

################################################################################
#   set number of test according to number of scripts
my $nb_test;
BEGIN{
  open(FH, '<', __FILE__) or die $!;
  while(<FH>){
    if( $_ !~ /^#/){ #if it is not a comment line
      if (index($_ , "ok(") != -1) { # This is a test line. We should avoid counting this line. It is why we remove 1 later
        $nb_test++;
      }
    }
  }
  $nb_test--;
  close(FH);
}
#
################################################################################

use Test::More tests => $nb_test ;

# Check if has to be run in Devel::Cover or not
my $script_prefix="";
if (exists $ENV{'HARNESS_PERL_SWITCHES'} ) {
  if ($ENV{'HARNESS_PERL_SWITCHES'} =~ m/Devel::Cover/) {
    $script_prefix="perl -MDevel::Cover ";
  }
}

# shared variables
my $input_folder = "t/scripts_output/in";
my $output_folder = "t/scripts_output/out";
my $outprefix = "tmp";
my $outtmp = "$outprefix.gff"; # path file where to save temporary output
my $config="agat_config.yaml";
my $script;
my $result;
my $result2;


# remove config in local folder if exists and potential tmp file already existing
cleaning_log();

# -------------------------- check agat_convert_bed2gff -------------------------

$script = $script_prefix."bin/agat_convert_bed2gff.pl";
$result = "$output_folder/agat_convert_bed2gff_1.gff";
system(" $script --bed $input_folder/test.bed -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
cleaning_log("test.bed");

# -------------------------- check agat_convert_embl2gff -------------------------

$script = $script_prefix."bin/agat_convert_embl2gff.pl";
$result = "$output_folder/agat_convert_embl2gff_1.gff";
system(" $script --embl $input_folder/agat_convert_embl2gff_1.embl -o $outtmp --emblmygff3 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
cleaning_log("agat_convert_embl2gff_1.embl");

# -------------------------- check agat_convert_genscan2gff -------------------------

$script = $script_prefix."bin/agat_convert_genscan2gff.pl";
$result = "$output_folder/agat_convert_genscan2gff_1.gff";
system(" $script --genscan $input_folder/test.genscan -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
cleaning_log("test.genscan");

# -------------------------- check agat_convert_mfannot2gff -------------------------

$script = $script_prefix."bin/agat_convert_mfannot2gff.pl";
$result = "$output_folder/agat_convert_mfannot2gff_1.gff";
system(" $script --mfannot $input_folder/test.mfannot -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
cleaning_log("test.mfannot");

$script = $script_prefix."bin/agat_convert_mfannot2gff.pl";
$result = "$output_folder/agat_convert_mfannot2gff_2.gff";
system(" $script --mfannot $input_folder/test.mfannot2 -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
cleaning_log();

# -------------------------- check agat_convert_minimap2_bam2gff.pl --------------------------
$script = $script_prefix."bin/agat_convert_minimap2_bam2gff.pl";
$result = "$output_folder/agat_convert_minimap2_bam2gff_1.gff";
system(" $script -i $input_folder/test_minimap2.sam -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
cleaning_log("test_minimap2.sam");

# -------------------------- check agat_convert_sp_gff2bed -------------------------

$script = $script_prefix."bin/agat_convert_sp_gff2bed.pl";
$result = "$output_folder/agat_convert_sp_gff2bed_1.gff";
system(" $script --gff $input_folder/1.gff -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
cleaning_log("1.gff");

# -------------------------- check agat_convert_sp_gff2gtf -------------------------
my $convert_sp_gff2gtf_folder = "$input_folder/agat_convert_sp_gff2gtf";

$script = $script_prefix."bin/agat_convert_sp_gff2gtf.pl";
$result = "$convert_sp_gff2gtf_folder/agat_convert_sp_gff2gtf_1.gtf";
system(" $script --gff $input_folder/1.gff --gtf_version 3 -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
cleaning_log("1.gff");

$result = "$convert_sp_gff2gtf_folder/agat_convert_sp_gff2gtf_2.gtf";
system(" $script --gff $convert_sp_gff2gtf_folder/stop_start_an_exon.gff --gtf_version 3 -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
cleaning_log("stop_start_an_exon.gff");

$result = "$convert_sp_gff2gtf_folder/agat_convert_sp_gff2gtf_3.gtf";
system(" $script --gff $convert_sp_gff2gtf_folder/stop_split_over_two_exons.gff --gtf_version 3 -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
cleaning_log("stop_split_over_two_exons.gff");

$result = "$convert_sp_gff2gtf_folder/result_issue_245.gtf";
system(" $script --gff $convert_sp_gff2gtf_folder/issue_245.gff --gtf_version 3 -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
cleaning_log("issue_245.gff");

# --------check agat_sp_to_tabulated.pl-------------

$script = $script_prefix."bin/agat_convert_sp_gff2tsv.pl";
$result = "$output_folder/agat_convert_sp_gff2tsv_1.tsv";
system(" $script --gff $input_folder/1.gff -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
cleaning_log("1.gff");

# -------------------------- check agat_convert_sp_gff2zff -------------------------

$script = $script_prefix."bin/agat_convert_sp_gff2zff.pl";
$result = "$output_folder/agat_convert_sp_gff2zff_1.gff";
system(" $script --gff $input_folder/1.gff --fasta $input_folder/1.fa -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outprefix.ann") == 0, "output $script");
cleaning_log("1.gff");

# -------------------------- check agat_convert_sp_gxf2gxf.pl -------------------------

# XXX No need to be tested, it is tested by gff_syntax tests

# -------------------------- check agat_sp_add_attribute_shortest_exon_size -------------------------

$script = $script_prefix."bin/agat_sp_add_attribute_shortest_exon_size.pl";
$result = "$output_folder/agat_sp_add_attribute_shortest_exon_size.gff";
system(" $script --gff $input_folder/1.gff -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
cleaning_log("1.gff");

# -------------------------- check agat_sp_add_attribute_shortest_intron_size -------------------------

$script = $script_prefix."bin/agat_sp_add_attribute_shortest_intron_size.pl";
$result = "$output_folder/agat_sp_add_attribute_shortest_intron_size.gff";
system(" $script --gff $input_folder/1.gff -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
cleaning_log("1.gff");

# -------------------------- check agat_sp_add_intergenic_regions -------------------------

$script = $script_prefix."bin/agat_sp_add_intergenic_regions.pl";
$result = "$output_folder/agat_sp_add_intergenic_regions_1.gff";
system(" $script --gff $input_folder/1.gff -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
cleaning_log("1.gff");

# -------------------------- check agat_sp_add_introns -------------------------

$script = $script_prefix."bin/agat_sp_add_introns.pl";
$result = "$output_folder/agat_sp_add_introns_1.gff";
system(" $script --gff $input_folder/1.gff -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
cleaning_log("1.gff");

# -------------------------- check agat_sp_add_splice_sites -------------------------

$script = $script_prefix."bin/agat_sp_add_splice_sites.pl";
$result = "$output_folder/agat_sp_add_splice_sites_1.gff";
system(" $script --gff $input_folder/0.gff -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
cleaning_log("0.gff");

# ---------------------- check agat_sp_add_start_and_stop ----------------------

$script = $script_prefix."bin/agat_sp_add_start_and_stop.pl";
$result = "$output_folder/agat_sp_add_start_and_stop_1.gff";
system(" $script --gff $input_folder/agat_sp_add_start_and_stop.gff --fasta $input_folder/1.fa --ni -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
cleaning_log("agat_sp_add_start_and_stop.gff");

$script = $script_prefix."bin/agat_sp_add_start_and_stop.pl";
$result = "$output_folder/agat_sp_add_start_and_stop_2.gff";
system(" $script --gff $input_folder/agat_sp_add_start_and_stop.gff --fasta $input_folder/1.fa -e --ni -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
cleaning_log("agat_sp_add_start_and_stop.gff");

# -------------------------- check agat_sp_statistics --------------------------

$script = $script_prefix."bin/agat_sp_alignment_output_style.pl";
$result = "$output_folder/agat_sp_alignment_output_style_1.gff";
system(" $script --gff $input_folder/0.gff -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
cleaning_log("0.gff");

# --------check agat_sp_clipN_seqExtremities_and_fixCoordinates.pl-------------

# XXX

# ------------------- check agat_sp_compare_two_annotations script-------------------
$script = $script_prefix."bin/agat_sp_compare_two_annotations.pl";
$result = "$output_folder/agat_sp_compare_two_annotations_1.txt";
system(" $script --gff1 $input_folder/1.gff  --gff2 $input_folder/1.gff -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp/report.txt") == 0, "output $script");
cleaning_log("1.gff");

$script = $script_prefix."bin/agat_sp_compare_two_annotations.pl";
$result = "$output_folder/agat_sp_compare_two_annotations_2.txt";
system(" $script --gff1 $input_folder/agat_sp_compare_two_annotations/file1.gff  --gff2 $input_folder/agat_sp_compare_two_annotations/file2.gff -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp/report.txt") == 0, "output $script");
cleaning_log("file1.gff");

$script = $script_prefix."bin/agat_sp_compare_two_annotations.pl";
$result = "$output_folder/agat_sp_compare_two_annotations_3.txt";
system(" $script --gff1 $input_folder/agat_sp_compare_two_annotations/file2.gff  --gff2 $input_folder/agat_sp_compare_two_annotations/file1.gff -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp/report.txt") == 0, "output $script");
cleaning_log("file2.gff");

# --------check agat_sp_compare_two_BUSCOs.pl -------------

# XXX

# ------------------- check agat_sp_complement_annotations script-------------------
$script = $script_prefix."bin/agat_sp_complement_annotations.pl";
$result = "$output_folder/agat_sp_complement_annotations_1.gff";
system(" $script --ref t/gff_syntax/in/25_test.gff  --add t/gff_syntax/in/9_test.gff -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
cleaning_log("25_test.gff");

$result = "$output_folder/agat_sp_complement_annotations_2.gff";
system(" $script --ref $input_folder/agat_sp_complement_annotations/agat_sp_complement_annotations_ref.gff  --add $input_folder/agat_sp_complement_annotations/agat_sp_complement_annotations_add.gff -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
cleaning_log("agat_sp_complement_annotations_ref.gff");

# --------check agat_sp_ensembl_output_style.pl-------------

$script = $script_prefix."bin/agat_sp_ensembl_output_style.pl";
$result = "$output_folder/agat_sp_ensembl_output_style_1.gff";
system("$script --gff $input_folder/0.gff -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
cleaning_log("0.gff");

# --------check agat_sp_extract_attributes.pl-------------

$script = $script_prefix."bin/agat_sp_extract_attributes.pl";
$result = "$output_folder/agat_sp_extract_attributes_1.txt";
system(" $script --gff $input_folder/1.gff --att protein_id -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result ${outprefix}_protein_id.gff") == 0, "output $script");
cleaning_log("1.gff");

# --------check agat_sp_extract_sequences.pl-------------

$script = $script_prefix."bin/agat_sp_extract_sequences.pl";
$result = "$input_folder/agat_sp_extract_sequences/agat_sp_extract_sequences_1.fa";
system(" $script --gff $input_folder/1.gff --fasta $input_folder/1.fa -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script test1");
cleaning_log("1.gff");

$script = $script_prefix."bin/agat_sp_extract_sequences.pl";
$result = "$input_folder/agat_sp_extract_sequences/agat_sp_extract_sequences_split.fa";
system(" $script --gff $input_folder/1.gff --fasta $input_folder/1.fa --split -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script test2");
cleaning_log("1.gff");

$script = $script_prefix."bin/agat_sp_extract_sequences.pl";
$result = "$input_folder/agat_sp_extract_sequences/agat_sp_extract_sequences_merge.fa";
system(" $script --gff $input_folder/1.gff --fasta $input_folder/1.fa -t exon --merge -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script test3");
cleaning_log("1.gff");

$script = $script_prefix."bin/agat_sp_extract_sequences.pl";
$result = "$input_folder/agat_sp_extract_sequences/agat_sp_extract_sequences_full.fa";
system(" $script --gff $input_folder/1.gff --fasta $input_folder/1.fa --full -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script test4");
cleaning_log("1.gff");

$script = $script_prefix."bin/agat_sp_extract_sequences.pl";
$result = "$input_folder/agat_sp_extract_sequences/agat_sp_extract_sequences_attributes_kept.fa";
system(" $script --gff $input_folder/1.gff --fasta $input_folder/1.fa --keep_attributes -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script test5");
cleaning_log("1.gff");

$script = $script_prefix."bin/agat_sp_extract_sequences.pl";
$result = "$input_folder/agat_sp_extract_sequences/agat_sp_extract_sequences_parent_attributes_kept.fa";
system(" $script --gff $input_folder/1.gff --fasta $input_folder/1.fa --keep_parent_attributes -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script test6");
cleaning_log("1.gff");

# --------check agat_sp_filter_by_locus_distance.pl-------------
$script = $script_prefix."bin/agat_sp_filter_by_locus_distance.pl";
$result = "$output_folder/agat_sp_filter_by_locus_distance_1.gff";
system(" $script --gff $input_folder/1.gff -o $outtmp 2>&1 1>/dev/null");
ok( system("diff $result $outtmp") == 0, "output $script");
cleaning_log("1.gff");

# --------check agat_sp_filter_by_mrnaBlastValue.pl-------------

# XXX

# --------check agat_sp_filter_by_ORF_size.pl-------------
$script = $script_prefix."bin/agat_sp_filter_by_ORF_size.pl";
$result = "$output_folder/agat_sp_filter_by_ORF_size_sup100.gff";
system(" $script --gff $input_folder/1.gff -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result ${outprefix}_sup100.gff") == 0, "output $script");
cleaning_log("1.gff");

# --------check agat_sp_filter_feature_by_attribute_presence.pl-------------
$script = $script_prefix."bin/agat_sp_filter_feature_by_attribute_presence.pl";
$result = "$output_folder/agat_sp_filter_feature_by_attribute_presence_1.gff";
system(" $script --gff $input_folder/1.gff -o $outtmp -a protein_id 2>&1 1>/dev/null");
ok( system("diff $result $outtmp") == 0, "output $script");
cleaning_log("1.gff");

# --------check agat_sp_filter_feature_by_attribute_value.pl-------------

$script = $script_prefix."bin/agat_sp_filter_feature_by_attribute_value.pl";
$result = "$output_folder/agat_sp_filter_feature_by_attribute_value_1.gff";
system(" $script --gff $input_folder/1.gff -o $outtmp --value Os01t0100100-01 -p level3 -a protein_id 2>&1 1>/dev/null");
ok( system("diff $result $outtmp") == 0, "output $script");
cleaning_log("1.gff");

# --------check agat_sp_filter_feature_from_keep_list.pl-------------

$script = $script_prefix."bin/agat_sp_filter_feature_from_keep_list.pl";
$result = "$output_folder/agat_sp_filter_feature_from_keep_list_1.gff";
system(" $script --gff $input_folder/1.gff -o $outtmp --kl $input_folder/keep_list.txt 2>&1 1>/dev/null");
ok( system("diff $result $outtmp") == 0, "output $script");
cleaning_log("1.gff");

# --------check agat_sp_filter_feature_from_kill_list.pl-------------

$script = $script_prefix."bin/agat_sp_filter_feature_from_kill_list.pl";
$result = "$output_folder/agat_sp_filter_feature_from_kill_list_1.gff";
system(" $script --gff $input_folder/1.gff -o $outtmp --kl $input_folder/kill_list.txt 2>&1 1>/dev/null");
ok( system("diff $result $outtmp") == 0, "output $script");
cleaning_log("1.gff");

# --------check agat_sp_filter_gene_by_intron_numbers.pl-------------

$script = $script_prefix."bin/agat_sp_filter_gene_by_intron_numbers.pl";
$result = "$output_folder/agat_sp_filter_gene_by_intron_numbers_1.gff";
system(" $script --gff $input_folder/1.gff -o $outtmp 2>&1 1>/dev/null");
ok( system("diff $result $outtmp") == 0, "output $script");
cleaning_log("1.gff");

# --------check agat_sp_filter_gene_by_length.pl-------------

$script = $script_prefix."bin/agat_sp_filter_gene_by_length.pl";
$result = "$output_folder/agat_sp_filter_gene_by_length_1.gff";
system(" $script --gff $input_folder/1.gff --size 1000 --test \"<\" -o $outtmp 2>&1 1>/dev/null");
ok( system("diff $result $outtmp") == 0, "output $script");
cleaning_log("1.gff");

# --------check agat_sp_filter_incomplete_gene_coding_models.pl-------------

$script = $script_prefix."bin/agat_sp_filter_incomplete_gene_coding_models.pl";
$result = "$output_folder/agat_sp_filter_incomplete_gene_coding_models_1.gff";
$result2 = "$output_folder/agat_sp_filter_incomplete_gene_coding_models_incomplete_1.gff";
system(" $script --gff $input_folder/1.gff --fasta $input_folder/1.fa -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
ok( system("diff $result2 $outprefix"."_incomplete.gff") == 0, "output $script");
cleaning_log("1.gff");

# --------check agat_sp_filter_record_by_coordinates.pl-------------

$script = $script_prefix."bin/agat_sp_filter_record_by_coordinates.pl";
$result = "$output_folder/agat_sp_filter_record_by_coordinates_1.gff";
system(" $script --gff $input_folder/1.gff --tsv $input_folder/coordinates.tsv -o $outtmp 2>&1 1>/dev/null");
ok( system("diff $result $outtmp/remaining.gff3") == 0, "output $script");
cleaning_log("1.gff");

# --------check agat_sp_fix_cds_phases.pl-------------

$script = $script_prefix."bin/agat_sp_fix_cds_phases.pl";
$result = "$output_folder/agat_sp_fix_cds_phases_1.gff";
system(" $script --gff $input_folder/1.gff --fasta $input_folder/1.fa -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
cleaning_log("1.gff");

# --------check agat_sp_fix_features_locations_duplicated.pl-------------
# removed because order can change. So not reproducible at 100%

$script = $script_prefix."bin/agat_sp_fix_features_locations_duplicated.pl";
$result = "$output_folder/agat_sp_fix_features_locations_duplicated_1.gff";
system(" $script --gff $input_folder/agat_sp_fix_features_locations_duplicated/test.gff -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
cleaning_log("test.gff");

# --------check agat_sp_fix_fusion.pl-------------

$script = $script_prefix."bin/agat_sp_fix_fusion.pl";
$result = "$output_folder/agat_sp_fix_fusion_1.txt";
system(" $script --gff $input_folder/1.gff --fasta $input_folder/1.fa -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff -b -I '^Job done in' -I '^Job done in' $result $outprefix-report.txt") == 0, "output $script");
cleaning_log("1.gff");

# --------agat_sp_fix_longest_ORF.pl-------------

$script = $script_prefix."bin/agat_sp_fix_longest_ORF.pl";
$result = "$output_folder/agat_sp_fix_longest_ORF_1.txt";
system(" $script --gff $input_folder/1.gff --fasta $input_folder/1.fa -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff -b -I '^Job done in' -I '^Job done in' $result $outprefix-report.txt") == 0, "output $script");
cleaning_log("1.gff");

# --------check agat_sp_fix_overlaping_genes.pl-------------

$script = $script_prefix."bin/agat_sp_fix_overlaping_genes.pl";
$result = "$output_folder/agat_sp_fix_overlaping_genes_1.gff";
system(" $script --gff $input_folder/genes_overlap.gff -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
cleaning_log("genes_overlap.gff");

$script = $script_prefix."bin/agat_sp_fix_overlaping_genes.pl";
$result = "$output_folder/agat_sp_fix_overlaping_genes_2.gff";
system(" $script --gff $input_folder/genes_overlap.gff -o $outtmp --merge 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
cleaning_log("genes_overlap.gff");

# --------check agat_sp_fix_small_exon_from_extremities.pl-------------

# XXX

# --------check agat_sp_flag_short_introns.pl-------------

# XXX

# --------check agat_sp_flag_premature_stop_codons.pl-------------
# I use result from another test because it shifted the annotation location, that allows to create pseudogenes because I use the original fasta not shifted
$script = $script_prefix."bin/agat_sp_flag_premature_stop_codons.pl";
$result = "$output_folder/agat_sp_flag_premature_stop_codons_1.gff";
system(" $script --gff $input_folder/prokka_fragmented_genes.gff --fasta $input_folder/prokka_cav_10DC88.fa -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
cleaning_log("prokka_fragmented_genes.gff");

# --------check agat_sp_functional_statistics.pl-------------

$script = $script_prefix."bin/agat_sp_functional_statistics.pl";
$result = "$output_folder/agat_sp_functional_statistics_1.txt";
system(" $script --gff $input_folder/function.gff -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $output_folder/agat_sp_functional_statistics/table_gene_mrna.txt $outtmp/gene\@mrna/table_per_feature_type.txt") == 0, "output $script");
ok( system("diff $output_folder/agat_sp_functional_statistics/table_repeat.txt $outtmp/repeat_region/table_per_feature_type.txt") == 0, "output $script");
cleaning_log("function.gff");

# --------check agat_sp_keep_longest_isoform.pl-------------

$script = $script_prefix."bin/agat_sp_keep_longest_isoform.pl";
$result = "$output_folder/agat_sp_keep_longest_isoform_1.gff";
system(" $script --gff $input_folder/1.gff -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
cleaning_log("1.gff");

# --------check agat_sp_kraken_assess_liftover.pl-------------

$script = $script_prefix."bin/agat_sp_kraken_assess_liftover.pl";
$result = "$output_folder/agat_sp_kraken_assess_liftover_1.gff";
system(" $script --gtf $input_folder/test_kraken.gtf -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
cleaning_log("test_kraken.gff");

# --------check agat_sp_list_short_introns.pl-------------

$script = $script_prefix."bin/agat_sp_list_short_introns.pl";
$result = "$output_folder/agat_sp_list_short_introns_1.txt";
system(" $script --gff $input_folder/1.gff -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
cleaning_log("1.gff");

# --------check agat_sp_load_function_from_protein_align.pl-------------

# XXX

# --------check agat_sp_manage_IDs.pl-------------

$script = $script_prefix."bin/agat_sp_manage_IDs.pl";
$result = "$output_folder/agat_sp_manage_IDs_1.gff";
system(" $script --gff $input_folder/1.gff --prefix NBIS --ensembl --tair --type_dependent -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
cleaning_log("1.gff");

# --------check agat_sp_manage_UTRs.pl-------------

$script = $script_prefix."bin/agat_sp_manage_UTRs.pl";
$result = "$output_folder/agat_sp_manage_UTRs_1.gff";
system(" $script --gff $input_folder/1.gff -b -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outprefix/1_bothSides_under5.gff") == 0, "output $script");
cleaning_log("1.gff");

# --------check agat_sp_manage_attributes.pl-------------

$script = $script_prefix."bin/agat_sp_manage_attributes.pl";
$result = "$output_folder/agat_sp_manage_attributes_1.gff";
system(" $script --gff $input_folder/1.gff --att protein_id -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
cleaning_log("1.gff");

# --------check agat_sp_manage_functional_annotation.pl-------------

$script = $script_prefix."bin/agat_sp_manage_functional_annotation.pl";
$result = "$output_folder/agat_sp_manage_functional_annotation_1.gff";
system(" $script --gff $input_folder/agat_sp_manage_functional_annotation/02413F.gff --db $input_folder/agat_sp_manage_functional_annotation/uniprot_sprot_test.fasta -b $input_folder/agat_sp_manage_functional_annotation/02413F_blast.out -i $input_folder/agat_sp_manage_functional_annotation/02413F_interpro.tsv --clean_name -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system( "diff $result $outtmp/02413F.gff" ) == 0, "output $script");
cleaning_log("02413F.gff");

# --------check agat_sp_manage_introns.pl-------------

$script = $script_prefix."bin/agat_sp_manage_introns.pl";
$result = "$output_folder/agat_sp_manage_introns_1.txt";
system(" $script --gff $input_folder/1.gff -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system( "diff  -b -I '^usage:' -I '^usage:' $result $outtmp/report.txt" ) == 0, "output $script");
cleaning_log("1.gff");

# ------------------- check agat_sp_merge_annotations script-------------------

$script = $script_prefix."bin/agat_sp_merge_annotations.pl";
$result = "$output_folder/agat_sp_merge_annotations_1.gff";
system(" $script --gff $input_folder/agat_sp_merge_annotations/file1.gff  --gff $input_folder/agat_sp_merge_annotations/file2.gff -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
cleaning_log("file1.gff");

$result = "$output_folder/agat_sp_merge_annotations_2.gff";
system(" $script --gff $input_folder/agat_sp_merge_annotations/fileA.gff  --gff $input_folder/agat_sp_merge_annotations/fileB.gff -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
cleaning_log("fileA.gff");

# ------------------- check agat_sp_move_attributes_within_records script-------------------

$script = $script_prefix."bin/agat_sp_move_attributes_within_records.pl";
$result = "$output_folder/agat_sp_move_attributes_within_records.gff";
system(" $script --gff $input_folder/agat_sp_move_attributes_within_records.gff --fp exon,CDS --fc mRNA -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
cleaning_log("agat_sp_move_attributes_within_records.gff");

# ------------------- check agat_sp_prokka_fragmented_gene_annotations script-------------------

$script = $script_prefix."bin/agat_sp_prokka_fix_fragmented_gene_annotations.pl";
$result = "$output_folder/agat_sp_prokka_fix_fragmented_gene_annotations_1.gff";
$result2 = "$output_folder/agat_sp_prokka_fix_fragmented_gene_annotations_1.fa";
system(" $script --gff $input_folder/prokka_cav_10DC88.gff --fasta $input_folder/prokka_cav_10DC88.fa --db $input_folder/prokka_bacteria_sprot.fa --skip_hamap --frags -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system( "diff $result $outtmp/prokka_cav_10DC88.gff" ) == 0, "output $script");
ok( system( "diff $result2 $outtmp/prokka_cav_10DC88.fa" ) == 0, "output $script");
cleaning_log("prokka_cav_10DC88.gff");

# ------------------- check agat_sp_sensitivity_specificity script-------------------

$script = $script_prefix."bin/agat_sp_sensitivity_specificity.pl";
$result = "$output_folder/agat_sp_sensitivity_specificity_1.txt";
system(" $script --gff1 $input_folder/1.gff --gff2 $input_folder/1.gff -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
cleaning_log("1_1");

$script = $script_prefix."bin/agat_sp_sensitivity_specificity.pl";
$result = "$output_folder/agat_sp_sensitivity_specificity_2.txt";
system(" $script --gff1 $input_folder/agat_sp_sensitivity_specificity/ref0.gff3 --gff2 $input_folder/agat_sp_sensitivity_specificity/query0.gff3 -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
cleaning_log("ref0.gff3");

$script = $script_prefix."bin/agat_sp_sensitivity_specificity.pl";
$result = "$output_folder/agat_sp_sensitivity_specificity_3.txt";
system(" $script --gff1 $input_folder/agat_sp_sensitivity_specificity/ref1.gff3 --gff2 $input_folder/agat_sp_sensitivity_specificity/query1.gff3 -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
cleaning_log("ref1.gff3");

$script = $script_prefix."bin/agat_sp_sensitivity_specificity.pl";
$result = "$output_folder/agat_sp_sensitivity_specificity_4.txt";
system(" $script --gff1 $input_folder/agat_sp_sensitivity_specificity/ref2.gff3 --gff2 $input_folder/agat_sp_sensitivity_specificity/query2.gff3 -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
cleaning_log("ref2.gff3");

$script = $script_prefix."bin/agat_sp_sensitivity_specificity.pl";
$result = "$output_folder/agat_sp_sensitivity_specificity_5.txt";
system(" $script --gff1 $input_folder/agat_sp_sensitivity_specificity/ref3.gff3 --gff2 $input_folder/agat_sp_sensitivity_specificity/query3.gff3 -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
cleaning_log("ref3.gff3");
# --------check agat_sp_split_by_level2_feature.pl-------------

$script = $script_prefix."bin/agat_sp_separate_by_record_type.pl";
$result = "$output_folder/agat_sp_separate_by_record_type_1.gff";
system(" $script --gff $input_folder/1.gff -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp/trna.gff") == 0, "output $script");
cleaning_log("1.gff");

# -------------------------- check agat_sp_statistics --------------------------

$script = $script_prefix."bin/agat_sp_statistics.pl";
$result = "$output_folder/agat_sp_statistics_1.txt";
system(" $script --gff $input_folder/1.gff -o $outtmp -r 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
cleaning_log("1.gff");

# --------check agat_sp_webApollo_compliant.pl-------------

$script = $script_prefix."bin/agat_sp_webApollo_compliant.pl";
$result = "$output_folder/agat_sp_webApollo_compliant_1.gff";
system(" $script --gff $input_folder/1.gff -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
cleaning_log("1.gff");

# --------check agat_sq_add_attributes_from_tsv.pl-------------

$script = $script_prefix."bin/agat_sq_add_attributes_from_tsv.pl";
$result = "$output_folder/agat_sq_add_attributes_from_tsv_1.gff";
system(" $script --gff $input_folder/agat_sq_add_attributes_from_tsv.gff --tsv $input_folder/agat_sq_add_attributes_from_tsv.tsv -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
cleaning_log("agat_sq_add_attributes_from_tsv.gff");

# --------check agat_sq_add_hash_tag.pl-------------

$script = $script_prefix."bin/agat_sq_add_hash_tag.pl";
$result = "$output_folder/agat_sq_add_hash_tag_1.gff";
system(" $script --gff $input_folder/1.gff -i 2 -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
cleaning_log("1.gff");

# --------check agat_sq_add_locus_tag.pl-------------

$script = $script_prefix."bin/agat_sq_add_locus_tag.pl";
$result = "$output_folder/agat_sq_add_locus_tag_1.gff";
system(" $script --gff $input_folder/1.gff -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
cleaning_log("1.gff");

# --------check agat_sq_keep_annotation_from_fastaSeq.pl-------------

# XXX


# --------check agat_sq_list_attributes.pl-------------

$script = $script_prefix."bin/agat_sq_list_attributes.pl";
$result = "$output_folder/agat_sq_list_attributes_1.txt";
system(" $script --gff $input_folder/1.gff -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff -b -I '^Job done in' $result $outtmp") == 0, "output $script");
cleaning_log("1.gff");

# --------check agat_sq_manage_ID.pl-------------

$script = $script_prefix."bin/agat_sq_manage_IDs.pl";
$result = "$output_folder/agat_sq_manage_IDs_1.gff";
system(" $script --gff $input_folder/1.gff -o $outtmp 2>&1 1>/dev/null");
#run test
#ok( system("diff -b -I '^usage:' -I '^usage:' $result $outtmp") == 0, "output $script");
ok( system("diff $result $outtmp") == 0, "output $script");
cleaning_log("1.gff");

# --------check agat_sq_mask.pl-------------

# XXX


# --------check agat_sq_remove_redundant_entries.pl-------------

$script = $script_prefix."bin/agat_sq_remove_redundant_entries.pl";
$result = "$output_folder/agat_sq_remove_redundant_entries_1.gff";
system(" $script --gff $input_folder/1.gff -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
cleaning_log("1.gff");

# --------check agat_sq_rename_seqid.pl-------------

$script = $script_prefix."bin/agat_sq_rename_seqid.pl";
$result = "$output_folder/agat_sq_rename_seqid_1.gff";
system(" $script --gff $input_folder/agat_sq_rename_seqid/rename_seqid.gff --tsv $input_folder/agat_sq_rename_seqid/rename_table.tsv -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
cleaning_log("rename_seqid.gff");

# --------check agat_sq_repeats_analyzer.pl-------------

# XXX

# --------check agat_sq_reverse_complement.pl-------------

$script = $script_prefix."bin/agat_sq_reverse_complement.pl";
$result = "$output_folder/agat_sq_reverse_complement_1.gff";
system(" $script --gff $input_folder/1.gff --fasta  $input_folder/1.fa -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
cleaning_log("1.gff");
unlink "1_rt.fa";

# --------check agat_sq_rfam_analyzer.pl-------------

# XXX

# --------check agat_sq_split.pl-------------

# XXX

# --------check agat_sq_stat_basic.pl-------------

$script = $script_prefix."bin/agat_sq_stat_basic.pl";
$result = "$output_folder/agat_sq_stat_basic_1.gff";
system(" $script --gff $input_folder/1.gff -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
cleaning_log("1.gff");

# --- convenient function ---

sub cleaning_log{
  my ($filename, $local_config)=@_;

  # REMOVE LOG folder if a file name is provided
  if($filename){
    my ($name, $path, $suffix) = fileparse($filename, qr/\.[^.]*/);
    if (-e "agat_log_$name"){
      remove_tree( "agat_log_$name" );
    }
  }

  # remove config
  # Si la variable $local_config est définie et vraie (c’est-à-dire qu’elle n’est pas undef, ni vide, ni 0), alors $config_to_remove prendra sa valeur. Sinon, il prendra la valeur de $config.
  my $config_to_remove = $local_config ? $local_config : $config;
  if (-e $config){
    unlink $config;
  }

  # the rest
  opendir(my $dh, ".") or die "Cannot open current directory: $!";
  while (my $file = readdir($dh)) {
    next if $file =~ /^\./;  # Skip hidden files and . / ..
    if ($file =~ /^\Q$outprefix\E/ or $file =~ /^\Q$outtmp\E/) {
      if ( -d $file ) {
        remove_tree($file);
        #print "remove_tree $file\n";
      } else {
        unlink $file;
        #print "unlink $file\n";
      }
    }
  }
  closedir($dh);
}