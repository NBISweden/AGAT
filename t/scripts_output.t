#!/usr/bin/env perl

use strict;
use warnings;
use File::Path;
use Test::More tests => 43;

=head1 DESCRIPTION

Test to verify the output of scripts using gff_syntax cases

=cut

# Check if has to be run in Devel::Cover or not
my $script_prefix="";
if (exists $ENV{'HARNESS_PERL_SWITCHES'} ) {
  if ($ENV{'HARNESS_PERL_SWITCHES'} =~ m/Devel::Cover/) {
    $script_prefix="perl -MDevel::Cover ";
  }
}

# shared variables
my $output_folder = "t/scripts_output";
my $outtmp = "tmp.gff"; # path file where to save temporary output
my $outprefix = "tmp";
my $script;
my $result;
my $result2;

# -------------------------- check agat_convert_bed2gff -------------------------

$script = $script_prefix."bin/agat_convert_bed2gff.pl";
$result = "$output_folder/agat_convert_bed2gff_1.gff";
system(" $script --bed $output_folder/test.bed -o $outtmp 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
unlink $outtmp;

# -------------------------- check agat_convert_embl2gff -------------------------

# Nothing yet

# -------------------------- check agat_convert_genscan2gff -------------------------

$script = $script_prefix."bin/agat_convert_genscan2gff.pl";
$result = "$output_folder/agat_convert_genscan2gff_1.gff";
system(" $script --genscan $output_folder/test.genscan -o $outtmp 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
unlink $outtmp;

# -------------------------- check agat_convert_mfannot2gff -------------------------

$script = $script_prefix."bin/agat_convert_mfannot2gff.pl";
$result = "$output_folder/agat_convert_mfannot2gff_1.gff";
system(" $script --mfannot $output_folder/test.mfannot -o $outtmp 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
unlink $outtmp;

# -------------------------- check agat_convert_sp_gff2bed -------------------------

$script = $script_prefix."bin/agat_convert_sp_gff2bed.pl";
$result = "$output_folder/agat_convert_sp_gff2bed_1.gff";
system(" $script --gff $output_folder/1.gff -o $outtmp 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
unlink $outtmp;

# -------------------------- check agat_convert_sp_gff2gtf -------------------------

$script = $script_prefix."bin/agat_convert_sp_gff2gtf.pl";
$result = "$output_folder/agat_convert_sp_gff2gtf_1.gff";
system(" $script --gff $output_folder/1.gff -o $outtmp 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
unlink $outtmp;

# --------check agat_sp_to_tabulated.pl-------------

$script = $script_prefix."bin/agat_convert_sp_gff2tsv.pl";
$result = "$output_folder/agat_convert_sp_gff2tsv_1.tsv";
system(" $script --gff $output_folder/1.gff -o $outtmp 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
unlink $outtmp;

# -------------------------- check agat_convert_sp_gff2zff -------------------------

$script = $script_prefix."bin/agat_convert_sp_gff2zff.pl";
$result = "$output_folder/agat_convert_sp_gff2zff_1.gff";
system(" $script --gff $output_folder/1.gff --fasta $output_folder/1.fa -o $outtmp 1>/dev/null");
#run test
ok( system("diff $result $outprefix.ann") == 0, "output $script");
unlink $outprefix.".ann";
unlink $outprefix.".dna";

# -------------------------- check agat_convert_sp_gxf2gxf.pl -------------------------

# No need to be tested, it is the same as agat_sp_gxf_to_gff3

# -------------------------- check agat_sp_add_introns -------------------------

$script = $script_prefix."bin/agat_sp_add_introns.pl";
$result = "$output_folder/agat_sp_add_introns_1.gff";
system(" $script --gff $output_folder/1.gff -o $outtmp 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
unlink $outtmp;

# ---------------------- check agat_sp_add_start_and_stop ----------------------

$script = $script_prefix."bin/agat_sp_add_start_and_stop.pl";
$result = "$output_folder/agat_sp_add_start_and_stop_1.gff";
system(" $script --gff $output_folder/1.gff --fasta $output_folder/1.fa -o $outtmp 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
unlink $outtmp;

# -------------------------- check agat_sp_statistics --------------------------

$script = $script_prefix."bin/agat_sp_alignment_output_style.pl";
$result = "$output_folder/agat_sp_alignment_output_style_1.gff";
system(" $script --gff t/gff_syntax/0_test.gff -o $outtmp 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
unlink $outtmp;

# --------check agat_sp_clipN_seqExtremities_and_fixCoordinates.pl-------------

# XXX

# ------------------- check agat_sp_complement_annotations script-------------------
$script = $script_prefix."bin/agat_sp_compare_two_annotations.pl";
$result = "$output_folder/agat_sp_compare_two_annotations_1.txt";
system(" $script --gff1 $output_folder/1.gff  --gff2 $output_folder/1.gff -o $outtmp 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
unlink $outtmp;

# ------------------- check agat_sp_complement_annotations script-------------------
$script = $script_prefix."bin/agat_sp_complement_annotations.pl";
$result = "$output_folder/agat_sp_complement_annotations_1.gff";
system(" $script --ref t/gff_syntax/25_test.gff  --add t/gff_syntax/9_test.gff -o $outtmp 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
unlink $outtmp;

# --------check agat_sp_ensembl_output_style.pl-------------

$script = $script_prefix."bin/agat_sp_ensembl_output_style.pl";
$result = "$output_folder/agat_sp_ensembl_output_style_1.gff";
system("$script --gff t/gff_syntax/0_test.gff -o $outtmp 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
unlink $outtmp;

# --------check agat_sp_extract_attributes.pl-------------

$script = $script_prefix."bin/agat_sp_extract_attributes.pl";
$result = "$output_folder/agat_sp_extract_attributes_1.txt";
system(" $script --gff $output_folder/1.gff --att protein_id -o $outtmp 1>/dev/null");
#run test
ok( system("diff $result ${outprefix}_protein_id.gff") == 0, "output $script");
unlink $outprefix."_protein_id.gff";

# --------check agat_sp_extract_sequences.pl-------------

$script = $script_prefix."bin/agat_sp_extract_sequences.pl";
$result = "$output_folder/agat_sp_extract_sequences_1.fa";
system(" $script --gff $output_folder/1.gff --fasta $output_folder/1.fa -o $outtmp 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
unlink $outtmp;

# --------check agat_sp_filter_by_locus_distance.pl-------------
$script = $script_prefix."bin/agat_sp_filter_by_locus_distance.pl";
$result = "$output_folder/agat_sp_filter_by_locus_distance_1.gff";
system(" $script --gff $output_folder/1.gff -o $outtmp 1>/dev/null");
ok( system("diff $result $outtmp") == 0, "output $script");
unlink $outtmp;

# --------check agat_sp_filter_by_mrnaBlastValue.pl-------------

# XXX

# --------check agat_sp_filter_by_ORF_size.pl-------------
$script = $script_prefix."bin/agat_sp_filter_by_ORF_size.pl";
$result = "$output_folder/agat_sp_filter_by_ORF_size_sup100.gff";
system(" $script --gff $output_folder/1.gff -o $outtmp 1>/dev/null");
#run test
ok( system("diff $result ${outprefix}_sup100.gff") == 0, "output $script");
unlink $outprefix."_sup100.gff";
unlink  $outprefix."_NOT_sup100.gff";

# --------check agat_sp_filter_feature_by_attribute_presence.pl-------------
$script = $script_prefix."bin/agat_sp_filter_feature_by_attribute_presence.pl";
$result = "$output_folder/agat_sp_filter_feature_by_attribute_presence_1.gff";
system(" $script --gff $output_folder/1.gff -o $outtmp -a protein_id 1>/dev/null");
ok( system("diff $result $outtmp") == 0, "output $script");
unlink $outtmp;
unlink $outprefix."_discarded.txt";
unlink $outprefix."_report.txt";


# --------check agat_sp_filter_feature_by_attribute_value.pl-------------
$script = $script_prefix."bin/agat_sp_filter_feature_by_attribute_value.pl";
$result = "$output_folder/agat_sp_filter_feature_by_attribute_value_1.gff";
system(" $script --gff $output_folder/1.gff -o $outtmp --value Os01t0100100-01 -p level3 -a protein_id 1>/dev/null");
ok( system("diff $result $outtmp") == 0, "output $script");
unlink $outtmp;
unlink $outprefix."_discarded.txt";
unlink $outprefix."_report.txt";

# --------check agat_sp_filter_incomplete_gene_coding_models.pl-------------

$script = $script_prefix."bin/agat_sp_filter_incomplete_gene_coding_models.pl";
$result = "$output_folder/agat_sp_filter_incomplete_gene_coding_models_1.gff";
$result2 = "$output_folder/agat_sp_filter_incomplete_gene_coding_models_incomplete_1.gff";
system(" $script --gff $output_folder/1.gff --fasta $output_folder/1.fa -o $outtmp 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
ok( system("diff $result2 $outprefix"."_incomplete.gff") == 0, "output $script");
unlink $outtmp;
unlink $outprefix."_incomplete.gff";

# --------check agat_sp_fix_cds_phases.pl-------------

$script = $script_prefix."bin/agat_sp_fix_cds_phases.pl";
$result = "$output_folder/agat_sp_fix_cds_phases_1.gff";
system(" $script --gff $output_folder/1.gff --fasta $output_folder/1.fa -o $outtmp 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
unlink $outtmp;

# --------check agat_sp_fix_features_locations_duplicated.pl-------------

# XXX

# --------check agat_sp_fix_fusion.pl-------------

$script = $script_prefix."bin/agat_sp_fix_fusion.pl";
$result = "$output_folder/agat_sp_fix_fusion_1.txt";
system(" $script --gff $output_folder/1.gff --fasta $output_folder/1.fa -o $outtmp 1>/dev/null");
#run test
ok( system("diff -b -I '^Job done in' -I '^Job done in' $result $outprefix-report.txt") == 0, "output $script");
unlink "$outprefix-report.txt";
unlink "$outprefix-all.gff";
unlink "$outprefix-intact.gff";
unlink "$outprefix-only_modified.gff";

# --------agat_sp_fix_longest_ORF.pl-------------

$script = $script_prefix."bin/agat_sp_fix_longest_ORF.pl";
$result = "$output_folder/agat_sp_fix_longest_ORF_1.txt";
system(" $script --gff $output_folder/1.gff --fasta $output_folder/1.fa -o $outtmp 1>/dev/null");
#run test
ok( system("diff -b -I '^Job done in' -I '^Job done in' $result $outprefix-report.txt") == 0, "output $script");
unlink "$outprefix-report.txt";
unlink "$outprefix-all.gff";
unlink "$outprefix-intact.gff";
unlink "$outprefix-only_modified.gff";

# --------check agat_sp_fix_overlaping_genes.pl-------------

# XXX

# --------check agat_sp_fix_small_exon_from_extremities.pl-------------

# XXX

# --------check agat_sp_flag_short_introns.pl-------------

# XXX

# --------check agat_sp_functional_statistics.pl-------------

$script = $script_prefix."bin/agat_sp_functional_statistics.pl";
$result = "$output_folder/agat_sp_functional_statistics_1.txt";
system(" $script --gff t/gff_syntax/10_test.gff -o $outtmp 1>/dev/null");
#run test
ok( system("diff $result $outtmp/report.txt") == 0, "output $script");
rmtree $outtmp;


# --------check agat_sp_gxf_to_gff3.pl-------------

# XXX NOT needed done by gff_syntax test

# --------check agat_sp_keep_longest_isoform.pl-------------

$script = $script_prefix."bin/agat_sp_keep_longest_isoform.pl";
$result = "$output_folder/agat_sp_keep_longest_isoform_1.gff";
system(" $script --gff $output_folder/1.gff -o $outtmp 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
unlink $outtmp;

# --------check agat_sp_kraken_assess_liftover.pl-------------

$script = $script_prefix."bin/agat_sp_kraken_assess_liftover.pl";
$result = "$output_folder/agat_sp_kraken_assess_liftover_1.gff";
system(" $script --gtf $output_folder/test_kraken.gtf -o $outtmp 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
unlink $outprefix;
unlink $outprefix."_report.txt";
unlink $outprefix."-geneMapped_plot.pdf";
unlink $outprefix."-geneMapped.txt";

# --------check agat_sp_list_short_introns.pl-------------

$script = $script_prefix."bin/agat_sp_list_short_introns.pl";
$result = "$output_folder/agat_sp_list_short_introns_1.txt";
system(" $script --gff $output_folder/1.gff -o $outtmp 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
unlink $outtmp;

# --------check agat_sp_load_function_from_protein_align.pl-------------

# XXX


# --------check agat_sp_manage_IDs.pl-------------

$script = $script_prefix."bin/agat_sp_manage_IDs.pl";
$result = "$output_folder/agat_sp_manage_IDs_1.gff";
system(" $script --gff $output_folder/1.gff --prefix NBIS --ensembl --tair --type_dependent -o $outtmp 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
unlink $outtmp;

# --------check agat_sp_manage_UTRs.pl-------------

$script = $script_prefix."bin/agat_sp_manage_UTRs.pl";
$result = "$output_folder/agat_sp_manage_UTRs_1.gff";
system(" $script --gff $output_folder/1.gff -b -o $outtmp 1>/dev/null");
#run test
ok( system("diff $result $outprefix/1_bothSides_under5.gff") == 0, "output $script");
rmtree $outprefix;

# --------check agat_sp_manage_attributes.pl-------------

$script = $script_prefix."bin/agat_sp_manage_attributes.pl";
$result = "$output_folder/agat_sp_manage_attributes_1.gff";
system(" $script --gff $output_folder/1.gff --att protein_id -o $outtmp 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
unlink $outtmp;


# --------check agat_sp_manage_functional_annotation.pl-------------

# XXX


# --------check agat_sp_manage_introns.pl-------------

$script = $script_prefix."bin/agat_sp_manage_introns.pl";
$result = "$output_folder/agat_sp_manage_introns_1.txt";
system(" $script --gff $output_folder/1.gff -o $outtmp 1>/dev/null");
#run test
ok( system( "diff $result $outtmp/report.txt" ) == 0, "output $script");
rmtree $outtmp;

# ------------------- check agat_sp_merge_annotations script-------------------
$script = $script_prefix."bin/agat_sp_merge_annotations.pl";
$result = "$output_folder/agat_sp_merge_annotations_1.gff";
system(" $script --gff t/gff_syntax/25_test.gff  --gff t/gff_syntax/9_test.gff -o $outtmp 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
unlink $outtmp;

# ------------------- check agat_sp_sensitivity_specificity script-------------------
$script = $script_prefix."bin/agat_sp_sensitivity_specificity.pl";
$result = "$output_folder/agat_sp_sensitivity_specificity_1.txt";
system(" $script --gff1 $output_folder/1.gff --gff2 $output_folder/1.gff -o $outtmp 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
unlink $outtmp;

# --------check agat_sp_split_by_level2_feature.pl-------------

$script = $script_prefix."bin/agat_sp_split_by_level2_feature.pl";
$result = "$output_folder/agat_sp_split_by_level2_feature_1.gff";
system(" $script --gff $output_folder/1.gff -o $outtmp 1>/dev/null");
#run test
ok( system("diff $result $outprefix/trna.gff") == 0, "output $script");
rmtree $outprefix;

# -------------------------- check agat_sp_statistics --------------------------

$script = $script_prefix."bin/agat_sp_statistics.pl";
$result = "$output_folder/agat_sp_statistics_1.txt";
system(" $script --gff $output_folder/1.gff -o $outtmp 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
unlink $outtmp;

# --------check agat_sp_webApollo_compliant.pl-------------

$script = $script_prefix."bin/agat_sp_webApollo_compliant.pl";
$result = "$output_folder/agat_sp_webApollo_compliant_1.gff";
system(" $script --gff $output_folder/1.gff -o $outtmp 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
unlink $outtmp;


# --------check agat_sq_add_hash_tag.pl-------------

$script = $script_prefix."bin/agat_sq_add_hash_tag.pl";
$result = "$output_folder/agat_sq_add_hash_tag_1.gff";
system(" $script --gff $output_folder/1.gff -i 2 -o $outtmp 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
unlink $outtmp;

# --------check agat_sq_add_locus_tag.pl-------------

$script = $script_prefix."bin/agat_sq_add_locus_tag.pl";
$result = "$output_folder/agat_sq_add_locus_tag_1.gff";
system(" $script --gff $output_folder/1.gff -o $outtmp 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
unlink $outtmp;

# --------check agat_sq_keep_annotation_from_fastaSeq.pl-------------

# XXX


# --------check agat_sq_list_attributes.pl-------------

$script = $script_prefix."bin/agat_sq_list_attributes.pl";
$result = "$output_folder/agat_sq_list_attributes_1.txt";
system(" $script --gff $output_folder/1.gff -o $outtmp 1>/dev/null");
#run test
ok( system("diff -b -I '^Job done in' $result $outtmp") == 0, "output $script");
unlink $outtmp;

# --------check agat_sq_manage_ID.pl-------------

$script = $script_prefix."bin/agat_sq_manage_ID.pl";
$result = "$output_folder/agat_sq_manage_ID_1.gff";
system(" $script --gff $output_folder/1.gff -o $outtmp 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
unlink $outtmp;

# --------check agat_sq_mask.pl-------------

# XXX


# --------check agat_sq_remove_redundant_entries.pl-------------

$script = $script_prefix."bin/agat_sq_remove_redundant_entries.pl";
$result = "$output_folder/agat_sq_remove_redundant_entries_1.gff";
system(" $script --gff $output_folder/1.gff -o $outtmp 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
unlink $outtmp;

# --------check agat_sq_repeats_analyzer.pl-------------

# XXX

# --------check agat_sq_rfam_analyzer.pl-------------

# XXX

# --------check agat_sq_split.pl-------------

# XXX

# --------check agat_sq_stat_basic.pl-------------

$script = $script_prefix."bin/agat_sq_stat_basic.pl";
$result = "$output_folder/agat_sq_stat_basic_1.gff";
system(" $script --gff $output_folder/1.gff -o $outtmp 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
unlink $outtmp;
