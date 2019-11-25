#!/usr/bin/env perl

use strict;
use warnings;
use File::Path;
use Test::More tests => 28;

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

# -------------------------- check agat_sp_add_introns --------------------------

$script = $script_prefix."bin/agat_sp_add_introns.pl";
$result = "$output_folder/agat_sp_add_introns_1.gff";
system(" $script --gff $output_folder/1.gff -o $outtmp 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
unlink $outtmp;


# -------------------------- check agat_sp_add_start_and_stop --------------------------
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

# XXX

# --------check agat_sp_extract_sequences.pl-------------

$script = $script_prefix."bin/agat_sp_extract_sequences.pl";
$result = "$output_folder/agat_sp_extract_sequences_1.fa";
system(" $script --gff $output_folder/1.gff --fasta $output_folder/1.fa -o $outtmp 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
unlink $outtmp;

# --------check agat_sp_extract_attributes.pl.pl-------------

# XXX

# --------check agat_sp_filter_by_ORF_size.pl-------------
# /!\ Two outputs with difficult characters as > <
#$script = $script_prefix."bin/agat_sp_filter_by_ORF_size.pl";
#$result = "$output_folder/agat_sp_filter_by_ORF_size\>100.gff";
#$result2 = "$output_folder/agat_sp_filter_by_ORF_size_NOT>100.gff";
#system(" $script --gff $output_folder/1.gff -o $outtmp ");
#run test
#ok( system("diff $result $outprefix\>100.gff") == 0, "output $script");
#unlink $outtmp;

# --------check agat_sp_filter_by_locus_distance.pl-------------
$script = $script_prefix."bin/agat_sp_filter_by_locus_distance.pl";
$result = "$output_folder/agat_sp_filter_by_locus_distance_1.gff";
system(" $script --gff $output_folder/1.gff -o $outtmp 1>/dev/null");
ok( system("diff $result $outtmp") == 0, "output $script");
unlink $outtmp;

# --------check agat_sp_filter_by_mrnaBlastValue.pl-------------

# XXX

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

# --------check agat_sp_fix_cds_frame.pl-------------

$script = $script_prefix."bin/agat_sp_fix_cds_frame.pl";
$result = "$output_folder/agat_sp_fix_cds_frame_1.gff";
system(" $script --gff $output_folder/1.gff --fasta $output_folder/1.fa -o $outtmp 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
unlink $outtmp;

# --------check agat_sp_fix_features_locations_duplicated.pl-------------

# XXX

# --------check agat_sp_fix_fusion.pl-------------

$script = $script_prefix."bin/agat_sp_fix_fusion.pl";
$result = "$output_folder/agat_sp_fix_fusion_1.txt"; # txt becaus
system(" $script --gff $output_folder/1.gff --fasta $output_folder/1.fa -o $outtmp 1>/dev/null");
#run test
ok( system("diff -b -I '^Job done in' -I '^Job done in' $result $outprefix-report.txt") == 0, "output $script");
unlink "$outprefix-report.txt";
unlink "$outprefix-all.gff";
unlink "$outprefix-intact.gff";
unlink "$outprefix-only_modified.gff";

# --------agat_sp_fix_longest_ORF.pl-------------

$script = $script_prefix."bin/agat_sp_fix_longest_ORF.pl";
$result = "$output_folder/agat_sp_fix_longest_ORF_1.txt"; # txt becaus
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
ok( system("diff $result $outprefix/report.txt") == 0, "output $script");
rmtree $outprefix;


# --------check agat_sp_gxf_to_gff3.pl-------------

# XXX NOT needed done by gff_syntax test

# --------check agat_sp_keep_longest_isoform.pl-------------

$script = $script_prefix."bin/agat_sp_keep_longest_isoform.pl";
$result = "$output_folder/agat_sp_keep_longest_isoform_1.gff";
system(" $script --gff $output_folder/1.gff -o $outtmp 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
unlink $outtmp;

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
system(" $script --gff $output_folder/1.gff --prefix NBIS -o $outtmp 1>/dev/null");
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

# XXX


# --------check agat_sp_manage_functional_annotation.pl-------------

# XXX


# --------check agat_sp_manage_introns.pl-------------

# XXX

# ------------------- check agat_sp_merge_annotations script-------------------
$script = $script_prefix."bin/agat_sp_merge_annotations.pl";
$result = "$output_folder/agat_sp_merge_annotations_1.gff";
system(" $script --gff t/gff_syntax/25_test.gff  --gff t/gff_syntax/9_test.gff -o $outtmp 1>/dev/null");
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
$result = "$output_folder/agat_sp_statistics_1.gff";
system(" $script --gff $output_folder/1.gff -o $outtmp 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
unlink $outtmp;

# --------check agat_sp_to_tabulated.pl-------------

$script = $script_prefix."bin/agat_sp_to_tabulated.pl";
$result = "$output_folder/agat_sp_to_tabulated_1.txt";
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
ok( system("diff -b -I '^Job done in' -I '^Job done in' $result $outtmp") == 0, "output $script");
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
