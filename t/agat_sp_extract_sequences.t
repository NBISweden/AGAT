#!/usr/bin/env perl
use strict;
use warnings;
use File::Path;
use FindBin qw($Bin);
use lib "$Bin/lib";
use File::Spec::Functions qw(catfile catdir);
use Cwd qw(abs_path);
use AGAT::TestUtilities qw(setup_tempdir check_diff script_prefix);
use Test::More;

my $script_prefix = script_prefix();
my $root = abs_path(catdir($Bin, '..'));
my $bin_dir = catdir($root, 'bin');
my $input_folder = catdir($Bin, 'scripts_output', 'in');
my $output_folder = catdir($Bin, 'scripts_output', 'out');
my $config = 'agat_config.yaml';

# --------check agat_sp_extract_sequences.pl-------------

my $script = $script_prefix . catfile($bin_dir, "agat_sp_extract_sequences.pl");
my $result = "$input_folder/agat_sp_extract_sequences/agat_sp_extract_sequences_1.fa";
{
    my $dir = setup_tempdir();
    my $outtmp = catfile($dir, 'tmp.gff');
    my $outprefix = catfile($dir, 'tmp');
    system(" $script --gff $input_folder/1.gff --fasta $input_folder/1.fa -o $outtmp 2>&1 1>/dev/null");
    #run test
    check_diff( $outtmp, $result, "output $script test1" );
}


$script = $script_prefix . catfile($bin_dir, "agat_sp_extract_sequences.pl");
$result = "$input_folder/agat_sp_extract_sequences/agat_sp_extract_sequences_split.fa";
{
    my $dir = setup_tempdir();
    my $outtmp = catfile($dir, 'tmp.gff');
    my $outprefix = catfile($dir, 'tmp');
    system(" $script --gff $input_folder/1.gff --fasta $input_folder/1.fa --split -o $outtmp 2>&1 1>/dev/null");
    #run test
    check_diff( $outtmp, $result, "output $script test2" );
}


$script = $script_prefix . catfile($bin_dir, "agat_sp_extract_sequences.pl");
$result = "$input_folder/agat_sp_extract_sequences/agat_sp_extract_sequences_merge.fa";
{
    my $dir = setup_tempdir();
    my $outtmp = catfile($dir, 'tmp.gff');
    my $outprefix = catfile($dir, 'tmp');
    system(" $script --gff $input_folder/1.gff --fasta $input_folder/1.fa -t exon --merge -o $outtmp 2>&1 1>/dev/null");
    #run test
    check_diff( $outtmp, $result, "output $script test3" );
}


$script = $script_prefix . catfile($bin_dir, "agat_sp_extract_sequences.pl");
$result = "$input_folder/agat_sp_extract_sequences/agat_sp_extract_sequences_full.fa";
{
    my $dir = setup_tempdir();
    my $outtmp = catfile($dir, 'tmp.gff');
    my $outprefix = catfile($dir, 'tmp');
    system(" $script --gff $input_folder/1.gff --fasta $input_folder/1.fa --full -o $outtmp 2>&1 1>/dev/null");
    #run test
    check_diff( $outtmp, $result, "output $script test4" );
}


$script = $script_prefix . catfile($bin_dir, "agat_sp_extract_sequences.pl");
$result = "$input_folder/agat_sp_extract_sequences/agat_sp_extract_sequences_attributes_kept.fa";
{
    my $dir = setup_tempdir();
    my $outtmp = catfile($dir, 'tmp.gff');
    my $outprefix = catfile($dir, 'tmp');
    system(" $script --gff $input_folder/1.gff --fasta $input_folder/1.fa --keep_attributes -o $outtmp 2>&1 1>/dev/null");
    #run test
    check_diff( $outtmp, $result, "output $script test5" );
}


$script = $script_prefix . catfile($bin_dir, "agat_sp_extract_sequences.pl");
$result = "$input_folder/agat_sp_extract_sequences/agat_sp_extract_sequences_parent_attributes_kept.fa";
{
    my $dir = setup_tempdir();
    my $outtmp = catfile($dir, 'tmp.gff');
    my $outprefix = catfile($dir, 'tmp');
    system(" $script --gff $input_folder/1.gff --fasta $input_folder/1.fa --keep_parent_attributes -o $outtmp 2>&1 1>/dev/null");
    #run test
    check_diff( $outtmp, $result, "output $script test6" );
}




done_testing();
