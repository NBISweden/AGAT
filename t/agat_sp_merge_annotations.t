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

# ------------------- check agat_sp_merge_annotations script-------------------

my $script = $script_prefix . catfile($bin_dir, "agat_sp_merge_annotations.pl");
my $result = "$output_folder/agat_sp_merge_annotations_1.gff";
{
    my $dir = setup_tempdir();
    my $outtmp = catfile($dir, 'tmp.gff');
    my $outprefix = catfile($dir, 'tmp');
    system(" $script --gff $input_folder/agat_sp_merge_annotations/file1.gff  --gff $input_folder/agat_sp_merge_annotations/file2.gff -o $outtmp 2>&1 1>/dev/null");
    check_diff( $outtmp, $result, "output $script" );
}


$result = "$output_folder/agat_sp_merge_annotations_2.gff";
{
    my $dir = setup_tempdir();
    my $outtmp = catfile($dir, 'tmp.gff');
    my $outprefix = catfile($dir, 'tmp');
    system(" $script --gff $input_folder/agat_sp_merge_annotations/fileA.gff  --gff $input_folder/agat_sp_merge_annotations/fileB.gff -o $outtmp 2>&1 1>/dev/null");
    check_diff( $outtmp, $result, "output $script" );
}


$result = "$output_folder/agat_sp_merge_annotations_3.gff";
{
    my $dir = setup_tempdir();
    my $outtmp = catfile($dir, 'tmp.gff');
    my $outprefix = catfile($dir, 'tmp');
    system(" $script --gff $input_folder/agat_sp_merge_annotations/test457_A.gff  --gff $input_folder/agat_sp_merge_annotations/test457_B.gff -o $outtmp 2>&1 1>/dev/null");
    check_diff( $outtmp, $result, "output $script" );
}




done_testing();
