#!/usr/bin/env perl
use strict;
use warnings;
use File::Path;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use File::Spec::Functions qw(catfile catdir);
use Cwd qw(abs_path);
use AGAT::TestUtilities; 
use Test::More;

my $script_prefix = script_prefix();
my $root = abs_path(catdir($Bin, '..', '..'));
my $bin_dir = catdir($root, 'bin');
my $input_folder = catdir($Bin, 'in');
my $output_folder = catdir($Bin, 'out');
my $config = 'agat_config.yaml';

# ------------------- check agat_sp_merge_annotations script-------------------

my $script = $script_prefix . catfile($bin_dir, "agat_sp_merge_annotations.pl");

{ ok(system("$script -h 1>\/dev\/null") == 0, "help $script"); }

{
    my $dir = setup_tempdir();
    my $outtmp = catfile($dir, 'tmp1.gff');
    my $result = "$output_folder/agat_sp_merge_annotations_1.gff";
    check_quiet_run(
        "$script --gff " . catfile($input_folder, 'agat_sp_merge_annotations/file1.gff') .
        " --gff " . catfile($input_folder, 'agat_sp_merge_annotations/file2.gff') .
        " -o $outtmp"
    );
    check_diff($outtmp, $result, 'output agat_sp_merge_annotations.pl case 1');
}

{
    my $dir = setup_tempdir();
    my $outtmp = catfile($dir, 'tmp2.gff');
    my $result = "$output_folder/agat_sp_merge_annotations_2.gff";
    check_quiet_run(
        "$script --gff " . catfile($input_folder, 'agat_sp_merge_annotations/fileA.gff') .
        " --gff " . catfile($input_folder, 'agat_sp_merge_annotations/fileB.gff') .
        " -o $outtmp"
    );
    check_diff($outtmp, $result, 'output agat_sp_merge_annotations.pl case 2');
}

{
    my $dir = setup_tempdir();
    my $outtmp = catfile($dir, 'tmp3.gff');
    my $result = "$output_folder/agat_sp_merge_annotations_3.gff";
    check_quiet_run(
        "$script --gff " . catfile($input_folder, 'agat_sp_merge_annotations/test457_A.gff') .
        " --gff " . catfile($input_folder, 'agat_sp_merge_annotations/test457_B.gff') .
        " -o $outtmp"
    );
    check_diff($outtmp, $result, 'output agat_sp_merge_annotations.pl case 3');
}

done_testing();