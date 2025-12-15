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

# ------------------- check agat_sp_compare_two_annotations script-------------------

my $script = $script_prefix . catfile($bin_dir, "agat_sp_compare_two_annotations.pl");

{ ok(system("$script -h 1>\/dev\/null") == 0, "help $script"); }


{
    my $dir = setup_tempdir();
    my $out = catfile($dir,'tmp');
    my $result = "$output_folder/agat_sp_compare_two_annotations_1.txt";
    check_quiet_run(
        "$script --gff1 " . catfile($input_folder,'1.gff') .
        " --gff2 " . catfile($input_folder,'1.gff') .
        " -o $out"
    );
    my $outtmp = catfile($out,'report.txt');
    check_diff($outtmp, $result, 'output agat_sp_compare_two_annotations.pl case1');
}

{
    my $dir = setup_tempdir();
    my $out = catfile($dir,'tmp');
    my $result = "$output_folder/agat_sp_compare_two_annotations_2.txt";
    check_quiet_run(
        "$script --gff1 " . catfile($input_folder,'agat_sp_compare_two_annotations','file1.gff') .
        " --gff2 " . catfile($input_folder,'agat_sp_compare_two_annotations','file2.gff') .
        " -o $out"
    );
    my $outtmp = catfile($out,'report.txt');
    check_diff($outtmp, $result, 'output agat_sp_compare_two_annotations.pl case2');
}

{
    my $dir = setup_tempdir();
    my $out = catfile($dir,'tmp');
    my $result = "$output_folder/agat_sp_compare_two_annotations_3.txt";
    check_quiet_run(
        "$script --gff1 " . catfile($input_folder,'agat_sp_compare_two_annotations','file2.gff') .
        " --gff2 " . catfile($input_folder,'agat_sp_compare_two_annotations','file1.gff') .
        " -o $out"
    );
    my $outtmp = catfile($out,'report.txt');
    check_diff($outtmp, $result, 'output agat_sp_compare_two_annotations.pl case3');
}

done_testing();
