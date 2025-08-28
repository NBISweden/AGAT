#!/usr/bin/env perl
use strict;
use warnings;
use File::Path;
use FindBin qw($Bin);
use lib "$Bin/lib";
use File::Spec::Functions qw(catfile catdir);
use Cwd qw(abs_path);
use AGAT::TestUtilities qw(setup_tempdir script_prefix check_quiet_and_normal_run); 
use Test::More;

my $script_prefix = script_prefix();
my $root = abs_path(catdir($Bin, '..'));
my $bin_dir = catdir($root, 'bin');
my $input_folder = catdir($Bin, 'scripts_output', 'in');
my $output_folder = catdir($Bin, 'scripts_output', 'out');
my $config = 'agat_config.yaml';

# ------------------- check agat_sp_compare_two_annotations script-------------------
my $script = $script_prefix . catfile($bin_dir, "agat_sp_compare_two_annotations.pl");
{ my $dir = setup_tempdir(); ok(system("$script -h 1>\/dev\/null") == 0, "help $script"); }

my $result = "$output_folder/agat_sp_compare_two_annotations_1.txt";
check_quiet_and_normal_run(
    $script,
    { gff1 => "$input_folder/1.gff", gff2 => "$input_folder/1.gff" },
    "$result.stdout",
    $result,
    '.gff/report.txt'
);


$script = $script_prefix . catfile($bin_dir, "agat_sp_compare_two_annotations.pl");
$result = "$output_folder/agat_sp_compare_two_annotations_2.txt";
check_quiet_and_normal_run(
    $script,
    { gff1 => "$input_folder/agat_sp_compare_two_annotations/file1.gff", gff2 => "$input_folder/agat_sp_compare_two_annotations/file2.gff" },
    "$result.stdout",
    $result,
    '.gff/report.txt'
);


$script = $script_prefix . catfile($bin_dir, "agat_sp_compare_two_annotations.pl");
$result = "$output_folder/agat_sp_compare_two_annotations_3.txt";
check_quiet_and_normal_run(
    $script,
    { gff1 => "$input_folder/agat_sp_compare_two_annotations/file2.gff", gff2 => "$input_folder/agat_sp_compare_two_annotations/file1.gff" },
    "$result.stdout",
    $result,
    '.gff/report.txt'
);




done_testing();
