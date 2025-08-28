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

# -------------------------- check agat_convert_mfannot2gff -------------------------

my $script = $script_prefix . catfile($bin_dir, "agat_convert_mfannot2gff.pl");
{ my $dir = setup_tempdir(); ok(system("$script -h 1>\/dev\/null") == 0, "help $script"); }

my $result = "$output_folder/agat_convert_mfannot2gff_1.gff";
check_quiet_and_normal_run(
    $script,
    { mfannot => "$input_folder/test.mfannot" },
    "$result.stdout",
    $result
);


$script = $script_prefix . catfile($bin_dir, "agat_convert_mfannot2gff.pl");
$result = "$output_folder/agat_convert_mfannot2gff_2.gff";
check_quiet_and_normal_run(
    $script,
    { mfannot => "$input_folder/test.mfannot2" },
    "$result.stdout",
    $result
);




done_testing();
