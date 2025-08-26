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

# --------check agat_sp_functional_statistics.pl-------------

my $script = $script_prefix . catfile($bin_dir, "agat_sp_functional_statistics.pl");
{ my $dir = setup_tempdir(); ok(system("$script -h 1>\/dev\/null") == 0, "help $script"); }

my $result  = "$output_folder/agat_sp_functional_statistics_1.txt";
my $result1 = "$output_folder/agat_sp_functional_statistics/table_gene_mrna.txt";
my $result2 = "$output_folder/agat_sp_functional_statistics/table_repeat.txt";
check_quiet_and_normal_run(
    $script,
    { gff => "$input_folder/function.gff" },
    "$result.stdout",
    [ $result1, $result2 ],
    [
        '.gff/gene@mrna/table_per_feature_type.txt',
        '.gff/repeat_region/table_per_feature_type.txt',
    ]
);
done_testing();
