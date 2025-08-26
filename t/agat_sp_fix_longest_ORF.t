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

# --------agat_sp_fix_longest_ORF.pl-------------

my $script = $script_prefix . catfile( $bin_dir, "agat_sp_fix_longest_ORF.pl" );
{ my $dir = setup_tempdir(); ok(system("$script -h 1>\/dev\/null") == 0, "help $script"); }

my $result = "$output_folder/agat_sp_fix_longest_ORF_1.txt";
check_quiet_and_normal_run(
    $script,
    { gff => "$input_folder/1.gff", fasta => "$input_folder/1.fa" },
    $result,
    '-report.txt'
);


done_testing();
