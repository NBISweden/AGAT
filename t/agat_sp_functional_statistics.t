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

# --------check agat_sp_functional_statistics.pl-------------

my $script = $script_prefix . catfile($bin_dir, "agat_sp_functional_statistics.pl");
{ my $dir = setup_tempdir(); ok(system("$script -h 1>\/dev\/null") == 0, "help $script"); }

my $result = "$output_folder/agat_sp_functional_statistics_1.txt";
{
    my $dir = setup_tempdir();
    my $outtmp = catfile($dir, 'tmp.gff');
    my $outprefix = catfile($dir, 'tmp');
    system(" $script --gff $input_folder/function.gff -o $outtmp 2>&1 1>/dev/null");
    check_diff( "$outtmp/gene\@mrna/table_per_feature_type.txt", "$output_folder/agat_sp_functional_statistics/table_gene_mrna.txt", "output $script" );
    check_diff( "$outtmp/repeat_region/table_per_feature_type.txt", "$output_folder/agat_sp_functional_statistics/table_repeat.txt", "output $script" );
}




done_testing();
