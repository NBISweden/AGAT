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

# --------check agat_sp_functional_statistics.pl-------------

my $script = $script_prefix . catfile($bin_dir, "agat_sp_functional_statistics.pl");

{ ok(system("$script -h 1>\/dev\/null") == 0, "help $script"); }

{
    my $dir = setup_tempdir();
    my $outtmp = catfile($dir,'tmp');
    my $result1 = "$output_folder/agat_sp_functional_statistics/table_gene_mrna.txt";
    my $result2 = "$output_folder/agat_sp_functional_statistics/table_repeat.txt";
    check_quiet_run(
        "$script --gff " . catfile($input_folder,'function.gff') . " -o $outtmp"
    );
    my $out1 = catfile($outtmp,'gene@mrna','table_per_feature_type.txt');
    my $out2 = catfile($outtmp,'repeat_region','table_per_feature_type.txt');
    check_diff($out1, $result1, 'functional stats gene@mrna');
    check_diff($out2, $result2, 'functional stats repeat_region');
}
done_testing();
