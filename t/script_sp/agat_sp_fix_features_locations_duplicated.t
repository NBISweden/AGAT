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

# --------check agat_sp_fix_features_locations_duplicated.pl-------------
# removed because order can change. So not reproducible at 100%

my $script = $script_prefix . catfile($bin_dir, "agat_sp_fix_features_locations_duplicated.pl");

{ ok(system("$script -h 1>\/dev\/null") == 0, "help $script"); }

{
    my $dir = setup_tempdir();
    my $outtmp = catfile($dir,'tmp.gff');
    my $result = "$output_folder/agat_sp_fix_features_locations_duplicated_1.gff";
    check_quiet_run(
        "$script --gff " . catfile($input_folder,'agat_sp_fix_features_locations_duplicated/test.gff') . " -o $outtmp"
    );
    check_diff($outtmp, $result, 'output fix_features_locations_duplicated');
}

done_testing();