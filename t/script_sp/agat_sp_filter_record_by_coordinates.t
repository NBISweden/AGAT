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

# --------check agat_sp_filter_record_by_coordinates.pl-------------

my $script = $script_prefix . catfile($bin_dir, "agat_sp_filter_record_by_coordinates.pl");

{ ok(system("$script -h 1>\/dev\/null") == 0, "help $script"); }

{
    my $dir = setup_tempdir();
    my $outtmp = catfile($dir,'tmp.gff');
    my $result = "$output_folder/agat_sp_filter_record_by_coordinates_1.gff";
    check_quiet_run(
        "$script --gff " . catfile($input_folder,'1.gff') .
        " --tsv " . catfile($input_folder,'coordinates.tsv') .
        " -o $outtmp"
    );
    check_diff("$outtmp/remaining.gff3", $result, 'output filter_record_by_coordinates');
}

done_testing();