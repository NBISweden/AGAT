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

# ---------------------- check agat_sp_add_start_and_stop ----------------------

my $script = $script_prefix . catfile($bin_dir, "agat_sp_add_start_and_stop.pl");

{ ok(system("$script -h 1>\/dev\/null") == 0, "help $script"); }


{
    my $dir = setup_tempdir();
    my $result = "$output_folder/agat_sp_add_start_and_stop_1.gff";
    my $outtmp = catfile($dir,'tmp.gff');
    check_quiet_run(
        "$script --gff " . catfile($input_folder,'agat_sp_add_start_and_stop.gff') .
        " --fasta " . catfile($input_folder,'1.fa') .
        " --ni -o $outtmp"
    );
    check_diff($outtmp, $result, 'output agat_sp_add_start_and_stop.pl case1');
}


{
    my $dir = setup_tempdir();
    my $result = "$output_folder/agat_sp_add_start_and_stop_2.gff";
    my $outtmp = catfile($dir,'tmp.gff');
    check_quiet_run(
        "$script --gff " . catfile($input_folder,'agat_sp_add_start_and_stop.gff') .
        " --fasta " . catfile($input_folder,'1.fa') .
        " -e --ni -o $outtmp"
    );
    check_diff($outtmp, $result, 'output agat_sp_add_start_and_stop.pl case2');
}

done_testing();
