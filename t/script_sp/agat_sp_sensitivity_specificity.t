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

# ------------------- check agat_sp_sensitivity_specificity script-------------------

my $script = $script_prefix . catfile($bin_dir, "agat_sp_sensitivity_specificity.pl");

{ ok(system("$script -h 1>\/dev\/null") == 0, "help $script"); }

{
    my $dir = setup_tempdir();
    my $outtmp = catfile($dir,'tmp1.txt');
    my $result = "$output_folder/agat_sp_sensitivity_specificity_1.txt";
    check_quiet_run(
        "$script --gff1 " . catfile($input_folder,'1.gff') .
        " --gff2 " . catfile($input_folder,'1.gff') .
        " -o $outtmp"
    );
    check_diff($outtmp, $result, 'output sensitivity_specificity case1');
}

{
    my $dir = setup_tempdir();
    my $outtmp = catfile($dir,'tmp2.txt');
    my $result = "$output_folder/agat_sp_sensitivity_specificity_2.txt";
    check_quiet_run(
        "$script --gff1 " . catfile($input_folder,'agat_sp_sensitivity_specificity/ref0.gff3') .
        " --gff2 " . catfile($input_folder,'agat_sp_sensitivity_specificity/query0.gff3') .
        " -o $outtmp"
    );
    check_diff($outtmp, $result, 'output sensitivity_specificity case2');
}

{
    my $dir = setup_tempdir();
    my $outtmp = catfile($dir,'tmp3.txt');
    my $result = "$output_folder/agat_sp_sensitivity_specificity_3.txt";
    check_quiet_run(
        "$script --gff1 " . catfile($input_folder,'agat_sp_sensitivity_specificity/ref1.gff3') .
        " --gff2 " . catfile($input_folder,'agat_sp_sensitivity_specificity/query1.gff3') .
        " -o $outtmp"
    );
    check_diff($outtmp, $result, 'output sensitivity_specificity case3');
}

{
    my $dir = setup_tempdir();
    my $outtmp = catfile($dir,'tmp4.txt');
    my $result = "$output_folder/agat_sp_sensitivity_specificity_4.txt";
    check_quiet_run(
        "$script --gff1 " . catfile($input_folder,'agat_sp_sensitivity_specificity/ref2.gff3') .
        " --gff2 " . catfile($input_folder,'agat_sp_sensitivity_specificity/query2.gff3') .
        " -o $outtmp"
    );
    check_diff($outtmp, $result, 'output sensitivity_specificity case4');
}

{
    my $dir = setup_tempdir();
    my $outtmp = catfile($dir,'tmp5.txt');
    my $result = "$output_folder/agat_sp_sensitivity_specificity_5.txt";
    check_quiet_run(
        "$script --gff1 " . catfile($input_folder,'agat_sp_sensitivity_specificity/ref3.gff3') .
        " --gff2 " . catfile($input_folder,'agat_sp_sensitivity_specificity/query3.gff3') .
        " -o $outtmp"
    );
    check_diff($outtmp, $result, 'output sensitivity_specificity case5');
}

done_testing();