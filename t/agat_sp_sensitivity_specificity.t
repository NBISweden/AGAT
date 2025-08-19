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

# ------------------- check agat_sp_sensitivity_specificity script-------------------

my $script = $script_prefix . catfile($bin_dir, "agat_sp_sensitivity_specificity.pl");
my $result = "$output_folder/agat_sp_sensitivity_specificity_1.txt";
{
    my $dir = setup_tempdir();
    my $outtmp   = catfile( $dir, 'tmp.gff' );
    system(" $script --gff1 $input_folder/1.gff --gff2 $input_folder/1.gff -o $outtmp 2>&1 1>/dev/null");
    #run test
    check_diff( $outtmp, $result, "output $script", "-I '^usage:'" );
}


$script = $script_prefix . catfile($bin_dir, "agat_sp_sensitivity_specificity.pl");
$result = "$output_folder/agat_sp_sensitivity_specificity_2.txt";
{
    my $dir = setup_tempdir();
    my $outtmp   = catfile( $dir, 'tmp.gff' );
    system(" $script --gff1 $input_folder/agat_sp_sensitivity_specificity/ref0.gff3 --gff2 $input_folder/agat_sp_sensitivity_specificity/query0.gff3 -o $outtmp 2>&1 1>/dev/null");
    #run test
    check_diff( $outtmp, $result, "output $script", "-I '^usage:'" );
}


$script = $script_prefix . catfile($bin_dir, "agat_sp_sensitivity_specificity.pl");
$result = "$output_folder/agat_sp_sensitivity_specificity_3.txt";
{
    my $dir = setup_tempdir();
    my $outtmp   = catfile( $dir, 'tmp.gff' );
    system(" $script --gff1 $input_folder/agat_sp_sensitivity_specificity/ref1.gff3 --gff2 $input_folder/agat_sp_sensitivity_specificity/query1.gff3 -o $outtmp 2>&1 1>/dev/null");
    #run test
    check_diff( $outtmp, $result, "output $script", "-I '^usage:'" );
}


$script = $script_prefix . catfile($bin_dir, "agat_sp_sensitivity_specificity.pl");
$result = "$output_folder/agat_sp_sensitivity_specificity_4.txt";
{
    my $dir = setup_tempdir();
    my $outtmp   = catfile( $dir, 'tmp.gff' );
    system(" $script --gff1 $input_folder/agat_sp_sensitivity_specificity/ref2.gff3 --gff2 $input_folder/agat_sp_sensitivity_specificity/query2.gff3 -o $outtmp 2>&1 1>/dev/null");
    #run test
    check_diff( $outtmp, $result, "output $script", "-I '^usage:'" );
}


$script = $script_prefix . catfile($bin_dir, "agat_sp_sensitivity_specificity.pl");
$result = "$output_folder/agat_sp_sensitivity_specificity_5.txt";
{
    my $dir = setup_tempdir();
    my $outtmp   = catfile( $dir, 'tmp.gff' );
    system(" $script --gff1 $input_folder/agat_sp_sensitivity_specificity/ref3.gff3 --gff2 $input_folder/agat_sp_sensitivity_specificity/query3.gff3 -o $outtmp 2>&1 1>/dev/null");
    #run test
    check_diff( $outtmp, $result, "output $script", "-I '^usage:'" );
}



done_testing();
