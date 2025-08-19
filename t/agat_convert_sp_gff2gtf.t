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

# -------------------------- check agat_convert_sp_gff2gtf -------------------------
my $convert_sp_gff2gtf_folder = "$input_folder/agat_convert_sp_gff2gtf";

my $script = $script_prefix . catfile($bin_dir, "agat_convert_sp_gff2gtf.pl");
my $result = "$convert_sp_gff2gtf_folder/agat_convert_sp_gff2gtf_1.gtf";
{
    my $dir = setup_tempdir();
    my $outtmp = catfile($dir, 'tmp.gff');
    my $outprefix = catfile($dir, 'tmp');
    system(" $script --gff $input_folder/1.gff --gtf_version 3 -o $outtmp 2>&1 1>/dev/null");
    #run test
    check_diff( $outtmp, $result, "output $script" );
}


$result = "$convert_sp_gff2gtf_folder/agat_convert_sp_gff2gtf_2.gtf";
{
    my $dir = setup_tempdir();
    my $outtmp = catfile($dir, 'tmp.gff');
    my $outprefix = catfile($dir, 'tmp');
    system(" $script --gff $convert_sp_gff2gtf_folder/stop_start_an_exon.gff --gtf_version 3 -o $outtmp 2>&1 1>/dev/null");
    #run test
    check_diff( $outtmp, $result, "output $script" );
}


$result = "$convert_sp_gff2gtf_folder/agat_convert_sp_gff2gtf_3.gtf";
{
    my $dir = setup_tempdir();
    my $outtmp = catfile($dir, 'tmp.gff');
    my $outprefix = catfile($dir, 'tmp');
    system(" $script --gff $convert_sp_gff2gtf_folder/stop_split_over_two_exons.gff --gtf_version 3 -o $outtmp 2>&1 1>/dev/null");
    #run test
    check_diff( $outtmp, $result, "output $script" );
}


$result = "$convert_sp_gff2gtf_folder/result_issue_245.gtf";
{
    my $dir = setup_tempdir();
    my $outtmp = catfile($dir, 'tmp.gff');
    my $outprefix = catfile($dir, 'tmp');
    system(" $script --gff $convert_sp_gff2gtf_folder/issue_245.gff --gtf_version 3 -o $outtmp 2>&1 1>/dev/null");
    #run test
    check_diff( $outtmp, $result, "output $script" );
}




done_testing();
