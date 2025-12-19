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

# -------------------------- check agat_convert_sp_gff2gtf -------------------------

my $convert_sp_gff2gtf_out_folder = "$output_folder/agat_convert_sp_gff2gtf";
my $convert_sp_gff2gtf_in_folder = "$input_folder/agat_convert_sp_gff2gtf";
my $script = $script_prefix . catfile($bin_dir, "agat_convert_sp_gff2gtf.pl");

{ ok(system("$script -h 1>\/dev\/null") == 0, "help $script"); }

{
    my $dir = setup_tempdir();
    my $outtmp = catfile($dir,'tmp1.gtf');
    my $result = "$convert_sp_gff2gtf_out_folder/agat_convert_sp_gff2gtf_1.gtf";
    check_quiet_run(
        "$script --gff " . catfile($input_folder,'1.gff') . " --gtf_version 3 -o $outtmp"
    );
    check_diff($outtmp, $result, 'gff2gtf case1');
}

{
    my $dir = setup_tempdir();
    my $outtmp = catfile($dir,'tmp2.gtf');
    my $result = "$convert_sp_gff2gtf_out_folder/agat_convert_sp_gff2gtf_2.gtf";
    check_quiet_run(
        "$script --gff " . catfile($convert_sp_gff2gtf_in_folder,'stop_start_an_exon.gff') . " --gtf_version 3 -o $outtmp"
    );
    check_diff($outtmp, $result, 'gff2gtf case2');
}

{
    my $dir = setup_tempdir();
    my $outtmp = catfile($dir,'tmp3.gtf');
    my $result = "$convert_sp_gff2gtf_out_folder/agat_convert_sp_gff2gtf_3.gtf";
    check_quiet_run(
        "$script --gff " . catfile($convert_sp_gff2gtf_in_folder,'stop_split_over_two_exons.gff') . " --gtf_version 3 -o $outtmp"
    );
    check_diff($outtmp, $result, 'gff2gtf case3');
}

{
    my $dir = setup_tempdir();
    my $outtmp = catfile($dir,'tmp4.gtf');
    my $result = "$convert_sp_gff2gtf_out_folder/result_issue_245.gtf";
    check_quiet_run(
        "$script --gff " . catfile($convert_sp_gff2gtf_in_folder,'issue_245.gff') . " --gtf_version 3 -o $outtmp"
    );
    check_diff($outtmp, $result, 'gff2gtf issue245');
}

done_testing();
