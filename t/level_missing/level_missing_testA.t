#!/usr/bin/env perl
use strict;
use warnings;
use Test::More;
use FindBin qw($Bin);
use File::Spec::Functions qw(catdir catfile);
use Cwd qw(abs_path);
use lib catdir($Bin, '..', 'lib');
use AGAT::TestUtilities qw(setup_tempdir check_diff script_prefix check_quiet_run);

my $script_prefix = script_prefix();
my $root = abs_path(catdir($Bin, '..', '..'));
my $script_agat = $script_prefix . catfile($root, 'bin', 'agat');
my $script = $script_prefix . catfile($root, 'bin', 'agat_convert_sp_gxf2gxf.pl');
my $input_folder = catdir($Bin, 'in');
my $output_folder = catdir($Bin, 'out');

# -------------------------- testA -------------------------
{
    my $dir = setup_tempdir();
    my $outtmp = catfile($dir, 'tmp.gff');
    my $result = catfile($output_folder, 'testA_output.gff');
    check_quiet_run("$script --gff " . catfile($input_folder, 'testA.gff') . " -o $outtmp");
    check_diff($outtmp, $result, 'output testA');
}
{
    my $dir = setup_tempdir();
    my $outtmp = catfile($dir, 'tmp.gff');
    my $result = catfile($output_folder, 'testA_output2.gff');
    system("$script_agat config --expose --locus_tag common_tag 2>&1 1>/dev/null");
    check_quiet_run("$script --gff " . catfile($input_folder, 'testA.gff') . " -o $outtmp");
    check_diff($outtmp, $result, 'output testA2');
}
{
    my $dir = setup_tempdir();
    my $outtmp = catfile($dir, 'tmp.gff');
    my $result = catfile($output_folder, 'testA_output3.gff');
    system("$script_agat config --expose --locus_tag gene_info 2>&1 1>/dev/null");
    check_quiet_run("$script --gff " . catfile($input_folder, 'testA.gff') . " -o $outtmp");
    check_diff($outtmp, $result, 'output testA3');
}
{
    my $dir = setup_tempdir();
    my $outtmp = catfile($dir, 'tmp.gff');
    my $result = catfile($output_folder, 'testA_output4.gff');
    system("$script_agat config --expose --locus_tag transcript_id 2>&1 1>/dev/null");
    check_quiet_run("$script --gff " . catfile($input_folder, 'testA.gff') . " -o $outtmp");
    check_diff($outtmp, $result, 'output testA4');
}

done_testing();
