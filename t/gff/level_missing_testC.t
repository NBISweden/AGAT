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
my $input_folder = catdir($Bin, 'level_missing', 'in');
my $output_folder = catdir($Bin, 'level_missing', 'out');

# -------------------------- testC -------------------------
{
    my $dir = setup_tempdir();
    my $outtmp = catfile($dir, 'tmp.gff');
    my $result = catfile($output_folder, 'testC_output.gff');
    check_quiet_run("$script --gff " . catfile($input_folder, 'testC.gff') . " -o $outtmp");
    check_diff($outtmp, $result, 'output testC');
}
{
    my $dir = setup_tempdir();
    my $outtmp = catfile($dir, 'tmp.gff');
    my $result = catfile($output_folder, 'testC_output2.gff');
    check_quiet_run("$script_agat config --expose --locus_tag locus_id");
    check_quiet_run("$script --gff " . catfile($input_folder, 'testC.gff') . " -o $outtmp");
    check_diff($outtmp, $result, 'output testC2');
}

done_testing();
