#!/usr/bin/env perl
use strict;
use warnings;
use Test::More;
use FindBin qw($Bin);
use File::Spec::Functions qw(catdir catfile);
use Cwd qw(abs_path);
use lib catdir($Bin, '..', 'lib');
use AGAT::TestUtilities qw(setup_tempdir check_diff script_prefix);

my $script_prefix = script_prefix();
my $root = abs_path(catdir($Bin, '..', '..'));
my $script_agat = $script_prefix . catfile($root, 'bin', 'agat');
my $script = $script_prefix . catfile($root, 'bin', 'agat_convert_sp_gxf2gxf.pl');
my $input_folder = catdir($Bin, 'level_missing', 'in');
my $output_folder = catdir($Bin, 'level_missing', 'out');

# -------------------------- testD -------------------------
{
    my $dir = setup_tempdir();
    my $outtmp = catfile($dir, 'tmp.gff');
    my $result = catfile($output_folder, 'testD_output.gff');
    system("$script --gff " . catfile($input_folder, 'testD.gff') . " -o $outtmp 2>&1 1>/dev/null");
    check_diff($outtmp, $result, 'output testD');
}
{
    my $dir = setup_tempdir();
    my $outtmp = catfile($dir, 'tmp.gff');
    my $result = catfile($output_folder, 'testD_output2.gff');
    system("$script_agat config --expose --locus_tag ID 2>&1 1>/dev/null");
    system("$script --gff " . catfile($input_folder, 'testD.gff') . " -o $outtmp 2>&1 1>/dev/null");
    check_diff($outtmp, $result, 'output testD2');
}

done_testing();
