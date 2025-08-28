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
my $input_folder = catdir($Bin, 'gff_other', 'in');
my $output_folder = catdir($Bin, 'gff_other', 'out');

# --------- Issue 448 bioperl adding extra empty ID attribute that mess up AGAT (only when input parsed with version 2 and 2.5) ----
{
    my $dir = setup_tempdir();
    my $pathtmp = catfile($dir, 'tmp.gff');
    my $correct_output = catfile($output_folder, 'issue448.gtf');
    check_quiet_run("$script_agat config --expose --output_format gtf");
    check_quiet_run("$script --g " . catfile($input_folder, 'issue448.gtf') . " -o $pathtmp");
    check_diff($pathtmp, $correct_output, 'issue448 check');
}
{
    my $dir = setup_tempdir();
    my $pathtmp = catfile($dir, 'tmp.gff');
    my $correct_output = catfile($output_folder, 'issue448.gff');
    check_quiet_run("$script --g " . catfile($input_folder, 'issue448.gtf') . " -o $pathtmp");
    check_diff($pathtmp, $correct_output, 'issue448 check');
}

done_testing();
