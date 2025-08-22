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
my $script = $script_prefix . catfile($root, 'bin', 'agat_convert_sp_gxf2gxf.pl');
my $input_folder = catdir($Bin, 'gff_other', 'in');
my $output_folder = catdir($Bin, 'gff_other', 'out');

# -------- test gzip file and contain fasta --------
{
    my $dir = setup_tempdir();
    my $pathtmp = catfile($dir, 'tmp.gff');
    my $correct_output = catfile($output_folder, 'zip_and_fasta_correct_output.gff');
    check_quiet_run("$script --gff " . catfile($input_folder, 'zip_and_fasta.gff.gz') . " -o $pathtmp");
    check_diff($pathtmp, $correct_output, 'zip_and_fasta check');
}

done_testing();
