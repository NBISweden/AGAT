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
my $script = $script_prefix . catfile($root, 'bin', 'agat_convert_sp_gff2gtf.pl');
my $input_folder = catdir($Bin, 'gff_other', 'in');
my $output_folder = catdir($Bin, 'gff_other', 'out');

# --------- Issue 457 multi-values attributes (gene_name "26266" "MT-TL1";) can be deflated to be compliant with GTF and CellRanger --------
{
    my $dir = setup_tempdir();
    my $pathtmp = catfile($dir, 'tmp.gff');
    my $correct_output = catfile($output_folder, 'issue457.gtf');
    system("$script_agat config --expose --deflate_attribute 2>&1 1>/dev/null");
    system("$script --gff " . catfile($input_folder, 'issue457.gff') . " -o $pathtmp  2>&1 1>/dev/null");
    check_diff($pathtmp, $correct_output, 'issue457 check');
}

done_testing();
