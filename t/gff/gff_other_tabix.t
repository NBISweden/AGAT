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
my $input_folder = catdir($Bin, 'gff_other', 'in');
my $output_folder = catdir($Bin, 'gff_other', 'out');

# -------- test tabix output sorting --------
{
    my $dir = setup_tempdir();
    my $pathtmp = catfile($dir, 'tmp.gff');
    my $correct_output = catfile($output_folder, '1_agat_tabix.gff');
    system("$script_agat config --expose --tabix 2>&1 1>/dev/null");
    system("$script --gff " . catfile($Bin, '..', 'scripts_output', 'in', '1.gff') . " -o $pathtmp 2>&1 1>/dev/null");
    check_diff($pathtmp, $correct_output, 'tabix check');
}

done_testing();
