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

# -------------------------- check agat_convert_genscan2gff -------------------------

my $script = $script_prefix . catfile($bin_dir, "agat_convert_genscan2gff.pl");
my $result = "$output_folder/agat_convert_genscan2gff_1.gff";
{
    my $dir = setup_tempdir();
    my $outtmp = catfile($dir, 'tmp.gff');
    system(" $script --genscan $input_folder/test.genscan -o $outtmp 2>&1 1>/dev/null");
    #run test
    check_diff( $outtmp, $result, "output $script" );
}




done_testing();
