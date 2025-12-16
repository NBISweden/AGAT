#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use File::Spec::Functions qw(catdir catfile);
use Cwd qw(abs_path);
use AGAT::TestUtilities; 
use Test::More;

my $script_prefix = script_prefix();
my $root = abs_path(catdir($Bin, '..', '..'));
my $bin_dir = catdir($root, 'bin');
my $input_folder = catdir($Bin, 'in');
my $output_folder = catdir($Bin, 'out');

# -------------------------- check agat_convert_sp_gff2tsv -------------------------

my $script = $script_prefix . catfile($bin_dir, 'agat_convert_sp_gff2tsv.pl');

{ ok(system("$script -h 1>/dev/null") == 0, "help $script"); }

{
    my $dir = setup_tempdir();
    my $outtmp = catfile($dir,'tmp.tsv');
    my $result = "$output_folder/agat_convert_sp_gff2tsv_1.tsv";
    check_quiet_run(
        "$script --gff " . catfile($input_folder,'1.gff') . " -o $outtmp"
    );
    check_diff($outtmp, $result, 'gff2tsv');
}

done_testing();
