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

# --------check agat_sp_separate_by_record_type.pl-------------

my $script = $script_prefix . catfile($bin_dir, 'agat_sp_separate_by_record_type.pl');

{ ok(system("$script -h 1>/dev/null") == 0, "help $script"); }

{
    my $dir = setup_tempdir();
    my $outtmp = catfile($dir, 'tmp1');
    my $result = "$output_folder/agat_sp_separate_by_record_type_1.gff";
    check_quiet_run(
        "$script --gff " . catfile($input_folder, '1.gff') ." -o $outtmp"
    );
    check_diff("$outtmp/trna.gff", $result, 'output agat_sp_separate_by_record_type.pl');
}

done_testing();
