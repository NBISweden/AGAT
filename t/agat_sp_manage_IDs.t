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

# --------check agat_sp_manage_IDs.pl-------------

my $script = $script_prefix . catfile($bin_dir, "agat_sp_manage_IDs.pl");
{ my $dir = setup_tempdir(); ok(system("$script -h 1>\/dev\/null") == 0, "help $script"); }

{
    my $dir = setup_tempdir();
    my $outtmp = catfile( $dir, 'tmp.gff' );
    ok( system(" $script --gff $input_folder/1.gff --ensembl -o $outtmp 2>&1 1>/dev/null") == 0,
        "output $script" );
}

{
    my $dir = setup_tempdir();
    ok( system(" $script --gff $input_folder/1.gff --tair --ensembl 2>/dev/null") != 0,
        "conflicting style flags" );
}

{
    my $dir = setup_tempdir();
    ok( system(" $script --gff $input_folder/1.gff --nb -1 2>/dev/null") != 0,
        "negative nb rejected" );
}

{
    my $dir = setup_tempdir();
    ok( system(" $script --gff $input_folder/1.gff --gap -2 2>/dev/null") != 0,
        "negative gap rejected" );
}

done_testing();
