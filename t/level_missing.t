#!/usr/bin/env perl

use strict;
use warnings;
use Test::More tests => 13;
use FindBin qw($Bin);
use lib "$Bin/lib";
use File::Spec::Functions qw(catdir catfile);
use Cwd qw(abs_path);
use AGAT::TestUtilities qw(setup_tempdir check_diff script_prefix);

=head1 DESCRIPTION

Test to verify the parser deals properly with the different flavor / bugged gff files.

=cut

# Check if has to be run in Devel::Cover or not
my $script_prefix = script_prefix();

# script to call to check the parser
my $root = abs_path(catdir($Bin, '..'));
my $script_agat = $script_prefix . catfile($root, 'bin', 'agat');
my $script = $script_prefix . catfile($root, 'bin', 'agat_convert_sp_gxf2gxf.pl');
my $input_folder = catdir($Bin, 'level_missing', 'in');
my $output_folder = catdir($Bin, 'level_missing', 'out');

# -------------------------- testA -------------------------

{
    my $dir = setup_tempdir();
    my $outtmp = catfile( $dir, 'tmp.gff' );
    my $result = catfile( $output_folder, 'testA_output.gff' );
    system("$script --gff " . catfile( $input_folder, 'testA.gff' ) . " -o $outtmp 2>&1 1>/dev/null");
    check_diff( $outtmp, $result, "output testA" );
}

{
    my $dir = setup_tempdir();
    my $outtmp = catfile( $dir, 'tmp.gff' );
    my $result = catfile( $output_folder, 'testA_output2.gff' );
    system("$script_agat config --expose --locus_tag common_tag 2>&1 1>/dev/null");
    system("$script --gff " . catfile( $input_folder, 'testA.gff' ) . " -o $outtmp 2>&1 1>/dev/null");
    check_diff( $outtmp, $result, "output testA2" );
}

{
    my $dir = setup_tempdir();
    my $outtmp = catfile( $dir, 'tmp.gff' );
    my $result = catfile( $output_folder, 'testA_output3.gff' );
    system("$script_agat config --expose --locus_tag gene_info 2>&1 1>/dev/null");
    system("$script --gff " . catfile( $input_folder, 'testA.gff' ) . " -o $outtmp 2>&1 1>/dev/null");
    check_diff( $outtmp, $result, "output testA3" );
}

{
    my $dir = setup_tempdir();
    my $outtmp = catfile( $dir, 'tmp.gff' );
    my $result = catfile( $output_folder, 'testA_output4.gff' );
    system("$script_agat config --expose --locus_tag transcript_id 2>&1 1>/dev/null");
    system("$script --gff " . catfile( $input_folder, 'testA.gff' ) . " -o $outtmp 2>&1 1>/dev/null");
    check_diff( $outtmp, $result, "output testA4" );
}

# -------------------------- testB -------------------------

{
    my $dir = setup_tempdir();
    my $outtmp = catfile( $dir, 'tmp.gff' );
    my $result = catfile( $output_folder, 'testB_output.gff' );
    system("$script --gff " . catfile( $input_folder, 'testB.gff' ) . " -o $outtmp 2>&1 1>/dev/null");
    check_diff( $outtmp, $result, "output testB" );
}

{
    my $dir = setup_tempdir();
    my $outtmp = catfile( $dir, 'tmp.gff' );
    my $result = catfile( $output_folder, 'testB_output2.gff' );
    system("$script_agat config --expose --locus_tag locus_id 2>&1 1>/dev/null");
    system("$script --gff " . catfile( $input_folder, 'testB.gff' ) . " -o $outtmp 2>&1 1>/dev/null");
    check_diff( $outtmp, $result, "output testB2" );
}

# -------------------------- testC -------------------------

{
    my $dir = setup_tempdir();
    my $outtmp = catfile( $dir, 'tmp.gff' );
    my $result = catfile( $output_folder, 'testC_output.gff' );
    system("$script --gff " . catfile( $input_folder, 'testC.gff' ) . " -o $outtmp 2>&1 1>/dev/null");
    check_diff( $outtmp, $result, "output testC" );
}

{
    my $dir = setup_tempdir();
    my $outtmp = catfile( $dir, 'tmp.gff' );
    my $result = catfile( $output_folder, 'testC_output2.gff' );
    system("$script_agat config --expose --locus_tag locus_id 2>&1 1>/dev/null");
    system("$script --gff " . catfile( $input_folder, 'testC.gff' ) . " -o $outtmp 2>&1 1>/dev/null");
    check_diff( $outtmp, $result, "output testC2" );
}

# -------------------------- testD -------------------------

{
    my $dir = setup_tempdir();
    my $outtmp = catfile( $dir, 'tmp.gff' );
    my $result = catfile( $output_folder, 'testD_output.gff' );
    system("$script --gff " . catfile( $input_folder, 'testD.gff' ) . " -o $outtmp 2>&1 1>/dev/null");
    check_diff( $outtmp, $result, "output testD" );
}

{
    my $dir = setup_tempdir();
    my $outtmp = catfile( $dir, 'tmp.gff' );
    my $result = catfile( $output_folder, 'testD_output2.gff' );
    system("$script_agat config --expose --locus_tag ID 2>&1 1>/dev/null");
    system("$script --gff " . catfile( $input_folder, 'testD.gff' ) . " -o $outtmp 2>&1 1>/dev/null");
    check_diff( $outtmp, $result, "output testD2" );
}

# -------------------------- testE -------------------------

{
    my $dir = setup_tempdir();
    my $outtmp = catfile( $dir, 'tmp.gff' );
    my $result = catfile( $output_folder, 'testE_output.gff' );
    system("$script --gff " . catfile( $input_folder, 'testE.gff' ) . " -o $outtmp 2>&1 1>/dev/null");
    check_diff( $outtmp, $result, "output testE" );
}

# -------------------------- testF -------------------------

{
    my $dir = setup_tempdir();
    my $outtmp = catfile( $dir, 'tmp.gff' );
    my $result = catfile( $output_folder, 'testF_output.gff' );
    system("$script --gff " . catfile( $input_folder, 'testF.gff' ) . " -o $outtmp 2>&1 1>/dev/null");
    check_diff( $outtmp, $result, "output testF" );
}

# -------------------------- testG -------------------------

{
    my $dir = setup_tempdir();
    my $outtmp = catfile( $dir, 'tmp.gff' );
    my $result = catfile( $output_folder, 'testG_output.gff' );
    system("$script --gff " . catfile( $input_folder, 'testG.gff' ) . " -o $outtmp 2>&1 1>/dev/null");
    check_diff( $outtmp, $result, "output testG" );
}
