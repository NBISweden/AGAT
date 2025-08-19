#!/usr/bin/env perl

use strict;
use warnings;
use Test::More tests => 9;
use FindBin qw($Bin);
use lib "$Bin/lib";
use File::Spec::Functions qw(catdir catfile);
use Cwd qw(abs_path);
use AGAT::TestUtilities qw(setup_tempdir check_diff script_prefix);

=head1 DESCRIPTION

Test to verify other GFF/GTF peculiarities e.g HEADER / FASTA / GZIP

=cut

# Check if has to be run in Devel::Cover or not
my $script_prefix = script_prefix();

# script to call to check the parser
my $script = "";
my $root = abs_path(catdir($Bin, '..'));
my $script_agat = $script_prefix . catfile($root, 'bin', 'agat');
my $input_folder = catdir($Bin, 'gff_other', 'in');
my $output_folder = catdir($Bin, 'gff_other', 'out');
# remove config in local folder if exists and potential tmp file already existing

# -------- test gzip file and contain fasta --------
{
    my $dir = setup_tempdir();
    my $pathtmp = catfile($dir, 'tmp.gff');
    $script = $script_prefix . catfile($root, 'bin', 'agat_convert_sp_gxf2gxf.pl');
    my $correct_output = "$output_folder/zip_and_fasta_correct_output.gff";

    system("$script --gff " . catfile($input_folder, 'zip_and_fasta.gff.gz') . " -o $pathtmp  2>&1 1>/dev/null");

    check_diff( $pathtmp, $correct_output, "zip_and_fasta check" );
}

# -------- test tabix output sorting --------
{
    my $dir = setup_tempdir();
    my $pathtmp = catfile($dir, 'tmp.gff');
    $script = $script_prefix . catfile($root, 'bin', 'agat_convert_sp_gxf2gxf.pl');
    my $correct_output = "$output_folder/1_agat_tabix.gff";

    system("$script_agat config --expose --tabix 2>&1 1>/dev/null");
    system("$script --gff " . catfile($Bin, 'scripts_output', 'in', '1.gff') . " -o $pathtmp 2>&1 1>/dev/null");

    check_diff( $pathtmp, $correct_output, "tabix check" );
}

# -------- Parent ID already used by same level feature --------
{
    my $dir = setup_tempdir();
    my $pathtmp = catfile($dir, 'tmp.gff');
    $script = $script_prefix . catfile($root, 'bin', 'agat_convert_sp_gxf2gxf.pl');
    my $correct_output = "$output_folder/issue329.gff";

    system("$script --gff " . catfile($input_folder, 'issue329.gff') . " -o $pathtmp 2>&1 1>/dev/null");

    check_diff( $pathtmp, $correct_output, "issue329 check" );
}

# -------- Issue 368 avoid seq_id 0 replaced by SEQ --------
{
    my $dir = setup_tempdir();
    my $pathtmp = catfile($dir, 'tmp.gff');
    $script = $script_prefix . catfile($root, 'bin', 'agat_convert_sp_gxf2gxf.pl');
    my $correct_output = "$output_folder/issue368.gff";

    system("$script --gff " . catfile($input_folder, 'issue368.gff') . " -o $pathtmp 2>&1 1>/dev/null");

    check_diff( $pathtmp, $correct_output, "issue368 check" );
}

# -------- Issue 389 very wierd --------
{
    my $dir = setup_tempdir();
    my $pathtmp = catfile($dir, 'tmp.gff');
    $script = $script_prefix . catfile($root, 'bin', 'agat_convert_sp_gxf2gxf.pl');
    my $correct_output = "$output_folder/issue389.gff";

    system("$script --gff " . catfile($input_folder, 'issue389.gff') . " -o $pathtmp 2>&1 1>/dev/null");

    check_diff( $pathtmp, $correct_output, "issue389 check" );
}

# --------- Issue 441 transccipt_id used for GTF while L2 L1 created from L3 with isoforms ----
{
    my $dir = setup_tempdir();
    my $pathtmp = catfile($dir, 'tmp.gff');
    $script = $script_prefix . catfile($root, 'bin', 'agat_convert_sp_gxf2gxf.pl');
    my $correct_output = "$output_folder/issue441.gtf";

    system("$script_agat config --expose --output_format gtf 2>&1 1>/dev/null");
    system("$script --g " . catfile($input_folder, 'issue441.gtf') . " -o $pathtmp  2>&1 1>/dev/null");

    check_diff( $pathtmp, $correct_output, "issue441 check" );
}

# --------- Issue 448 bioperl adding extra empty ID attribute that mess up AGAT (only when input parsed with version 2 and 2.5)  ----
{
    my $dir = setup_tempdir();
    my $pathtmp = catfile($dir, 'tmp.gff');
    $script = $script_prefix . catfile($root, 'bin', 'agat_convert_sp_gxf2gxf.pl');
    my $correct_output = "$output_folder/issue448.gtf";

    system("$script_agat config --expose --output_format gtf 2>&1 1>/dev/null");
    system("$script --g " . catfile($input_folder, 'issue448.gtf') . " -o $pathtmp  2>&1 1>/dev/null");

    check_diff( $pathtmp, $correct_output, "issue448 check" );
}

{
    my $dir = setup_tempdir();
    my $pathtmp = catfile($dir, 'tmp.gff');
    $script = $script_prefix . catfile($root, 'bin', 'agat_convert_sp_gxf2gxf.pl');
    my $correct_output = "$output_folder/issue448.gff";

    system("$script --g " . catfile($input_folder, 'issue448.gtf') . " -o $pathtmp  2>&1 1>/dev/null");

    check_diff( $pathtmp, $correct_output, "issue448 check" );
}

# --------- Issue 457 multi-values attributes (gene_name "26266" "MT-TL1";) can be deflated to be compliant with GTF and CellRanger
{
    my $dir = setup_tempdir();
    my $pathtmp = catfile($dir, 'tmp.gff');
    $script = $script_prefix . catfile($root, 'bin', 'agat_convert_sp_gff2gtf.pl');
    my $correct_output = "$output_folder/issue457.gtf";

    system("$script_agat config --expose --deflate_attribute 2>&1 1>/dev/null");
    system("$script --gff " . catfile($input_folder, 'issue457.gff') . " -o $pathtmp  2>&1 1>/dev/null");

    check_diff( $pathtmp, $correct_output, "issue457 check" );
}
