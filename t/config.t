#!/usr/bin/env perl

use strict;
use warnings;
use Test::More tests => 3;
use FindBin qw($Bin);
use lib "$Bin/lib";
use File::Spec::Functions qw(catdir catfile);
use Cwd qw(abs_path);
use AGAT::TestUtilities qw(setup_tempdir check_diff script_prefix);

=head1 DESCRIPTION

Test to verify other type of check. e.g. configuation call.

=cut

# Check if has to be run in Devel::Cover or not
my $script_prefix = script_prefix();

# script to call to check the parser
my $script = "";
my $root = abs_path(catdir($Bin, '..'));
my $output_folder = catdir($Bin, 'config', 'out');
my $config = 'agat_config.yaml';

# ---- test gzip file and contain fasta ----
$script = $script_prefix . catfile($root, 'bin', 'agat');
my $correct_output = "$output_folder/$config";
{
    my $dir = setup_tempdir();
    system("$script config -e \\
                                                               --no-log \\
                                                               --no-create_l3_for_l2_orphan  \\
                                                               --throw_fasta  \\
                                                               --no-tabix  \\
                                                               --output_format GTF  \\
                                                               --gff_output_version 2  \\
                                                               --gtf_output_version 2  \\
                                                               --debug  \\
                                                               --deflate_attribute  \\
                                                               --no-check_all_level1_locations  \\
                                                               --no-check_identical_isoforms  \\
                                                               --no-check_utrs  \\
                                                               --no-check_sequential  \\
                                                               --no-check_all_level3_locations  \\
                                                               --no-remove_orphan_l1  \\
                                                               --no-check_cds  \\
                                                               --no-check_l1_linked_to_l2  \\
                                                               --verbose 3  \\
                                                               --no-check_l2_linked_to_l3  \\
                                                               --force_gff_input_version 3  \\
                                                               --no-check_all_level2_locations  \\
                                                               --locus_tag test1,test2 \\
                                                               --no-progress_bar  \\
                                                               --no-check_exons  \\
                                                               --merge_loci \\
                                                               --prefix_new_id nbisTEST \\
                                                               --clean_attributes_from_template  ");
    check_diff( $config, $correct_output, 'modif config check' );
}

# ----- Test Rename config file ----
{
    my $dir = setup_tempdir();
    my $new_config_name = 'agat_config_renamed.yml';
    system("$script config -e --output $new_config_name --locus_tag Name");
    ok( -e $new_config_name, 'rename agat config file check' );
}

# ----- Test use a renamed config file ----
{
    my $dir = setup_tempdir();
    my $new_config_name = 'agat_config_renamed.yml';
    system("$script config -e --output $new_config_name --locus_tag Name");
    system(catfile($root, 'bin', 'agat_convert_sp_gxf2gxf.pl') .
           " --gff " . catfile($Bin, 'gff_syntax', 'in', '28_test.gff') .
           " -c $new_config_name -o tmp.gff  2>&1 1>/dev/null");
    check_diff( 'tmp.gff', catfile($Bin, 'gff_syntax', 'out', '28_correct_output.gff'),
                'Use custom agat config file check' );
}

