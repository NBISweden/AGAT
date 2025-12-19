#!/usr/bin/env perl

use strict;
use warnings;
use Test::More;
use FindBin qw($Bin);
use File::Spec::Functions qw(catdir catfile);
use lib catdir($Bin, '..', 'lib');
use Cwd qw(abs_path);
use AGAT::TestUtilities qw(setup_tempdir check_diff script_prefix check_quiet_run);
use YAML qw(LoadFile);

=head1 DESCRIPTION

Test to verify other type of check. e.g. configuration call.

=cut

# Check if has to be run in Devel::Cover or not
my $script_prefix = script_prefix();
my $root = abs_path(catdir($Bin, '..', '..'));
my $script_agat = $script_prefix . catfile($root, 'bin', 'agat');
my $script = $script_prefix . catfile($root, 'bin', 'agat_convert_sp_gxf2gxf.pl');
my $output_folder = catdir($Bin, 'out');

# ---- test expose config and change value ----
{
    my $dir     = setup_tempdir();
    my $config  = catfile( $dir, 'agat_config.yaml' );
    my $correct_output = catfile( $output_folder, 'agat_config.yaml' );
    system("$script_agat config -e \\
                               --no-check_all_level1_locations  \\
                               --no-check_all_level2_locations  \\
                               --no-check_all_level3_locations  \\
                               --no-check_cds  \\
                               --no-check_exons  \\
                               --no-check_identical_isoforms  \\
                               --no-check_l1_linked_to_l2  \\
                               --no-check_l2_linked_to_l3  \\
                               --no-check_sequential  \\
                               --no-check_utrs  \\
                               --clean_attributes_from_template \\
                               --cpu 2 \\
                               --no-create_l3_for_l2_orphan  \\
                               --debug  \\
                               --deflate_attribute  \\
                               --force  \\
                               --force_gff_input_version 3  \\
                               --gff_output_version 2  \\
                               --gtf_output_version 2  \\
                               --locus_tag test1,test2 \\
                               --no-log \\
                               --merge_loci \\
                               --output_format GTF  \\
                               --prefix_new_id nbisTEST \\
                               --no-progress_bar  \\
                               --no-remove_orphan_l1  \\
                               --no-tabix  \\
                               --throw_fasta  \\
                               --no-url_encode_out \\
                               --verbose 3  \\
                               --output $config 2>&1 1>/dev/null");
    check_diff( $config, $correct_output, 'modif config check' );
}


done_testing();