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
                               --no-log \\
                               --no-create_l3_for_l2_orphan  \\
                               --throw_fasta  \\
                               --cpu 2 \\
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
                               --clean_attributes_from_template \\
                               --output $config 2>&1 1>/dev/null");
    check_diff( $config, $correct_output, 'modif config check' );
}

# ----- Test Rename config file  output ----
{
    my $dir = setup_tempdir();
    my $new_config_name = catfile( $dir, 'agat_config_renamed.yml' );
    system("$script_agat config -e --output $new_config_name --locus_tag Name 2>&1 1>/dev/null");
    ok( -e $new_config_name, 'rename agat config file check' );

# ----- Test use a selected config file ----
# ----- also test if a script uses the config file ----
    my $tmp = catfile( $dir, 'tmp.gff' );
    check_quiet_run( $script .
          " --gff " . catfile( $Bin, '..', 'gff_syntax', 'in', '28_test.gff' ) .
          " --config $new_config_name -o $tmp"
    );
    check_diff(
        $tmp,
        catfile( $Bin, '..', 'gff_syntax', 'out', '28_correct_output.gff' ),
        'Use custom agat config file check'
    );
}

# Test CLI overrides config
{
    my $dir     = setup_tempdir();
    my $config  = catfile( $dir, 'agat_config_cli_override.yaml' );
    system("$script_agat config -e --output $config --verbose 0 2>&1 1>/dev/null");
    ok( -e $config, 'config generated for CLI override test' );
    my $cfg = LoadFile($config);
    is( $cfg->{verbose} + 0, 0, 'CLI --verbose 0 overrides default verbose' );
}

# Test default log existence when running a simple script
# in input is 1.gff called the log output is called agat_log_1
{   
    my $dir     = setup_tempdir();
    my $tmp = catfile( $dir, 'tmp.gff' );
    check_quiet_run( $script . " --gff " . catfile( $Bin, '..', 'gff_syntax', 'in', '28_test.gff' ) .
                     " -o $tmp" );
    my $log_file = catfile( $dir, 'agat_log_28_test' );
    ok( -e $log_file, 'log file created at default location by script' );
}

# Test no log existence when running a simple script
{
    my $dir     = setup_tempdir();
    my $tmp = catfile( $dir, 'tmp.gff' );
    my $config  = catfile( $dir, 'agat_config_nolog.yaml' );
    system("$script_agat config -e --no-log --output $config 2>&1 1>/dev/null");
    check_quiet_run( $script . " --gff " . catfile( $Bin, '..', 'gff_syntax', 'in', '28_test.gff' ) .
                     " --config $config -o $tmp" );
    my $log_file = catfile( $dir, 'agat_log_28_test' );
    ok( ! -e $log_file, "log file '$log_file' was NOT created (expected)" );
}

# Test CLI option --verbose 0 overrides default config value (redundant check with script)
{
    my $dir     = setup_tempdir();
    my $config  = catfile( $dir, 'agat_config_verbose0.yaml' );
    system("$script_agat config -e --output $config --verbose 0 2>&1 1>/dev/null");
    my $cfg = LoadFile($config);
    is( $cfg->{verbose} + 0, 0, 'verbose set to 0 in generated config' );

    # Ensure a script runs quietly when using verbose 0
    my $tmp = catfile( $dir, 'tmp.gff' );
    check_quiet_run( $script . " --gff " . catfile( $Bin, '..', 'gff_syntax', 'in', '28_test.gff' ) .
                     " --config $config -o $tmp" );
}

done_testing();