#!/usr/bin/env perl

use strict;
use warnings;
use Test::More tests => 1;
use Bio::Tools::GFF;
use AGAT::Omniscient;
use AGAT::OmniscientTool;

=head1 DESCRIPTION

Test to verify other type of check. e.g. configuation call.

=cut

# Check if has to be run in Devel::Cover or not
my $script_prefix="";
if (exists $ENV{'HARNESS_PERL_SWITCHES'} ) {
  if ($ENV{'HARNESS_PERL_SWITCHES'} =~ m/Devel::Cover/) {
    $script_prefix="perl -MDevel::Cover ";
  }
}

# script to call to check the parser
my $script = "";
my $output_folder = "t/other";

# ---- test gzip file and contain fasta ----
$script = $script_prefix."bin/agat";
my $correct_output = "$output_folder/config.yaml";

system("$script config -e \\
								--no-log \\
								--no-create_l3_for_l2_orphan  \\
								--throw_fasta  \\
								--no-tabix  \\
								--gff_output_version 2  \\
								--debug  \\
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
								--prefix_new_id nbisTEST   ");

#run test
ok( system("diff config.yaml $correct_output") == 0, "modif config check");
# remove file created for the test
unlink "config.yaml";
