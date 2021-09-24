#!/usr/bin/env perl

use strict;
use warnings;
use Test::More tests => 6;
use AGAT::Omniscient;
use AGAT::OmniscientTool;

=head1 DESCRIPTION

Test to verify independent functions

=cut

# Check if has to be run in Devel::Cover or not
my $script_prefix="";
if (exists $ENV{'HARNESS_PERL_SWITCHES'} ) {
  if ($ENV{'HARNESS_PERL_SWITCHES'} =~ m/Devel::Cover/) {
    $script_prefix="perl -MDevel::Cover ";
  }
}


my ($hash_omniscient, $hash_mRNAGeneLink) = slurp_gff3_file_JD({ input => "t/scripts_output/1.gff",
                                                                  verbose => 0});




#run remove_l2_related_feature test
my $feature2 = @{$hash_omniscient->{"level2"}{"mrna"}{"gene:os01g0100100"}}[0];
my $nb_gene = scalar keys %{$hash_omniscient->{"level1"}{"gene"}};
remove_l2_related_feature($hash_omniscient,$feature2, 0);
my $nb_gene2 = scalar keys %{$hash_omniscient->{"level1"}{"gene"}};
ok(  $nb_gene != $nb_gene2, "remove_l2_related_feature");

#run remove_omniscient_elements_from_level1_id_list test
remove_omniscient_elements_from_level1_id_list($hash_omniscient,["gene:os01g0100200"]);
my $nb_gene3 = scalar keys %{$hash_omniscient->{"level1"}{"gene"}};
ok(  $nb_gene2 != $nb_gene3, "remove_omniscient_elements_from_level1_id_list");

 
#run remove_omniscient_elements_from_level2_feature_list test
$feature2 = @{$hash_omniscient->{"level2"}{"mrna"}{"gene:os01g0100300"}}[0];
remove_omniscient_elements_from_level2_feature_list($hash_omniscient, [$feature2]);
my $nb_gene4 = scalar keys %{$hash_omniscient->{"level1"}{"gene"}};
ok(  $nb_gene3 != $nb_gene4, "remove_omniscient_elements_from_level2_feature_list");

# run group_features_from_omniscient test
ok(  group_features_from_omniscient($hash_omniscient), "group_features_from_omniscient");

  
# run group_l1IDs_from_omniscient test
ok( group_l1IDs_from_omniscient($hash_omniscient), "group_l1IDs_from_omniscient");

# run get_longest_cds_start_end
ok( get_longest_cds_start_end($hash_omniscient,"gene:os01g0100466"), "get_longest_cds_start_end");

