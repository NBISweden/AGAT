#!/usr/bin/env perl

use strict;
use warnings;
use Test::More tests => 16;
use Bio::Tools::GFF;
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


# remove config in local folder if exists
unlink "config.yaml"; 

# get standard config
my $config = get_agat_config();
$config->{verbose}=0;
$config->{progress_bar}=0;

# Get one omnisceint to work with
my ($hash_omniscient, $hash_mRNAGeneLink) = slurp_gff3_file_JD({ input => "t/scripts_output/1.gff",
                                                                config => $config
                                                                });

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


# run create_omniscient
my $feature1 = $hash_omniscient->{"level1"}{"gene"}{"gene:os01g0100500"};
$feature2 = @{$hash_omniscient->{"level2"}{"mrna"}{"gene:os01g0100500"}}[0];
my $l2_ID = lc($feature2->_tag_value("ID"));
my $feature3 = @{$hash_omniscient->{"level3"}{"cds"}{$l2_ID}}[0];
ok( create_omniscient([$feature1], [$feature2], [$feature3]), "create_omniscient");

# run get_feature_l2_from_id_l2_l1
my $l1_ID = lc($feature1->_tag_value("ID"));
ok( get_feature_l2_from_id_l2_l1($hash_omniscient, $l2_ID, $l1_ID), "get_feature_l2_from_id_l2_l1");

# run location_overlap_update
ok(  location_overlap_update(["locA",100,200],["locB",150,300]), "location_overlap_update");

# run check_record_positions
ok( check_record_positions($hash_omniscient, $l1_ID), "check_record_positions");

# run l1_has_l3_type
ok( l1_has_l3_type ($hash_omniscient, $feature1, "cds", undef), "l1_has_l3_type");

# run l2_has_cds
ok( l2_has_cds ($hash_omniscient, $feature2), "l2_has_cds");

# run get_cds_from_l2
ok( get_cds_from_l2 ($hash_omniscient, $feature2), "get_cds_from_l2");

# run l2_has_cds
ok( l2_has_cds ($hash_omniscient, $feature2), "l2_has_cds");

# run check_if_feature_overlap
ok( check_features_overlap($feature1, $feature2), "check_features_overlap");

# run is_single_exon_gene
$feature1 = $hash_omniscient->{"level1"}{"gene"}{"gene:os01g0100650"};
ok( is_single_exon_gene($hash_omniscient, $feature1), "is_single_exon_gene");
