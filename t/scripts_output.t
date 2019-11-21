#!/usr/bin/env perl

use strict;
use warnings;

use Test::Simple tests => 2;

=head1 DESCRIPTION

Test to verify the output of scripts using gff_syntax cases

=cut

# shared variables
my $output_folder = "t/scripts_output";
my $outtmp = "tmp.gff"; # path file where to save temporary output

# ------------------- check agat_sp_merge_annotations script-------------------
my $script = "bin/agat_sp_merge_annotations.pl";
my $result = "$output_folder/agat_sp_merge_annotations_1.gff";
system("$script --gff t/gff_syntax/25_test.gff  --gff t/gff_syntax/9_test.gff -o $outtmp ");
#run test
ok( system("diff $outtmp $outtmp") == 0, "output $script");



# -------------------------- check agat_sp_statistics --------------------------

my $script = "bin/agat_sp_statistics.pl";
my $result = "$output_folder/agat_sp_statistics_1.gff";
system("$script --gff t/gff_syntax/0_test.gff -o $outtmp ");
#run test
ok( system("diff $outtmp $outtmp") == 0, "output $script");


# ---------------------------- remove $outtmp file -----------------------------
unlink $outtmp;
