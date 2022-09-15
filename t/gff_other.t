#!/usr/bin/env perl

use strict;
use warnings;
use Test::More tests => 2;
use Bio::Tools::GFF;
use AGAT::Omniscient;
use AGAT::OmniscientTool;

=head1 DESCRIPTION

Test to verify other GFF/GTF peculiarities e.g HEADER / FASTA / GZIP

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

# script to call to check the parser
my $script = "";
my $input_folder = "t/gff_other/in";
my $output_folder = "t/gff_other/out";
my $pathtmp = "tmp.gff"; # path file where to save temporary output

# -------- test gzip file and contain fasta --------
$script = $script_prefix."bin/agat_convert_sp_gxf2gxf.pl";
my $correct_output = "$output_folder/zip_and_fasta_correct_output.gff";

unlink "config.yaml"; # remove config in local folder if exists
system("$script --gff $input_folder/zip_and_fasta.gff.gz -o $pathtmp  2>&1 1>/dev/null");

#run test
ok( system("diff $pathtmp $correct_output") == 0, "zip_and_fasta check");
unlink $pathtmp;

# -------- test tabix output sorting --------
$script = $script_prefix."bin/agat_convert_sp_gxf2gxf.pl";
$correct_output = "$output_folder/1_agat_tabix.gff";

system("agat config --expose --tabix 2>&1 1>/dev/null");
system("$script --gff t/scripts_output/1.gff -o $pathtmp  2>&1 1>/dev/null");

#run test
ok( system("diff $pathtmp $correct_output") == 0, "tabix check");
unlink $pathtmp;
unlink "config.yaml";
