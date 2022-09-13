#!/usr/bin/env perl

use strict;
use warnings;
use Test::More tests => 13;

=head1 DESCRIPTION

Test to verify the parser deals properly with the different flavor / bugged gff files.

=cut

# Check if has to be run in Devel::Cover or not
my $script_prefix="";
if (exists $ENV{'HARNESS_PERL_SWITCHES'} ) {
  if ($ENV{'HARNESS_PERL_SWITCHES'} =~ m/Devel::Cover/) {
    $script_prefix="perl -MDevel::Cover ";
  }
}

# script to call to check the parser
my $script = $script_prefix."bin/agat_convert_sp_gxf2gxf.pl";
my $input_folder = "t/level_missing";
my $output_folder = "t/level_missing/output";
my $outtmp = "tmp.gff"; # path file where to save temporary output
my $result;
unlink "config.yaml"; # remove custom config if exists

# -------------------------- testA -------------------------

$result = "$output_folder/testA_output.gff";
system(" $script --gff $input_folder/testA.gff -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
unlink $outtmp;

$result = "$output_folder/testA_output2.gff";
system("agat config --expose --locus_tag common_tag 2>&1 1>/dev/null");
system(" $script --gff $input_folder/testA.gff -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
unlink $outtmp;
unlink "config.yaml";

$result = "$output_folder/testA_output3.gff";
system("agat config --expose --locus_tag gene_info 2>&1 1>/dev/null");
system(" $script --gff $input_folder/testA.gff -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
unlink $outtmp;
unlink "config.yaml";

$result = "$output_folder/testA_output4.gff";
system("agat config --expose --locus_tag transcript_id 2>&1 1>/dev/null");
system(" $script --gff $input_folder/testA.gff -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
unlink $outtmp;
unlink "config.yaml";

# -------------------------- testB -------------------------

$result = "$output_folder/testB_output.gff";
system(" $script --gff $input_folder/testB.gff -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
unlink $outtmp;

$result = "$output_folder/testB_output2.gff";
system("agat config --expose --locus_tag locus_id 2>&1 1>/dev/null");
system(" $script --gff $input_folder/testB.gff -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
unlink $outtmp;
unlink "config.yaml";

# -------------------------- testC -------------------------

$result = "$output_folder/testC_output.gff";
system(" $script --gff $input_folder/testC.gff -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
unlink $outtmp;

$result = "$output_folder/testC_output2.gff";
system("agat config --expose --locus_tag locus_id 2>&1 1>/dev/null");
system(" $script --gff $input_folder/testC.gff -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
unlink $outtmp;
unlink "config.yaml";

# -------------------------- testD -------------------------

$result = "$output_folder/testD_output.gff";
system(" $script --gff $input_folder/testD.gff -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
unlink $outtmp;

$result = "$output_folder/testD_output2.gff";
system("agat config --expose --locus_tag ID 2>&1 1>/dev/null");
system(" $script --gff $input_folder/testD.gff -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
unlink $outtmp;
unlink "config.yaml";

# -------------------------- testE -------------------------

$result = "$output_folder/testE_output.gff";
system(" $script --gff $input_folder/testE.gff -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
unlink $outtmp;


# -------------------------- testF -------------------------

$result = "$output_folder/testF_output.gff";
system(" $script --gff $input_folder/testF.gff -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
unlink $outtmp;

# -------------------------- testG -------------------------

$result = "$output_folder/testG_output.gff";
system(" $script --gff $input_folder/testG.gff -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output $script");
unlink $outtmp;
