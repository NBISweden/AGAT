#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use File::Path qw(remove_tree); # to remove directory easily (tmp directory)
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
my $script_agat = $script_prefix."bin/agat";
my $script = $script_prefix."bin/agat_convert_sp_gxf2gxf.pl";
my $input_folder = "t/level_missing/in";
my $output_folder = "t/level_missing/out";
my $outprefix = "tmp";
my $outtmp = "$outprefix.gff"; # path file where to save temporary output
my $result;
my $config="agat_config.yaml";

# remove config in local folder if exists
cleaning_log();

# -------------------------- testA -------------------------

$result = "$output_folder/testA_output.gff";
system(" $script --gff $input_folder/testA.gff -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output testA");
cleaning_log("testA.gff");

$result = "$output_folder/testA_output2.gff";
system("$script_agat config --expose --locus_tag common_tag 2>&1 1>/dev/null");
system(" $script --gff $input_folder/testA.gff -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output testA2");
cleaning_log("testA.gff");

$result = "$output_folder/testA_output3.gff";
system("$script_agat config --expose --locus_tag gene_info 2>&1 1>/dev/null");
system(" $script --gff $input_folder/testA.gff -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output testA3");
cleaning_log("testA.gff");

$result = "$output_folder/testA_output4.gff";
system("$script_agat config --expose --locus_tag transcript_id 2>&1 1>/dev/null");
system(" $script --gff $input_folder/testA.gff -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output testA4");
cleaning_log("testA.gff");

# -------------------------- testB -------------------------

$result = "$output_folder/testB_output.gff";
system(" $script --gff $input_folder/testB.gff -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output testB");
cleaning_log("testB.gff");

$result = "$output_folder/testB_output2.gff";
system("$script_agat config --expose --locus_tag locus_id 2>&1 1>/dev/null");
system(" $script --gff $input_folder/testB.gff -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output testB2");
cleaning_log("testB.gff");

# -------------------------- testC -------------------------

$result = "$output_folder/testC_output.gff";
system(" $script --gff $input_folder/testC.gff -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output testC");
cleaning_log("testC.gff");

$result = "$output_folder/testC_output2.gff";
system("$script_agat config --expose --locus_tag locus_id 2>&1 1>/dev/null");
system(" $script --gff $input_folder/testC.gff -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output testC2");
cleaning_log("testC.gff");

# -------------------------- testD -------------------------

$result = "$output_folder/testD_output.gff";
system(" $script --gff $input_folder/testD.gff -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output testD");
cleaning_log("testD.gff");

$result = "$output_folder/testD_output2.gff";
system("$script_agat config --expose --locus_tag ID 2>&1 1>/dev/null");
system(" $script --gff $input_folder/testD.gff -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output testD2");
cleaning_log("testD.gff");

# -------------------------- testE -------------------------

$result = "$output_folder/testE_output.gff";
system(" $script --gff $input_folder/testE.gff -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output testE");
cleaning_log("testE.gff");


# -------------------------- testF -------------------------

$result = "$output_folder/testF_output.gff";
system(" $script --gff $input_folder/testF.gff -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output testF");
cleaning_log("testF.gff");

# -------------------------- testG -------------------------

$result = "$output_folder/testG_output.gff";
system(" $script --gff $input_folder/testG.gff -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output testG");
cleaning_log("testG.gff");

# --- convenient function ---

sub cleaning_log{
  my ($filename, $local_config)=@_;

  # REMOVE LOG folder if a file name is provided
  if($filename){
    my ($name, $path, $suffix) = fileparse($filename, qr/\.[^.]*/);
    if (-e "agat_log_$name"){
      remove_tree( "agat_log_$name" );
    }
  }

  # remove config
  # Si la variable $local_config est définie et vraie (c’est-à-dire qu’elle n’est pas undef, ni vide, ni 0), alors $config_to_remove prendra sa valeur. Sinon, il prendra la valeur de $config.
  my $config_to_remove = $local_config ? $local_config : $config;
  if (-e $config_to_remove){
    unlink $config_to_remove;
  }

  # the rest
  opendir(my $dh, ".") or die "Cannot open current directory: $!";
  while (my $file = readdir($dh)) {
    next if $file =~ /^\./;  # Skip hidden files and . / ..
    if ($file =~ /^\Q$outprefix\E/ or $file =~ /^\Q$outtmp\E/) {
      if ( -d $file ) {
        remove_tree($file);
        #print "remove_tree $file\n";
      } else {
        unlink $file;
        #print "unlink $file\n";
      }
    }
  }
  closedir($dh);
}