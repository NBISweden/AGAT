#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use File::Path qw(remove_tree); # to remove directory easily (tmp directory)
use Test::More tests => 8;

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

# script to call to check the parser
my $script = "";
my $script_agat = $script_prefix."bin/agat";
my $input_folder = "t/gff_other/in";
my $output_folder = "t/gff_other/out";
my $outprefix = "tmp";
my $outtmp = "$outprefix.gff"; # path file where to save temporary output
my $config="agat_config.yaml";

# remove config in local folder if exists and potential tmp file already existing
cleaning_log();

# -------- test gzip file and contain fasta --------
$script = $script_prefix."bin/agat_convert_sp_gxf2gxf.pl";
my $correct_output = "$output_folder/zip_and_fasta_correct_output.gff";

system("$script --gff $input_folder/zip_and_fasta.gff.gz -o $outtmp  2>&1 1>/dev/null");

#run test
ok( system("diff $outtmp $correct_output") == 0, "zip_and_fasta check");
cleaning_log("zip_and_fasta.gff.gz");

# -------- test tabix output sorting --------
$script = $script_prefix."bin/agat_convert_sp_gxf2gxf.pl";
$correct_output = "$output_folder/1_agat_tabix.gff";

system("$script_agat config --expose --tabix 2>&1 1>/dev/null");
system("$script --gff t/scripts_output/in/1.gff -o $outtmp 2>&1 1>/dev/null");

#run test
ok( system("diff $outtmp $correct_output") == 0, "tabix check");
cleaning_log("1.gff");

# -------- Parent ID already used by same level feature --------
$script = $script_prefix."bin/agat_convert_sp_gxf2gxf.pl";
$correct_output = "$output_folder/issue329.gff";

system("$script --gff $input_folder/issue329.gff -o $outtmp 2>&1 1>/dev/null");

#run test
ok( system("diff $outtmp $correct_output") == 0, "issue329 check");
cleaning_log("issue329.gff");

# -------- Issue 368 avoid seq_id 0 replaced by SEQ --------
$script = $script_prefix."bin/agat_convert_sp_gxf2gxf.pl";
$correct_output = "$output_folder/issue368.gff";

system("$script --gff $input_folder/issue368.gff -o $outtmp 2>&1 1>/dev/null");

#run test
ok( system("diff $outtmp $correct_output") == 0, "issue368 check");
cleaning_log("issue368.gff");

# -------- Issue 389 very wierd --------
$script = $script_prefix."bin/agat_convert_sp_gxf2gxf.pl";
$correct_output = "$output_folder/issue389.gff";

system("$script --gff $input_folder/issue389.gff -o $outtmp 2>&1 1>/dev/null");

#run test
ok( system("diff $outtmp $correct_output") == 0, "issue389 check");
cleaning_log("issue389.gff");

# --------- Issue 441 transccipt_id used for GTF while L2 L1 created from L3 with isoforms ----
$script = $script_prefix."bin/agat_convert_sp_gxf2gxf.pl";
$correct_output = "$output_folder/issue441.gtf";

system("$script_agat config --expose --output_format gtf 2>&1 1>/dev/null");
system("$script --g $input_folder/issue441.gtf -o $outtmp  2>&1 1>/dev/null");

ok( system("diff $outtmp $correct_output") == 0, "issue441 check");
cleaning_log("issue441.gff");

# --------- Issue 448 bioperl adding extra empty ID attribute that mess up AGAT (only when input parsed with version 2 and 2.5)  ----
$script = $script_prefix."bin/agat_convert_sp_gxf2gxf.pl";
$correct_output = "$output_folder/issue448.gtf";

system("$script_agat config --expose --output_format gtf 2>&1 1>/dev/null");
system("$script --g $input_folder/issue448.gtf -o $outtmp  2>&1 1>/dev/null");

ok( system("diff $outtmp $correct_output") == 0, "issue441 check");
cleaning_log("issue448.gff");

$correct_output = "$output_folder/issue448.gff";
system("$script --g $input_folder/issue448.gtf -o $outtmp  2>&1 1>/dev/null");

ok( system("diff $outtmp $correct_output") == 0, "issue441 check");
cleaning_log("issue448.gff");

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
  if (-e $config){
    unlink $config;
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