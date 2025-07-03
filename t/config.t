#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use File::Path qw(remove_tree); # to remove directory easily (tmp directory)
use Test::More tests => 3;

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
my $output_folder = "t/config/out";
my $config="agat_config.yaml";
my $outprefix = "tmp";
my $outtmp = "$outprefix.gff"; # path file where to save temporary output

# remove config in local folder if exists
cleaning_log();

# ---- test gzip file and contain fasta ----
$script = $script_prefix."bin/agat";
my $correct_output = "$output_folder/$config";

system("$script config -e \\
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
								--clean_attributes_from_template  ");

#run test
ok( system("diff $config $correct_output") == 0, "modif config check");
# remove file created for the test
unlink $config;

# ----- Test Rename config file ----

my $new_config_name = "agat_config_renamed.yml";
system("$script config -e --output $new_config_name --locus_tag Name");

#run test 
ok( system("if [ -e $new_config_name ];then exit 0;fi") == 0, "rename agat config file check");

# ----- Test use a renamed config file ----

system("bin/agat_convert_sp_gxf2gxf.pl --gff t/gff_syntax/in/28_test.gff -c $new_config_name -o $outtmp  2>&1 1>/dev/null");
#run test 
ok( system("diff tmp.gff t/gff_syntax/out/28_correct_output.gff") == 0, "Use custom agat config file check");

# remove file created for the test
cleaning_log("28_test.gff", $new_config_name);

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