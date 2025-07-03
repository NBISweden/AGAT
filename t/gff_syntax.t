#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use File::Path qw(remove_tree); # to remove directory easily (tmp directory)
use Test::More tests => 47;

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
my $outprefix = "tmp";
my $outtmp = "$outprefix.gff"; # path file where to save temporary output
my $expected_output_path = "t/gff_syntax/out/";
my $input_path = "t/gff_syntax/in/";
my $config="agat_config.yaml";

# remove config in local folder if exists
cleaning_log();

# Loop over test
my $dir = "t/gff_syntax/in"; # folder where the test files are
opendir my $dh, $dir or die "Could not open '$dir' for reading: $!\n";
my @files = readdir $dh;
foreach my $file (sort { (($a =~ /^(\d+)/)[0] || 0) <=> (($b =~ /^(\d+)/)[0] || 0) } @files) {

  # skip few cases to be faster
  #my ($value) = $file =~ m/^(\d+)_.*/;
  #my $analyzed = 44;
  #if(defined($value) and $value < $analyzed){next;}
  

  # for all test files
  if ( $file =~ m/test.gff$/ ){

    # skip cases
    if ( ($file =~ m/^24_*/) or  ($file =~ m/^26_*/) ){ #Skip those cases
        next;
    }

    # case do not merge loci 8,32,34,36
    if ($file =~ m/^8_/ or $file =~ m/^33_/ or $file =~ m/^34_/ or $file =~ m/^36_/){
        system("$script --gff $input_path/$file -o $outtmp  2>&1 1>/dev/null");
    }
		# peculiar cases with locus_tag Name
    elsif($file =~ m/^28_/ or $file =~ m/^45_/ or $file =~ m/^46_/){
        system("$script_agat config --expose --locus_tag Name 2>&1 1>/dev/null"); # set special config for the test
        system("$script --gff $input_path/$file -o $outtmp  2>&1 1>/dev/null");
    }
    # standard cases
    else{
			system("$script_agat config --expose --merge_loci 2>&1 1>/dev/null"); # set special config for the test
      system("$script --gff $input_path/$file -o $outtmp  2>&1 1>/dev/null");
    }

    my @splitname = split /_/, $file;
    my $correct_output = $expected_output_path."/".$splitname[0]."_correct_output.gff";

    #run test
    ok( system("diff $outtmp $correct_output") == 0, "parse $file");

    # cleaning 
    #if($value and $value == $analyzed){exit;}
    cleaning_log($file);
  }
}
closedir($dh);

#  ---------------------------- special tests ----------------------------
my $result = "$expected_output_path/stop_start_an_exon_correct_output.gff";
system(" $script --gff $input_path/stop_start_an_exon.gtf -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output stop_start_an_exon");
cleaning_log("stop_start_an_exon.gtf");

$result = "$expected_output_path/stop_split_over_two_exons_correct_output.gff";
system(" $script --gff $input_path/stop_split_over_two_exons.gtf -o $outtmp 2>&1 1>/dev/null");
#run test
ok( system("diff $result $outtmp") == 0, "output stop_split_over_two_exons_correct_output");
cleaning_log("stop_split_over_two_exons.gtf");

# --- convenient function ---

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
