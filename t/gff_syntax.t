#!/usr/bin/env perl

use strict;
use warnings;
use Test::More tests => 35;

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
my $pathtmp = "tmp.gff"; # path file where to save temporary output
#my $pathtmp2 = "tmp2.gff"; # path file where to save temporary output
my $dir = "t/gff_syntax"; # folder where the test files are
opendir my $dh, $dir or die "Could not open '$dir' for reading: $!\n";
my @files = readdir $dh;
foreach my $file (sort { (($a =~ /^(\d+)/)[0] || 0) <=> (($b =~ /^(\d+)/)[0] || 0) } @files) {

  # for all test files
  if ( $file =~ m/test.gff$/ ){

    # skip cases
    if ( ($file =~ m/^24_*/) or  ($file =~ m/^26_*/) ){ #Skip those cases
        next;
    }

    # case do not merge loci 8,32,34,36
    if ($file =~ m/^8_/ or $file =~ m/^33_/ or $file =~ m/^34_/ or $file =~ m/^36_/){
        system("$script --gff t/gff_syntax/$file -o $pathtmp  2>&1 1>/dev/null");
    }
		# peculiar case 28
    elsif($file =~ m/^28_/){
        system("$script --gff t/gff_syntax/$file -c Name -o $pathtmp  2>&1 1>/dev/null");
    }

    # standard cases
    else{
      system("$script --gff t/gff_syntax/$file  --merge_loci -o $pathtmp  2>&1 1>/dev/null");
    }

    my @splitname = split /_/, $file;
    my $correct_output = $splitname[0]."_correct_output.gff";

    #run test
    ok( system("diff $pathtmp t/gff_syntax/$correct_output") == 0, "parse $file");

    #--------------- rerun on the first output #---------------

     # peculiar case
     #if ($file =~ m/^8_/){
    #     system("$script --gff $pathtmp -o $pathtmp2 ");
     #}
     #elsif($file =~ m/^28_/){
    #     system("$script --gff $pathtmp -c Name -o $pathtmp2 ");
     #}
     # standard cases
     #else{
    #   system("$script --gff $pathtmp --merge_loci -o $pathtmp2 ");
     #}

     #run test
     #ok( system("diff $pathtmp2 t/gff_syntax/$correct_output") == 0, "parse2 $file");
  }
}
closedir($dh);


unlink $pathtmp;
#unlink $pathtmp2;
