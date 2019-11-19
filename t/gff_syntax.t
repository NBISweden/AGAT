#!/usr/bin/env perl

use strict;
use warnings;

use Test::Simple tests => 62; # half of file to test but each tested twice

=head1 DESCRIPTION

Test to verify the parser deals properly with the different flavor / bugged gff files.

=cut

# script to call to check the parser
my $handler_script = "bin/agat_sp_gxf_to_gff3.pl";
my $pathtmp = "tmp.gff"; # path file where to save temporary output
my $pathtmp2 = "tmp2.gff"; # path file where to save temporary output
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

    # peculiar case
    if ($file =~ m/^8_/){
        system("$handler_script --gff t/gff_syntax/$file -o $pathtmp ");
    }
    elsif($file =~ m/^28_/){
        system("$handler_script --gff t/gff_syntax/$file -c Name -o $pathtmp ");
    }
    # standard cases
    else{
      system("$handler_script --gff t/gff_syntax/$file  --merge_loci -o $pathtmp ");
    }

    my @splitname = split /_/, $file;
    my $correct_output = $splitname[0]."_correct_output.gff";

    #run test
    ok( system("diff $pathtmp t/gff_syntax/$correct_output") == 0, "parse $file");

    #--------------- rerun on the first output #---------------

     # peculiar case
     if ($file =~ m/^8_/){
         system("$handler_script --gff $pathtmp -o $pathtmp2 ");
     }
     elsif($file =~ m/^28_/){
         system("$handler_script --gff $pathtmp -c Name -o $pathtmp2 ");
     }
     # standard cases
     else{
       system("$handler_script --gff $pathtmp --merge_loci -o $pathtmp2 ");
     }

     #run test
     ok( system("diff $pathtmp2 t/gff_syntax/$correct_output") == 0, "parse2 $file");
  }
}
closedir($dh);


unlink $pathtmp;
unlink $pathtmp2;
