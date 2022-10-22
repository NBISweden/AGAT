#!/usr/bin/env perl

use strict;
use warnings;

=head1 DESCRIPTION

Test to see if all script in the bin file can be compiled without error.

=cut

################################################################################
#   set number of test according to number of scripts
my $nb_test;
BEGIN{
  opendir (DIR, "bin") or die $!;
  while (my $file = readdir(DIR)) {
     # Use a regular expression to ignore files beginning with a period
     next if ($file =~ m/^\./);

     #add exe file
     $nb_test++;
  }
  closedir(DIR);
}
#
################################################################################

use Test::More tests => $nb_test ;

# remove config in local folder if exists
unlink "config.yaml"; 

# foreach script in the bin, let run the test
opendir (DIR, "bin") or die $!;
while (my $file = readdir(DIR)) {
   # Use a regular expression to ignore files beginning with a period
   next if ($file =~ m/^\./);
   print "bin/$file -h 1>/dev/null\n";
   #run test - check the script can run calling the help.
   ok( system("bin/$file -h 1>/dev/null") == 0, "test $file")
}
closedir(DIR);
