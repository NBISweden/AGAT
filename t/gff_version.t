#!/usr/bin/env perl

use strict;
use warnings;
use AGAT::AGAT;
use Test::More tests => 6; # half of file to test but each tested twice

=head1 DESCRIPTION

Test to verify the method detecting the gff parser version to use with bioperl.

=cut

# remove config in local folder if exists
unlink "config.yaml";

# Loop over test
opendir (DIR, "t/gff_version/in") or die $!;
while (my $file = readdir(DIR)) {

  # for all test files
  if ( $file =~ m/test.gff$/ ){

    # format detected by select_gff_format
    my $format = select_gff_format("t/gff_version/in/$file");
    # first character is the format of the file expected
    my $firstchar = substr $file, 0, 1;

    ok( $firstchar == $format, "Good version found for $file");
  }
}
closedir(DIR);
