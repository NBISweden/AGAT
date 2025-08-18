#!/usr/bin/env perl

use strict;
use warnings;
use AGAT::AGAT;
use Test::More tests => 6; # half of file to test but each tested twice
use FindBin qw($Bin);
use lib "$Bin/lib";
use File::Spec::Functions qw(catdir catfile);
use AGAT::TestUtilities qw(setup_tempdir check_diff);

=head1 DESCRIPTION

Test to verify the method detecting the gff parser version to use with bioperl.

=cut

# remove config in local folder if exists
my $config="agat_config.yaml";

# Loop over test
my $input_dir = catdir($Bin, 'gff_version', 'in');
opendir (DIR, $input_dir) or die $!;
while (my $file = readdir(DIR)) {

  # for all test files
  if ( $file =~ m/test.gff$/ ){
    {
        my $dir = setup_tempdir();
        # format detected by select_gff_format
        my $format = select_gff_format(catfile($input_dir, $file));
        # first character is the format of the file expected
        my $firstchar = substr $file, 0, 1;

        ok( $firstchar == $format, "Good version found for $file" );
    }
  }
}
closedir(DIR);
