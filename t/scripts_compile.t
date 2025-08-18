#!/usr/bin/env perl

use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/lib";
use File::Spec::Functions qw(catdir catfile);
use Cwd qw(abs_path);
use AGAT::TestUtilities qw(setup_tempdir check_diff);

my $root    = abs_path(catdir($Bin, '..'));
my $BIN_DIR = catdir($root, 'bin');

=head1 DESCRIPTION

Test to see if all script in the bin file can be compiled without error.

=cut

use Test::More;

# foreach script in the bin, let run the test
opendir (DIR, $BIN_DIR) or die $!;
while (my $file = readdir(DIR)) {
   # Use a regular expression to ignore files beginning with a period
   next if ($file =~ m/^\./);
   my $path = catfile($BIN_DIR, $file);
   {
       my $dir = setup_tempdir();
       print "$path -h 1>/dev/null\n";
       #run test - check the script can run calling the help.
       ok( system("$path -h 1>/dev/null") == 0, "test $file" );
   }
}
closedir(DIR);

done_testing();
