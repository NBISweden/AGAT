#!/usr/bin/env perl
use strict;
use warnings;
use Test::More;
use File::Temp qw(tempdir);
use Cwd qw(getcwd);
use AGAT::AGAT;
use AGAT::Config;

# Work in a clean temporary directory without a local config file
my $orig = getcwd();
my $dir  = tempdir(CLEANUP => 1);
chdir $dir;

# Simulate a quiet run by setting global verbosity
$AGAT::AGAT::CONFIG = { verbose => 0 };

my $output = '';
{
    open my $fh, '>', \$output;
    local *STDOUT = $fh;
    get_config({ type => 'original' });
}

is( $output, '', 'get_config respects global verbosity when none provided' );

chdir $orig;

done_testing();
