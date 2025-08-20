#!/usr/bin/env perl
use strict;
use warnings;
use Test::More;
use Getopt::Long::Descriptive;
use AGAT::CLI::Common qw(common_spec);

my @args = (
    '--config', 'foo.yaml', '--output', 'bar',
    '--log', 'baz.log', '--verbose', '2', '--debug', '--quiet',
    '--extra', 'val'
);
local @ARGV = @args;
my ($opt) = describe_options('test %o', common_spec());

is_deeply(
    { %{$opt} },
    { config => 'foo.yaml', out => 'bar', log => 'baz.log', verbose => 2, debug => 1, quiet => 1 },
    'parsed values'
);
is_deeply( \@ARGV, ['--extra', 'val'], 'remaining args preserved' );

done_testing;
