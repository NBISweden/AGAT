#!/usr/bin/env perl
use strict;
use warnings;
use Test::More;
use Getopt::Long::Descriptive;
use AGAT::AGAT qw(common_spec);

my @args = (
    '--config', 'foo.yaml', '--output', 'bar',
    '--log', 'baz.log', '--quiet',
    '--extra', 'val'
);
local @ARGV = @args;
my ($opt) = describe_options('test %o', common_spec());

is_deeply(
    { %{$opt} },
    {
        config       => 'foo.yaml',
        out          => 'bar',
        log          => 'baz.log',
        verbose      => 0,
        debug        => 0,
        progress_bar => 0,
        quiet        => 1,
    },
    'parsed values',
);
is_deeply( \@ARGV, ['--extra', 'val'], 'remaining args preserved' );

done_testing;
