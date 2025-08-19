#!/usr/bin/env perl
use strict;
use warnings;
use Test::More tests => 3;
use AGAT::AGAT;

{
    my @args = (
        '--config', 'foo.yaml', '--output', 'bar',
        '--verbose', '--debug', '--help', '--extra', 'val'
    );
    local @ARGV = @args;
    my $opts = parse_common_options();
    my $orig = delete $opts->{argv};
    is_deeply(
        $opts,
        { config => 'foo.yaml', output => 'bar', verbose => 1, debug => 1, help => 1 },
        'parsed values'
    );
    is_deeply( \@ARGV, [ '--extra', 'val' ], 'remaining args preserved' );
    is_deeply( $orig, \@args, 'original argv captured' );
}

