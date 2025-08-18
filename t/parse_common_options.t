#!/usr/bin/env perl
use strict;
use warnings;
use Test::More tests => 2;
use AGAT::AGAT;

{
    local @ARGV = ('--config','foo.yaml','--output','bar','--verbose','--extra','val');
    my $opts = parse_common_options();
    is_deeply($opts, {config=>'foo.yaml', output=>'bar', verbose=>1}, 'parsed values');
    is_deeply(\@ARGV, ['--extra','val'], 'remaining args preserved');
}
