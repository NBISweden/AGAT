#!/usr/bin/env perl
use strict;
use warnings;
use Test::More;
use File::Temp qw(tempdir);
use File::Copy qw(copy);
use File::Spec;
use AGAT::AGAT;

# Use a temporary config file to test layering
my $tmpdir = tempdir(CLEANUP => 1);
my $config_path = File::Spec->catfile($tmpdir, 'agat_config.yaml');
copy('share/agat_config.yaml', $config_path);

# Ensure config file sets verbose to 2
open my $fh, '<', $config_path or die $!;
my @lines = <$fh>;
close $fh;
for (@lines) { s/^verbose: .*/verbose: 2/; }
open $fh, '>', $config_path or die $!;
print {$fh} @lines;
close $fh;

{
    local @ARGV = ('--config', $config_path);
    my $cfg = resolve_common_options();
    is($cfg->{verbose}, 2, 'config file value used');
}

{
    local @ARGV = ('--config', $config_path, '--verbose', '0');
    my $cfg = resolve_common_options();
    is($cfg->{verbose}, 0, 'command line overrides config');
}

done_testing;
