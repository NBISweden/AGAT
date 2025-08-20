#!/usr/bin/env perl
use strict;
use warnings;
use Test::More;
use File::Temp qw(tempdir);
use File::Copy qw(copy);
use File::Spec;
use Getopt::Long::Descriptive;
use AGAT::CLI::Common qw(common_spec resolve_config);

# prepare temporary config file
my $tmpdir      = tempdir(CLEANUP => 1);
my $config_path = File::Spec->catfile($tmpdir, 'agat_config.yaml');
copy('share/agat_config.yaml', $config_path);

# 1. config file value used
{
    local @ARGV = ('--config', $config_path);
    my ($opt) = describe_options('test %o', common_spec());
    my $cfg = resolve_config($opt);
    is($cfg->{verbose}, 1, 'config file value used');
}

# 2. CLI overrides config
{
    local @ARGV = ('--config', $config_path, '--verbose', '2');
    my ($opt) = describe_options('test %o', common_spec());
    my $cfg = resolve_config($opt);
    is($cfg->{verbose}, 2, 'command line overrides config');
}

# 3. CLI log path mapped to log_path
{
    local @ARGV = ('--log', 'foo.log');
    my ($opt) = describe_options('test %o', common_spec());
    my $cfg = resolve_config($opt);
    is($cfg->{log_path}, 'foo.log', 'CLI log path mapped to log_path');
}

# 4. default log path when none provided
{
    local @ARGV = ();
    my ($opt) = describe_options('test %o', common_spec());
    my $cfg = resolve_config($opt);
    my ($file) = $0 =~ /([^\\\/]+)$/;
    is($cfg->{log_path}, "$file.agat.log", 'default log path derived from script name');
}

# 5. quiet profile overrides config
{
    local @ARGV = ('--config', $config_path, '--quiet');
    my ($opt) = describe_options('test %o', common_spec());
    my $cfg = resolve_config($opt);
    is($cfg->{verbose}, 0, 'quiet sets verbose to 0');
    is($cfg->{progress_bar}, 0, 'quiet disables progress bar');
}

done_testing;
