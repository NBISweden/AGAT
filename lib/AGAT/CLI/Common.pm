package AGAT::CLI::Common;

use strict;
use warnings;
use Exporter 'import';
use AGAT::AGAT qw(resolve_common_options);

our @EXPORT_OK = qw(common_spec resolve_config);

sub common_spec {
    return (
        [ 'config|c=s',          'Configuration file' ],
        [ 'out|o|output=s',      'Output GFF3 file' ],
        [ 'log=s',               'Log file path' ],
        [ 'verbose|v=i',         'Verbosity level' ],
        [ 'debug|d',             'Enable debug output' ],
        [ 'progress_bar|progressbar!', 'Show progress bar', { default => undef } ],
        [ 'quiet|q',             'Disable progress bar and verbose output',
            { implies => { debug => 0, verbose => 0, progress_bar => 0 } } ],
        [ 'help|h',              'Show this help', { shortcircuit => 1 } ],
        { getopt_conf => ['pass_through'] },
    );
}

sub resolve_config {
    my ($opt) = @_;
    my %cli = %{ $opt || {} };
    $cli{output} = delete $cli{out} if exists $cli{out};
    return AGAT::AGAT::resolve_common_options( \%cli );
}

1;

__END__

=head1 NAME

AGAT::CLI::Common - Shared command line option specifications for AGAT scripts

=head1 SYNOPSIS

  use AGAT::CLI::Common qw(common_spec resolve_config);
  my ($opt,$usage) = describe_options('%c %o',
                                     [ 'input=s', 'Input file', { required => 1 } ],
                                     common_spec());
  my $config = resolve_config($opt);

=head1 DESCRIPTION

Provides reusable option descriptors and helpers to merge command line
options with configuration defaults.

=cut
