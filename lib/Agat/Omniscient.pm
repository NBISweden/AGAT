#!/usr/bin/perl -w

package Agat::Omniscient;

use strict;
use warnings;
use Exporter;

use Agat::OmniscientI;
use Agat::OmniscientO;
use Agat::OmniscientTool;
use Agat::OmniscientStat;
use Agat::PlotR;

our $VERSION     = 0.01;
our @ISA         = qw(Exporter);
our @EXPORT      = qw(get_agat_header);
sub import {
  Agat::Omniscient->export_to_level(1, @_); # to be able to load the EXPORT functions when direct call; (normal case)
}

=head1 SYNOPSIS

  Meta package for conveniency. It allows to call all packages needed in once to deal with Omniscient data structure.

=head1 DESCRIPTION

    Omniscient packages are non-OO packages use to handle any kind of gtf/gff data.

=head1 CONTACT

    jacques.dainat@nbis.se (Jacques Dainat)

=cut

# Provide meta information
sub get_agat_header{

  my ($verbose) = @_;

  my $header = qq{
 ------------------------------------------------------------------------------
|   Another GFF Analysis Toolkit (AGAT) - Version: $VERSION                        |
|   https://github.com/NBISweden/AGAT                                          |
|   National Bioinformatics Infrastructure Sweden (NBIS) - www.nbis.se         |
 ------------------------------------------------------------------------------
  };

return $header;

}
1;
