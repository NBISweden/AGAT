#!/usr/bin/perl -w

package AGAT::SeqFeatureLite;
use parent 'Bio::SeqFeature::Lite';

use strict;
use warnings;
use Exporter;

our @ISA = qw(Exporter);
our @EXPORT = qw(_tag_value);

sub _tag_value {
    my ($self, $tag) = @_;
    return unless $self->has_tag($tag);
    my ($val) = $self->get_tag_values($tag);
    return $val;
}

1;