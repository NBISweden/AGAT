#!/usr/bin/perl -w

package AGAT::SeqFeatureLite;
use parent 'Bio::SeqFeature::Lite';

use strict;
use warnings;

# convenience method to get single tag value
sub _tag_value {
    my ($self, $tag) = @_;
    return unless $self->has_tag($tag);
    my ($val) = $self->get_tag_values($tag);
    return $val;
}

# frame methods call frame mehtods of parent class
# no need to override them here
sub frame {
    return shift->phase;
}

1;