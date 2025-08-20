#!/usr/bin/env perl
use strict;
use warnings;
use FindBin;
our ($FILE, $EXPECTED);
$FILE = 'stop_split_over_two_exons.gtf';
$EXPECTED = 'stop_split_over_two_exons_correct_output.gff';
require "$FindBin::Bin/case_base.pl";
