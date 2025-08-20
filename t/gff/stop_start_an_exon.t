#!/usr/bin/env perl
use strict;
use warnings;
use FindBin;
our ($FILE, $EXPECTED);
$FILE = 'stop_start_an_exon.gtf';
$EXPECTED = 'stop_start_an_exon_correct_output.gff';
require "$FindBin::Bin/case_base.pl";
