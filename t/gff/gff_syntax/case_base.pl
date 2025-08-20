#!/usr/bin/env perl

use strict;
use warnings;
use Test::More tests => 1;
use File::Basename qw(dirname);
use File::Spec::Functions qw(catdir catfile);
use Cwd qw(abs_path);

BEGIN {
    our $CASE_DIR = dirname(abs_path(__FILE__));
    unshift @INC, "$CASE_DIR/../../lib";
}

use AGAT::TestUtilities qw(setup_tempdir check_diff script_prefix);

our ($FILE, $EXPECTED);  # input and expected filenames relative to t/gff/gff_syntax
our $CASE_DIR;           # populated in BEGIN

my $file = $FILE // die 'No $FILE provided';
my $expected_override = $EXPECTED; # may be undef

my $script_prefix = script_prefix();
my $root = abs_path(catdir($CASE_DIR, '..', '..', '..'));
my $script_agat = $script_prefix . catfile($root, 'bin', 'agat');
my $script = $script_prefix . catfile($root, 'bin', 'agat_convert_sp_gxf2gxf.pl');
my $expected_output_path = catdir($CASE_DIR, 'out');
my $input_path = catdir($CASE_DIR, 'in');

my $dir = setup_tempdir();
my $pathtmp = catfile($dir, 'tmp.gff');

if ( $file =~ m/\.gff$/ ) {
    if ($file =~ m/^8_/ or $file =~ m/^33_/ or $file =~ m/^34_/ or $file =~ m/^36_/) {
        system("$script --gff " . catfile($input_path, $file) . " -o $pathtmp  2>&1 1>/dev/null");
    }
    elsif ($file =~ m/^28_/ or $file =~ m/^45_/ or $file =~ m/^46_/) {
        system("$script_agat config --expose --locus_tag Name 2>&1 1>/dev/null");
        system("$script --gff " . catfile($input_path, $file) . " -o $pathtmp  2>&1 1>/dev/null");
    }
    else {
        system("$script_agat config --expose --merge_loci 2>&1 1>/dev/null");
        system("$script --gff " . catfile($input_path, $file) . " -o $pathtmp  2>&1 1>/dev/null");
    }
}
else {
    system("$script --gff " . catfile($input_path, $file) . " -o $pathtmp 2>&1 1>/dev/null");
}

my $correct_output;
if (defined $expected_override) {
    $correct_output = catfile($expected_output_path, $expected_override);
} else {
    my ($num) = $file =~ /^(\d+)_/;
    $correct_output = catfile($expected_output_path, "${num}_correct_output.gff");
}

check_diff($pathtmp, $correct_output, "parse $file");

1;
