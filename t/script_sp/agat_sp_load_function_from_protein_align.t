#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use File::Spec::Functions qw(catdir catfile);
use Cwd qw(abs_path);
use AGAT::TestUtilities; 
use Test::More;

my $script_prefix = script_prefix();
my $root = abs_path(catdir($Bin, '..', '..'));
my $bin_dir = catdir($root, 'bin');

# --------check agat_sp_load_function_from_protein_align.pl-------------

my $script = $script_prefix . catfile($bin_dir, 'agat_sp_load_function_from_protein_align.pl');

{ system("$script -h 1>/dev/null"); ok(1, "help $script"); }

done_testing();
