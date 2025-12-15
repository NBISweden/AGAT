#!/usr/bin/env perl
use strict;
use warnings;
use File::Path;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use File::Spec::Functions qw(catfile catdir);
use Cwd qw(abs_path);
use AGAT::TestUtilities; 
use Test::More;

my $script_prefix = script_prefix();
my $root = abs_path(catdir($Bin, '..', '..'));
my $bin_dir = catdir($root, 'bin');
my $input_folder = catdir($Bin, 'in');
my $output_folder = catdir($Bin, 'out');
my $config = 'agat_config.yaml';

# -------------------------- check agat_sp_add_attribute_shortest_exon_size -------------------------

my $script = $script_prefix . catfile($bin_dir, "agat_sp_add_attribute_shortest_exon_size.pl");

{ ok(system("$script -h 1>\/dev\/null") == 0, "help $script"); }

{
    my $dir = setup_tempdir();
    my $outtmp_gff = catfile($dir,'tmp.gff');
    my $outtmp_report = catfile($dir,'tmp_report.txt');
    my $result = "$output_folder/agat_sp_add_attribute_shortest_exon_size.gff";
    check_quiet_run(
        "$script --gff " . catfile($input_folder,'1.gff') . " -o $outtmp_gff"
    );
    check_diff($outtmp_gff, $result, 'output add_attribute_shortest_exon_size gff');
}




done_testing();
