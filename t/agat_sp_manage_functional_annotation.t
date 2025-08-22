#!/usr/bin/env perl
use strict;
use warnings;
use File::Path;
use FindBin qw($Bin);
use lib "$Bin/lib";
use File::Spec::Functions qw(catfile catdir);
use Cwd qw(abs_path);
use AGAT::TestUtilities qw(setup_tempdir check_diff script_prefix check_quiet_run);
use Test::More;

my $script_prefix = script_prefix();
my $root = abs_path(catdir($Bin, '..'));
my $bin_dir = catdir($root, 'bin');
my $input_folder = catdir($Bin, 'scripts_output', 'in');
my $output_folder = catdir($Bin, 'scripts_output', 'out');
my $config = 'agat_config.yaml';

# --------check agat_sp_manage_functional_annotation.pl-------------

my $script = $script_prefix . catfile($bin_dir, "agat_sp_manage_functional_annotation.pl");
{ my $dir = setup_tempdir(); ok(system("$script -h 1>\/dev\/null") == 0, "help $script"); }

my $result = "$output_folder/agat_sp_manage_functional_annotation_1.gff";
{
    my $dir = setup_tempdir();
    my $outtmp = catfile($dir, 'tmp.gff');
    my $outprefix = catfile($dir, 'tmp');
    check_quiet_run(" $script --gff $input_folder/agat_sp_manage_functional_annotation/02413F.gff --db $input_folder/agat_sp_manage_functional_annotation/uniprot_sprot_test.fasta -b $input_folder/agat_sp_manage_functional_annotation/02413F_blast.out -i $input_folder/agat_sp_manage_functional_annotation/02413F_interpro.tsv --clean_name -o $outtmp");
    check_diff( "$outtmp/02413F.gff", $result, "output $script" );
}




done_testing();
