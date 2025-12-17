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

# --------check agat_sp_filter_incomplete_gene_coding_models.pl-------------

my $script = $script_prefix . catfile($bin_dir, "agat_sp_filter_incomplete_gene_coding_models.pl");

{ ok(system("$script -h 1>\/dev\/null") == 0, "help $script"); }

{
  my $dir = setup_tempdir();
  my $outtmp = catfile($dir,'tmp.gff');
  check_quiet_run(
    "$script --gff " . catfile($input_folder,'1.gff') .
    " --fasta " . catfile($input_folder,'1.fa') .
    " -o $outtmp"
  );
  my $out_complete   = catfile('tmp.gff');
  my $out_incomplete = catfile('tmp_incomplete.gff');
  my $result_complete   = "$output_folder/agat_sp_filter_incomplete_gene_coding_models_1.gff";
  my $result_incomplete = "$output_folder/agat_sp_filter_incomplete_gene_coding_models_incomplete_1.gff";
  check_diff($out_complete,   $result_complete,   'filter incomplete models complete');
  check_diff($out_incomplete, $result_incomplete, 'filter incomplete models incomplete');
}

done_testing();