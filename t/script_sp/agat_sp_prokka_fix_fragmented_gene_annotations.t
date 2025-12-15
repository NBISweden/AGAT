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

# ------------------- check agat_sp_prokka_fix_fragmented_gene_annotations script-------------------

my $script = $script_prefix . catfile($bin_dir, "agat_sp_prokka_fix_fragmented_gene_annotations.pl");

{ ok(system("$script -h 1>\/dev\/null") == 0, "help $script"); }

{
    my $dir = setup_tempdir();
    my $outtmp = catfile($dir, 'tmp');
    my $result_gff = "$output_folder/agat_sp_prokka_fix_fragmented_gene_annotations_1.gff";
    my $result_fa  = "$output_folder/agat_sp_prokka_fix_fragmented_gene_annotations_1.fa";
    check_quiet_run(
        "$script --gff " . catfile($input_folder,'prokka_cav_10DC88.gff') .
        " --fasta " . catfile($input_folder,'prokka_cav_10DC88.fa') .
        " --db " . catfile($input_folder,'prokka_bacteria_sprot.fa') .
        " --skip_hamap --frags -o $outtmp"
    );
    my $out_gff = catfile($outtmp,'prokka_cav_10DC88.gff');
    my $out_fa  = catfile($outtmp,'prokka_cav_10DC88.fa');
    check_diff($out_gff, $result_gff, 'output prokka_fix_fragmented gff');
    check_diff($out_fa,  $result_fa,  'output prokka_fix_fragmented fasta');
}

done_testing();