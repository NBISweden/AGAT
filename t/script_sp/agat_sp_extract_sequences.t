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

# --------check agat_sp_extract_sequences.pl-------------

my $script = $script_prefix . catfile($bin_dir, "agat_sp_extract_sequences.pl");

{ ok(system("$script -h 1>\/dev\/null") == 0, "help $script"); }

{
    my $dir = setup_tempdir();
    my $outtmp = catfile($dir,'tmp1.fa');
    my $result = "$output_folder/agat_sp_extract_sequences/agat_sp_extract_sequences_1.fa";
    check_quiet_run(
        "$script --gff " . catfile($input_folder,'1.gff') .
        " --fasta " . catfile($input_folder,'1.fa') .
        " -o $outtmp"
    );
    check_diff($outtmp, $result, 'extract_sequences basic');
}

{
    my $dir = setup_tempdir();
    my $outtmp = catfile($dir,'tmp2.fa');
    my $result = "$output_folder/agat_sp_extract_sequences/agat_sp_extract_sequences_split.fa";
    check_quiet_run(
        "$script --gff " . catfile($input_folder,'1.gff') .
        " --fasta " . catfile($input_folder,'1.fa') .
        " --split -o $outtmp"
    );
    check_diff($outtmp, $result, 'extract_sequences split');
}

{
    my $dir = setup_tempdir();
    my $outtmp = catfile($dir,'tmp3.fa');
    my $result = "$output_folder/agat_sp_extract_sequences/agat_sp_extract_sequences_merge.fa";
    check_quiet_run(
        "$script --gff " . catfile($input_folder,'1.gff') .
        " --fasta " . catfile($input_folder,'1.fa') .
        " -t exon --merge -o $outtmp"
    );
    check_diff($outtmp, $result, 'extract_sequences merge exon');
}

{
    my $dir = setup_tempdir();
    my $outtmp = catfile($dir,'tmp4.fa');
    my $result = "$output_folder/agat_sp_extract_sequences/agat_sp_extract_sequences_full.fa";
    check_quiet_run(
        "$script --gff " . catfile($input_folder,'1.gff') .
        " --fasta " . catfile($input_folder,'1.fa') .
        " --full -o $outtmp"
    );
    check_diff($outtmp, $result, 'extract_sequences full');
}

{
    my $dir = setup_tempdir();
    my $outtmp = catfile($dir,'tmp5.fa');
    my $result = "$output_folder/agat_sp_extract_sequences/agat_sp_extract_sequences_attributes_kept.fa";
    check_quiet_run(
        "$script --gff " . catfile($input_folder,'1.gff') .
        " --fasta " . catfile($input_folder,'1.fa') .
        " --keep_attributes -o $outtmp"
    );
    check_diff($outtmp, $result, 'extract_sequences keep_attributes');
}

{
    my $dir = setup_tempdir();
    my $outtmp = catfile($dir,'tmp6.fa');
    my $result = "$output_folder/agat_sp_extract_sequences/agat_sp_extract_sequences_parent_attributes_kept.fa";
    check_quiet_run(
        "$script --gff " . catfile($input_folder,'1.gff') .
        " --fasta " . catfile($input_folder,'1.fa') .
        " --keep_parent_attributes -o $outtmp"
    );
    check_diff($outtmp, $result, 'extract_sequences keep_parent_attributes');
}

done_testing();