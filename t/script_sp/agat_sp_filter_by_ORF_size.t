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

# --------check agat_sp_filter_by_ORF_size.pl-------------
my $script = $script_prefix . catfile($bin_dir, "agat_sp_filter_by_ORF_size.pl");

{ ok(system("$script -h 1>\/dev\/null") == 0, "help $script"); }


{
    my $dir = setup_tempdir();
    my $result = "$output_folder/agat_sp_filter_by_ORF_size_sup100.gff";
    # The script creates suffixed outputs from the -o prefix
    my $prefix = catfile($dir,'tmp.gff');
    my $outpass = catfile($dir,'tmp_sup100.gff');
    check_quiet_run(
        "$script --gff " . catfile($input_folder,'1.gff') .
        " -o $prefix"
    );
    check_diff($outpass, $result, 'output agat_sp_filter_by_ORF_size.pl sup100');
}

done_testing();