#!/usr/bin/env perl

use strict;
use warnings;
use Test::More tests => 47;
use FindBin qw($Bin);
use lib "$Bin/lib";
use File::Spec::Functions qw(catdir catfile);
use Cwd qw(abs_path);
use AGAT::TestUtilities qw(setup_tempdir check_diff script_prefix);

=head1 DESCRIPTION

Test to verify the parser deals properly with the different flavor / bugged gff files.

=cut

# Check if has to be run in Devel::Cover or not
my $script_prefix = script_prefix();

# script to call to check the parser
my $root = abs_path(catdir($Bin, '..'));
my $script_agat = $script_prefix . catfile($root, 'bin', 'agat');
my $script = $script_prefix . catfile($root, 'bin', 'agat_convert_sp_gxf2gxf.pl');
my $expected_output_path = catdir($Bin, 'gff_syntax', 'out');
my $input_path = catdir($Bin, 'gff_syntax', 'in');

# remove config in local folder if exists

# Loop over test
my $dir = $input_path; # folder where the test files are
opendir my $dh, $dir or die "Could not open '$dir' for reading: $!\n";
my @files = readdir $dh;
foreach my $file (sort { (($a =~ /^(\d+)/)[0] || 0) <=> (($b =~ /^(\d+)/)[0] || 0) } @files) {

  # skip few cases to be faster
  #my ($value) = $file =~ m/^(\d+)_.*/;
  #my $analyzed = 44;
  #if(defined($value) and $value < $analyzed){next;}
  

  # for all test files
  if ( $file =~ m/test.gff$/ ){

    # skip cases
    if ( ($file =~ m/^24_*/) or  ($file =~ m/^26_*/) ){ #Skip those cases
        next;
    }

    my $dir = setup_tempdir();
    my $pathtmp = catfile($dir, 'tmp.gff');

    # case do not merge loci 8,32,34,36
    if ($file =~ m/^8_/ or $file =~ m/^33_/ or $file =~ m/^34_/ or $file =~ m/^36_/){
        system("$script --gff " . catfile($input_path, $file) . " -o $pathtmp  2>&1 1>/dev/null");
    }
                # peculiar cases with locus_tag Name
    elsif($file =~ m/^28_/ or $file =~ m/^45_/ or $file =~ m/^46_/){
        system("$script_agat config --expose --locus_tag Name 2>&1 1>/dev/null"); # set special config for the test
        system("$script --gff " . catfile($input_path, $file) . " -o $pathtmp  2>&1 1>/dev/null");
    }
    # standard cases
    else{
        system("$script_agat config --expose --merge_loci 2>&1 1>/dev/null"); # set special config for the test
        system("$script --gff " . catfile($input_path, $file) . " -o $pathtmp  2>&1 1>/dev/null");
    }

    my @splitname = split /_/, $file;
    my $correct_output = catfile($expected_output_path, "$splitname[0]_correct_output.gff");

    #run test
    check_diff( $pathtmp, $correct_output, "parse $file" );

    # cleaning
    #if($value and $value == $analyzed){exit;}
  }
}
closedir($dh);

#  ---------------------------- special tests ----------------------------
{
    my $dir = setup_tempdir();
    my $pathtmp = catfile($dir, 'tmp.gff');
    my $result = catfile($expected_output_path, 'stop_start_an_exon_correct_output.gff');
    system(" $script --gff " . catfile($input_path, 'stop_start_an_exon.gtf') . " -o $pathtmp 2>&1 1>/dev/null");
    check_diff( $pathtmp, $result, "output stop_start_an_exon" );
}

{
    my $dir = setup_tempdir();
    my $pathtmp = catfile($dir, 'tmp.gff');
    my $result = catfile($expected_output_path, 'stop_split_over_two_exons_correct_output.gff');
    system(" $script --gff " . catfile($input_path, 'stop_split_over_two_exons.gtf') . " -o $pathtmp 2>&1 1>/dev/null");
    check_diff( $pathtmp, $result, "output stop_split_over_two_exons_correct_output" );
}

