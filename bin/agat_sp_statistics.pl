#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long qw(:config no_auto_abbrev);
use Pod::Usage;
use IO::File;
use Try::Tiny;
use List::MoreUtils qw(uniq);
use AGAT::AGAT;

my $header = get_agat_header();
my $config;
my $gff = undef;
my $opt_output = undef;
my $opt_yaml = undef;
my $opt_percentile = 90;
my $opt_genomeSize = undef;
my $opt_plot = undef;
my $opt_raw = undef;
my $opt_verbose = 0;
my $opt_help= 0;

if ( !GetOptions(
    'c|config=s'   => \$config,
    "h|help"       => \$opt_help,
    'o|output=s'   => \$opt_output,
    'percentile=i' => \$opt_percentile,
    'yaml!'        => \$opt_yaml,
    'r|raw!'       => \$opt_raw,
    'd|p!'         => \$opt_plot,
    'v|verbose'    => \$opt_verbose,
    'g|f|gs=s'     => \$opt_genomeSize,
    "gff|i=s"      => \$gff))

{
    pod2usage( { -message => "Failed to parse command line",
                 -verbose => 1,
                 -exitval => 1 } );
}

# Print Help and exit
if ($opt_help) {
    pod2usage( { -verbose => 99,
                 -exitval => 0,
                 -message => "$header\n" } );
}

if ( ! (defined($gff)) ){
    pod2usage( {
           -message => "$header\nAt least 1 parameter is mandatory:\nInput reference gff file (--gff) \n\n",
           -verbose => 0,
           -exitval => 1 } );
}

# --- Manage config ---
initialize_agat({ config_file_in => $config, input => $gff });

#### IN / OUT
my $out = prepare_fileout($opt_output);
if(defined($opt_yaml)){
  if( defined($opt_output)){
    $opt_yaml = $opt_output.".yaml";
  }
}

# Manage raw data
if($opt_raw){
  if ($opt_output){
    $opt_raw = $opt_output."_raw_data";
  }
  else{
    $opt_raw = "raw_data";

    if (-f $opt_raw){
      print "Cannot create a directory with the name $opt_raw because a file with this name already exists.\n";exit();
    }
    if (-d $opt_raw){
      print "Cannot create a directory with the name $opt_raw because a folder with this name already exists.\n";exit();
    }
  }
}

#Manage plot folder output
if($opt_plot){

  # Check if dependencies for plot are available
  if ( ! may_i_plot() ) {
    $opt_plot = undef;
  }
  else{
    if ($opt_output){
      $opt_plot = $opt_output."_distribution_plots";
    }
    else{
      $opt_plot = "distribution_plots";

      if (-f $opt_plot){
        print "Cannot create a directory with the name $opt_plot because a file with this name already exists.\n";exit();
      }
      if (-d $opt_plot){
        print "The default output directory $opt_plot use to save the distribution plots already exists. Please give me another folder name.\n";exit();
      }
    }
  }
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>     MAIN     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

######################
### Parse GFF input #
my ($hash_omniscient, $hash_mRNAGeneLink) =  slurp_gff3_file_JD({ input => $gff });
### END Parse GFF input #
#########################

##############
# STATISTICS #
print "Compute statistics\n";
print_omniscient_statistics ({ input   => $hash_omniscient,
															 genome  => $opt_genomeSize,
                               percentile => $opt_percentile,
															 output  => $out,
                               yaml    => $opt_yaml,
                               raw     => $opt_raw,
															 distri  => $opt_plot,
															 isoform => 1,
															 verbose => $opt_verbose
														 });
# END STATISTICS #
##################

print "Bye Bye.\n";

__END__

=head1 NAME

agat_sp_statistics.pl

=head1 DESCRIPTION

The script provides exhaustive statistics of a gft/gff file.
/!\ If you have isoforms in your file, even if correct, some values calculated
might sounds incoherent: e.g. total length mRNA can be superior than the genome size.
Because all isoforms length is added... It is why by default
we always compute the statistics twice when there are isoforms, once with the
isoforms, once without (In that case we keep the longest isoform per locus).

=head1 SYNOPSIS

    agat_sp_statistics.pl --gff file.gff  [ -o outfile ]
    agat_sp_statistics.pl --help

=head1 OPTIONS

=over 8

=item B<--gff> or B<-i>

Input GTF/GFF file.

=item B<--gs>, B<-f> or B<-g>

This option inform about the genome size in oder to compute more statistics. You can give the size in Nucleotide or directly the fasta file.

=item B<-d> or B<-p>

Bolean - When this option is used, an histogram of distribution of the features will be printed in pdf files in a folder with distribution_plots suffix. (d means distribution, p means plot).

=item B<-v> or B<--verbose>

Verbose option. To modify verbosity. Default is 1. 0 is quiet, 2 and 3 are increasing verbosity.

=item B<--output> or B<-o>

File where will be written the result. If no output file is specified, the output will be written to STDOUT.

=item B<--percentile>

Integer - Percentile to compute. Default is 90.

=item B<-r> or B<--raw>

Bolean - When this option is used, the raw data (same as used to create histogram of distribution of the features) are printed in a dedicated folder with raw_data suffix.

=item B<--yaml>

Bolean - When this option is activated, a second output will be printed either in STDOUT if no output provided or in <output.yaml> (a .yaml suffix is added to the --output value provided).

=item B<-c> or B<--config>

String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any, 
otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

=item B<-h> or B<--help>

Display this helpful text.

=back

=head1 FEEDBACK

=head2 Did you find a bug?

Do not hesitate to report bugs to help us keep track of the bugs and their
resolution. Please use the GitHub issue tracking system available at this
address:

            https://github.com/NBISweden/AGAT/issues

 Ensure that the bug was not already reported by searching under Issues.
 If you're unable to find an (open) issue addressing the problem, open a new one.
 Try as much as possible to include in the issue when relevant:
 - a clear description,
 - as much relevant information as possible,
 - the command used,
 - a data sample,
 - an explanation of the expected behaviour that is not occurring.

=head2 Do you want to contribute?

You are very welcome, visit this address for the Contributing guidelines:
https://github.com/NBISweden/AGAT/blob/master/CONTRIBUTING.md

=cut

AUTHOR - Jacques Dainat
