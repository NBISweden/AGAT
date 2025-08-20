#!/usr/bin/env perl

# script similar to agat_sp_gxf_to_gff3.pl

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use File::Basename;
use AGAT::AGAT;

my $header = get_agat_header();
my $config;
my $start_run = time();
my $opt_gfffile;
my $opt_output;
my $opt_help;

my $common = parse_common_options() || {};
$config     = $common->{config};
$opt_output = $common->{output};
$opt_help   = $common->{help};

# OPTION MANAGMENT
my @copyARGV=@ARGV;
if ( !GetOptions( 'g|gxf|gtf|gff=s'          => \$opt_gfffile,
                  'c|config=s'               => \$config,
                  'o|output=s'               => \$opt_output,
                  'h|help!'                  => \$opt_help ) )
{
    pod2usage( { -message => 'Failed to parse command line',
                 -verbose => 1,
                 -exitval => 1 } );
}

# Print Help and exit
if ($opt_help) {
    pod2usage( { -verbose => 99,
                 -exitval => 0,
                 -message => "$header\n" } );
}

if (! defined($opt_gfffile) ){
    pod2usage( {
           -message => "$header\nAt least 1 parameter is mandatory:\n --gff (Input reference gff file).\n\n".
           "Invoke the help for more information (--help).\n",
           -verbose => 0,
           -exitval => 1 } );
}

# --- Manage config ---
$config = get_agat_config({config_file_in => $config});

my $log;
my $log_name = get_log_path($common, $config);
open($log, '>', $log_name) or die "Can not open $log_name for printing: $!";
dual_print($log, $header, 0);

######################
# Manage output file #
my $gffout = prepare_gffout($config, $opt_output);

                #####################
                #     MAIN          #
                #####################

######################
### Parse GFF input #
my ($hash_omniscient, $hash_mRNAGeneLink) = slurp_gff3_file_JD({
                                                                input  => $opt_gfffile,
															    config => $config
																 });
print ("GFF3 file parsed\n");

###
# Print result
print_omniscient( {omniscient => $hash_omniscient,
                   output => $gffout
                } );

my $end_run = time();
my $run_time = $end_run - $start_run;
print "usage: $0 @copyARGV\n";
print "Job done in $run_time seconds\n";

__END__
=head1 NAME

agat_convert_sp_gxf2gxf.pl

=head1 DESCRIPTION

This script fixes and/or standardizes any GTF/GFF file into full sorted GTF/GFF file.
It AGAT parser removes duplicate features, fixes duplicated IDs, adds missing ID and/or Parent attributes,
deflates factorized attributes (attributes with several parents are duplicated with uniq ID),
add missing features when possible (e.g. add exon if only CDS described, add UTR if CDS and exon described),
fix feature locations (e.g. check exon is embedded in the parent features mRNA, gene), etc...

All AGAT's scripts with the _sp_ prefix use the AGAT parser, before to perform any supplementary task.
So, it is not necessary to run this script prior the use of any other _sp_ script. 

=head1 SYNOPSIS

    agat_convert_sp_gxf2gxf.pl -g infile.gff [ -o outfile ]
    agat_convert_sp_gxf2gxf.pl --help

=head1 OPTIONS

=over 8

=item B<-g>, B<--gtf>, B<--gff> or B<--gxf>

String - Input GTF/GFF file. Compressed file with .gz extension is accepted.

=item B<-o> or B<--output>

String - Output GFF file. If no output file is specified, the output will be
written to STDOUT.

=item B<-c> or B<--config>

String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any, 
otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

=item B<-h> or B<--help>

Boolean - Display this helpful text.

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
