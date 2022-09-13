#!/usr/bin/env perl

# script similar to agat_sp_gxf_to_gff3.pl

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use File::Basename;
use Bio::Tools::GFF;
use AGAT::Omniscient;

my $header = get_agat_header();
my $config = get_agat_config();
my $start_run = time();
my $opt_gfffile;
my $opt_verbose = 1;
my $opt_output;
my $opt_help;

# OPTION MANAGMENT
my @copyARGV=@ARGV;
if ( !GetOptions( 'g|gxf|gtf|gff=s'          => \$opt_gfffile,
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

This script fixes and/or standardizes any GTF/GFF file into full sorted GFF3 file.
The output GFF syntax is shaped by bioperl and choose among the versions
1,2,2.5 (GTF equivalent) and 3. For a correct GTF file, it is recommended to use
agat_convert_sp_gff2gtf.pl

Without specifying an input GTF/GFF version, the Omniscient parser will first detect
automtically the most appropriate GFF parser to use from bioperl (GFF1,GFF2,GFF3)
in order to read you file properly.
Then the Omniscient parser removes duplicate features, fixes duplicated IDs,
adds missing ID and/or Parent attributes, deflates factorized attributes
(attributes with several parents are duplicated with uniq ID), add missing features
when possible (e.g. add exon if only CDS described, add UTR if CDS and exon described),
fix feature locations (e.g. check exon is embedded in the parent features mRNA, gene), etc...
All AGAT's scripts with the _sp_ prefix use the same parser, before to perform supplement tasks.
With that script you can tuned the Omniscient parser behaviour. I.e. you can decide
to merge loci that have an overlap at their CDS features (Only one top feature
is kept (gene), and the mRNA features become isoforms). This is not activated by
default in case you are working on a prokaryote annotation that often have overlaping
loci.
The Omniscient parser defines relationship between features using 3 levels.
e.g Level1=gene; Level2=mRNA,tRNA; Level3=exon,cds,utr.
The feature type information is stored within the 3rd column of a GTF/GFF file.
The parser need to know to which level a feature type is part of. This information
is stored by default in a json file coming with the tool. We have implemented the
most common feature types met in gff/gtf files. If a feature type is not yet handle
by the parser it will throw a warning. You can easily inform the parser how
to handle it (level1, level2 or level3) by modifying the appropriate json file.
How to access the json files? Easy just use the --expose option and the json files
will appear in the working folder. By default, the Omniscient parser use
the json files from the working directory when any.

Omniscient parser phylosophy:

 Parse by Parent/child relationship
   ELSE Parse by a comomn tag  (an attribute value shared by feature that must be grouped together.
        By default we are using locus_tag and gene_id as locus tag, but you can specify the one of your choice
     ELSE Parse sequentially (features are grouped in a bucket, and the bucket change at each level2 feature met, and bucket(s) are linked to the first l1 top feature met)

=head1 SYNOPSIS

    agat_convert_sp_gxf2gxf.pl -g infile.gff [ -o outfile ]
    agat_convert_sp_gxf2gxf.pl --help

=head1 OPTIONS

=over 8

=item B<-g>, B<--gtf>, B<--gff> or B<--gxf>

String - Input GTF/GFF file. Compressed file with .gz extension is accepted.

=item B<-v>

Integer - Verbose option. To modify verbosity. Default is 1. 0 is quiet, 2 and 3 are increasing verbosity.

=item B<-o> or B<--output>

String - Output GFF file.  If no output file is specified, the output will be
written to STDOUT.

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
