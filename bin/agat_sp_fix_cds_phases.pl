#!/usr/bin/env perl

use strict;
use warnings;
use Bio::DB::Fasta;
use AGAT::AGAT;

my $header    = get_agat_header();
my $start_run = time();
my ( $opt, $usage, $config ) = AGAT::AGAT::describe_script_options(
    $header,
    [ 'gff|g=s',   'Input reference gff file',   { required => 1 } ],
    [ 'fasta|fa=s', 'Input fasta file', { required => 1 } ],
);

my $opt_gfffile = $opt->gff;
my $opt_fasta   = $opt->fasta;
my $opt_output  = $config->{output};
my $opt_verbose = $config->{verbose};

my $log;
if ( my $log_name = $config->{log_path} ) {
    open( $log, '>', $log_name )
      or die "Can not open $log_name for printing: $!";
    dual_print( $log, $header, 0 );
}

######################
# Manage output file #
my $gffout = prepare_gffout($config, $opt_output);

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>     MAIN     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

######################
### Parse GFF input #
my ($hash_omniscient, $hash_mRNAGeneLink) = slurp_gff3_file_JD({ input => $opt_gfffile,
                                                                 config => $config
                                                            });
dual_print($log, "GFF3 file parsed\n", $opt_verbose);

####################
# index the genome #
my $db = Bio::DB::Fasta->new($opt_fasta);
dual_print($log, "Fasta file parsed\n", $opt_verbose);

###
# Fix frame
fil_cds_frame($hash_omniscient, $db, $opt_verbose);

###
# Print result
print_omniscient( {omniscient => $hash_omniscient, output => $gffout} );

my $end_run  = time();
my $run_time = $end_run - $start_run;
dual_print( $log, "Job done in $run_time seconds\n", $opt_verbose );

close $log if $log;
__END__

=head1 NAME

agat_sp_fix_cds_phases.pl

=head1 DESCRIPTION

This script aims to fix the CDS phases.
The script is compatible with incomplete gene models (Missing start, CDS
multiple of 3 or not, i.e. with offset of 1 or 2) and + and - strand.

How it works?  

AGAT uses the fasta sequence to verify the CDS frame.
In case the CDS start by a start codon the phase of the first CDS piece is set to 0.
In the case there is no start codon and: 
  - If there is only one stop codon in the sequence and it is located at the last position, the phase of the first CDS piece is set to 0.
  - If there is no stop codon, the phase of the first CDS piece is set to 0 (because sequence can be translated without premature stop codon).
  - If there is/are stop codon(s) in the middle of the sequence we re-execute the check with an offset of +2 nucleotides:
    - If there is only one stop codon in the sequence and it is located at the last position, the phase of the first CDS piece is set to 0.
    - If there is no stop codon, the phase of the first CDS piece is set to 0 (because sequence can be translated without premature stop codon).
    - If there is/are stop codon(s) in the middle of the sequence we re-execute the check with an offset of +1 nucleotide:
        - If there is only one stop codon in the sequence and it is located at the last position, the phase of the first CDS piece is set to 0.
        - If there is no stop codon, the phase of the first CDS piece is set to 0 (because sequence can be translated without premature stop codon).
        - If there is/are still stop codon(s) we keep original phase and throw a warning. In this last case it means we never succeded to make a translation without premature stop codon in all the 3 possible phases.
Then in case of CDS made of multiple CDS pieces (i.e. discontinuous feature), the rest of the CDS pieces will be checked accordingly to the first CDS piece.

What is a phase?  

For features of type "CDS", the phase indicates where the next codon begins
relative to the 5' end (where the 5' end of the CDS is relative to the strand
of the CDS feature) of the current CDS feature. For clarification the 5' end
for CDS features on the plus strand is the feature's start and and the 5' end
for CDS features on the minus strand is the feature's end. The phase is one of
the integers 0, 1, or 2, indicating the number of bases forward from the start
of the current CDS feature the next codon begins. A phase of "0" indicates that
a codon begins on the first nucleotide of the CDS feature (i.e. 0 bases forward),
a phase of "1" indicates that the codon begins at the second nucleotide of this
CDS feature and a phase of "2" indicates that the codon begins at the third
nucleotide of this region. Note that "Phase" in the context of a GFF3 CDS
feature should not be confused with the similar concept of frame that is also a
common concept in bioinformatics. Frame is generally calculated as a value for
a given base relative to the start of the complete open reading frame (ORF) or
the codon (e.g. modulo 3) while CDS phase describes the start of the next codon
relative to a given CDS feature.  
The phase is REQUIRED for all CDS features.

=head1 SYNOPSIS

    agat_sp_fix_cds_phases.pl --gff infile.gff -f fasta [ -o outfile ]
    agat_sp_fix_cds_phases.pl --help

=head1 OPTIONS

=over 8

=item B<-g>, B<--gff> or B<-ref>

Input GTF/GFF file.

=item B<-f>, B<--fa> or B<--fasta>

Input fasta file.

=item B<-v> or B<--verbose>

Add verbosity.

=item B<-o> or B<--output>

Output GFF file.  If no output file is specified, the output will be
written to STDOUT.

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
