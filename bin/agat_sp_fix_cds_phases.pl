#!/usr/bin/env perl

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use Bio::DB::Fasta;
use AGAT::AGAT;

start_script();
my $header = get_agat_header();

# ---------------------------- OPTIONS ----------------------------
my $opt_fasta = undef;
my $opt_gff;
my $opt_output;
my $opt_help = 0;

# OPTION MANAGEMENT: partition @ARGV into shared vs script options via library
my ($shared_argv, $script_argv) = split_argv_shared_vs_script(\@ARGV);

# Parse script-specific options from its own list
my $script_parser = Getopt::Long::Parser->new;
$script_parser->configure('bundling','no_auto_abbrev');
if ( ! $script_parser->getoptionsfromarray(
    $script_argv,
    'g|gff=s'          => \$opt_gff,
    'o|output=s'       => \$opt_output,
    'fasta|fa=s'       => \$opt_fasta,
    'h|help!'          => \$opt_help )  )
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

if (! defined($opt_gff) or ! defined($opt_fasta)){
    pod2usage( {
           -message => "$header\nAt least 2 parameters are mandatory:\nInput reference gff file (-g) and Input fasta file (--fasta).\n\n".
           "Ouptut is optional. Look at the help documentation to know more.\n",
           -verbose => 0,
           -exitval => 1 } );
}

# --- Manage config ---
my $shared_opts = parse_shared_options($shared_argv);
initialize_agat({ config_file_in => ( $shared_opts->{config} ), input => $opt_gff, shared_opts => $shared_opts });

# ----------------------------------------------------------------------------

######################
# Manage output file #
my $gffout = prepare_gffout( $opt_output );

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>     MAIN     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

######################
### Parse GFF input #
my ($hash_omniscient) = slurp_gff3_file_JD({ input => $opt_gff });

####################
# index the genome #
my $db = Bio::DB::Fasta->new($opt_fasta);
dual_print1 "Fasta file parsed\n";

###
# Fix frame
fil_cds_frame($hash_omniscient, $db);

###
# Print result
print_omniscient( {omniscient => $hash_omniscient, output => $gffout} );

# --- final messages ---
end_script();
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

=item B<-o> or B<--output>

Output GFF file.  If no output file is specified, the output will be
written to STDOUT.


=item B<-h> or B<--help>

Display this helpful text.

=back

=head1 SHARED OPTIONS

Shared options are defined in the AGAT configuration file and can be overridden via the command line for this script only.
Common shared options are listed below; for the full list, please refer to the AGAT agat_config.yaml.

=over 8

=item B<--config>

String - Path to a custom AGAT configuration file.  
By default, AGAT uses `agat_config.yaml` from the working directory if present, otherwise the default file shipped with AGAT
(available locally via `agat config --expose`).

=item B<--cpu>, B<--core>, B<--job> or B<--thread>

Integer - Number of parallel processes to use for file input parsing (via forking).

=item B<-v> or B<--verbose>

Integer - Verbosity, choice are 0,1,2,3,4. 0 is quiet, 1 is normal, 2,3,4 is more verbose. Default 1.

=back

=head1 FEEDBACK

For questions, suggestions, or general discussions about AGAT, please use the AGAT community forum:
https://github.com/NBISweden/AGAT/discussions

=head1 BUG REPORTING

Bug reports should be submitted through the AGAT GitHub issue tracker:
https://github.com/NBISweden/AGAT/issues

=cut

AUTHOR - Jacques Dainat
