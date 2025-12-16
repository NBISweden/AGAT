#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Pod::Usage;
use Getopt::Long;
use Bio::DB::Fasta;
use IO::File ;
use AGAT::AGAT;
start_script();
my $header = get_agat_header();
# ---------------------------- OPTIONS ----------------------------
my $opt_gfffile=undef;
my $opt_fastafile=undef;
my $outfile=undef;
my $opt_help = 0;

# OPTION MANAGEMENT: partition @ARGV into shared vs script options via library
my ($shared_argv, $script_argv) = split_argv_shared_vs_script(\@ARGV);

# Parse script-specific options from its own list
my $script_parser = Getopt::Long::Parser->new;
$script_parser->configure('bundling','no_auto_abbrev');
if ( ! $script_parser->getoptionsfromarray(
  $script_argv,
  'file|input|gff=s' => \$opt_gfffile,
  'f|fasta=s'        => \$opt_fastafile,
  'o|output=s'       => \$outfile,
  'h|help!'          => \$opt_help )  )
{
    pod2usage( { -message => "$header\nFailed to parse command line",
         -verbose => 1,
         -exitval => 1 } );
}

if ($opt_help) {
    pod2usage( { -verbose => 99,
                 -exitval => 0,
                 -message => "$header \n" } );
}

if ((!defined($opt_gfffile)) ){
   pod2usage( { -message => 'at least 2 parameters are mandatory',
                 -verbose => 0,
                 -exitval => 2 } );
}

# Parse shared options (CPU, config, etc.)
my ($shared_opts) = parse_shared_options($shared_argv);

# --- Manage config ---
initialize_agat({ config_file_in => $shared_opts->{config}, input => $opt_gfffile, shared_opts => $shared_opts });

# ------------------------------------------------------------------------------

# Manage input fasta file
my $format = $CONFIG->{force_gff_input_version};
if(! $format ){ $format = select_gff_format($opt_gfffile); }
my $inputfh = open_maybe_gz($opt_gfffile);
my $ref_in = AGAT::BioperlGFF->new(-fh => $inputfh, -gff_version => $format);

# Manage Output
my $gffout = prepare_gffout( $outfile );

#### read fasta
my $nbFastaSeq=0;
my $db = Bio::DB::Fasta->new($opt_fastafile);
dual_print1 "Fasta file parsed\n";

# get all seq id from fasta and convert to hash
my @ids      = $db->get_all_primary_ids;
my %hash_id;
$hash_id{ lc( $_ ) }++ for (@ids);

#time to calcul progression
my $startP=time;
my $cpt_removed=0;
my %seqNameSeen;
my $cpt_kept=0;
while (my $feature = $ref_in->next_feature() ) {
  my $seq_id = lc($feature->seq_id);
  if( exists_keys( \%hash_id, ( $seq_id ) ) ){
    $gffout->write_feature($feature);
    # to count number of sequence with annotation
    if(! exists_keys(\%seqNameSeen, ( $seq_id ) ) ){
      $seqNameSeen{$seq_id}++;
    }
    $cpt_kept++;
  }
  else{
    dual_print2 "SequenceID ".$feature->seq_id." is absent from the fasta file\n";
    $cpt_removed++;
  }
}

dual_print1 "We removed $cpt_removed annotations.\n";
my $nbSeqWithAnnotation = scalar keys %seqNameSeen;
dual_print1 "We kept $cpt_kept annotations that are linked to $nbSeqWithAnnotation sequences.\n";

# print fasta in asked and any
write_fasta($gffout, $ref_in);

# --- final messages ---
end_script();

# ---------------------------- FUNCTIONS ----------------------------

__END__

=head1 NAME

agat_sq_filter_feature_from_fasta.pl

=head1 DESCRIPTION

This script is a kind of annotation filter by sequence name.
It goes through the gff annotation features and remove those that are not linked to a sequence from the fasta file provided.
The match between sequence name in the fasta file and the 1st column of the gff3 file is case sensitive.

=head1 SYNOPSIS

    agat_sq_filter_feature_from_fasta.pl --gff <gff_file.gff> --fasta <fasta_file.fa> [-o <output file>]
    agat_sq_filter_feature_from_fasta.pl --help

=head1 OPTIONS

=over 8

=item B<--gff>, B<--file> or B<--input>

STRING: Input GTF/GFF file.

=item B<-f> or B<--fasta>

STRING: fasta file.


=item B<-o> or B<--output>

STRING: Output file.  If no output file is specified, the output will be written to STDOUT. The result is in tabulate format.

=item B<--help> or B<-h>

Display this helpful text.

=back

=head1 SHARED OPTIONS

Shared options are defined in the AGAT configuration file and can be overridden via the command line for this script only.
Common shared options are listed below; for the full list, please refer to the AGAT agat_config.yaml.
Note: For _sq_ scripts, only the following options are supported: verbose, output_format, gff_output_version, gtf_output_version, progress_bar, and tabix.

=over 8

=item B<--config>

String - Path to a custom AGAT configuration file.  
By default, AGAT uses `agat_config.yaml` from the working directory if present, otherwise the default file shipped with AGAT
(available locally via `agat config --expose`).

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
