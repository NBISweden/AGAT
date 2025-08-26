#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Bio::DB::Fasta;
use IO::File;
use AGAT::AGAT;

my $header = get_agat_header();
my $start_run = time();

my ( $opt, $usage, $config ) = AGAT::AGAT::describe_script_options(
    $header,
    [ 'gff|file|input=s',  'Input reference gff file',   { required => 1 } ],
    [ 'fasta|f=s',         'Input reference fasta file', { required => 1 } ],
);

my $opt_gfffile   = $opt->gff;
my $opt_fastafile = $opt->fasta;
my $outfile       = $config->{output};
my $opt_verbose   = $config->{verbose};

my $log;
if ( my $log_name = $config->{log_path} ) {
    open( $log, '>', $log_name )
      or die "Can not open $log_name for printing: $!";
    dual_print( $log, $header,  3 );
}

# Manage input fasta file
my $format = $config->{force_gff_input_version};
if(! $format ){ $format = select_gff_format($opt_gfffile); }
my $ref_in = AGAT::BioperlGFF->new(-file => $opt_gfffile, -gff_version => $format);

# Manage Output
my $gffout = prepare_gffout($config, $outfile);

#### read fasta
my $nbFastaSeq=0;
my $db = Bio::DB::Fasta->new($opt_fastafile);
dual_print($log, "Fasta file parsed\n");

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
    warn "SequenceID ".$feature->seq_id." is absent from the fasta file\n" if $opt_verbose;
    $cpt_removed++;
  }
}

dual_print($log, "We removed $cpt_removed annotations.\n");
my $nbSeqWithAnnotation = scalar keys %seqNameSeen;
dual_print(
    $log,
    "We kept $cpt_kept annotations that are linked to $nbSeqWithAnnotation sequences.\n");
my $end_run = time();
my $run_time = $end_run - $start_run;
dual_print($log, "Job done in $run_time seconds\n");

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

=item B<-v> or B<--verbose>

For verbosity

=item B<-o> or B<--output>

STRING: Output file.  If no output file is specified, the output will be written to STDOUT. The result is in tabulate format.

=item B<-c> or B<--config>

String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any, 
otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

=item B<--help> or B<-h>

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
