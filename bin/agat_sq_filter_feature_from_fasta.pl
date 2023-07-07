#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Pod::Usage;
use Getopt::Long;
use Bio::DB::Fasta;
use IO::File ;
use AGAT::AGAT;

my $header = get_agat_header();
my $config;
my $start_run = time();
my $opt_gfffile=undef;
my $verbose=undef;
my $opt_fastafile=undef;
my $outfile=undef;
my $opt_help = 0;


Getopt::Long::Configure ('bundling');
if ( !GetOptions ('file|input|gff=s' => \$opt_gfffile,
      'f|fasta=s' => \$opt_fastafile,
      'o|output=s' => \$outfile,
      'v|verbose!' => \$verbose,
      'c|config=s'               => \$config,
      'h|help!'         => \$opt_help )  )
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

# Manage input fasta file
my $format = $config->{gff_output_version};
if(! $format ){ $format = select_gff_format($opt_gfffile); }
my $ref_in = Bio::Tools::GFF->new(-file => $opt_gfffile, -gff_version => $format);

# Manage Output
my $gffout = prepare_gffout($config, $outfile);

#### read fasta
my $nbFastaSeq=0;
my $db = Bio::DB::Fasta->new($opt_fastafile);
print ("Fasta file parsed\n");

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
    print "SequenceID ".$feature->seq_id." is absent from the fasta file\n" if($verbose);
    $cpt_removed++;
  }
}

print "We removed $cpt_removed annotations.\n";
my $nbSeqWithAnnotation = scalar keys %seqNameSeen;
print "We kept $cpt_kept annotations that are linked to $nbSeqWithAnnotation sequences.\n";
my $end_run = time();
my $run_time = $end_run - $start_run;
print "Job done in $run_time seconds\n";

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
