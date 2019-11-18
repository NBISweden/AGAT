#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Pod::Usage;
use Getopt::Long;
use Bio::DB::Fasta;
use IO::File ;
use Bio::Tools::GFF;
use Agat::Omniscient;

my $header = get_agat_header();
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

my $ostream     = IO::File->new();

# Manage input fasta file
my $format = select_gff_format($opt_gfffile);
my $ref_in = Bio::Tools::GFF->new(-file => $opt_gfffile, -gff_version => $format);

# Manage Output
my $gffout;
if ($outfile) {
  open(my $fh, '>', $outfile) or die "Could not open file '$outfile' $!";
  $gffout= Bio::Tools::GFF->new(-fh => $fh, -gff_version => 3 );
}
else{
  $gffout = Bio::Tools::GFF->new(-fh => \*STDOUT, -gff_version => 3);
}


#### read fasta
my $nbFastaSeq=0;
my $db = Bio::DB::Fasta->new($opt_fastafile);
print ("Genome fasta parsed\n");

#time to calcul progression
my $startP=time;
my $cpt_removed=0;
my %seqNameSeen;
my $cpt_kept=0;
while (my $feature = $ref_in->next_feature() ) {

  if($db->seq($feature->seq_id)){
    $gffout->write_feature($feature);
    # to count number of sequence with annotation
    if(! exists_keys(\%seqNameSeen, ($feature->seq_id))){
      $seqNameSeen{$feature->seq_id}++;
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

agat_sq_keep_annotation_from_fastaSeq.pl

=head1 DESCRIPTION

This script is a kind of annotation filter by sequence name.
It goes through the gff annotation features and remove those that are not linked to a sequence from the fasta file provided.
The match between sequence name in the fasta file and the 1st column of the gff3 file is case sensitive.

=head1 SYNOPSIS

    agat_sq_keep_annotation_from_fastaSeq.pl --gff <gff_file.gff> --fasta <fasta_file.fa> [-o <output file>]
    agat_sq_keep_annotation_from_fastaSeq.pl --help

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

=cut
