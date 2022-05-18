#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Pod::Usage;
use Getopt::Long;
use Bio::DB::Fasta;
use IO::File ;
use Bio::Tools::GFF;
use File::Basename;
use AGAT::Omniscient;

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

if ((!defined($opt_gfffile) or !defined($opt_fastafile) ) ){
   pod2usage( { -message => "At least 2 parameters are mandatory:\n * A gff/gtf file (--gff)\n * A fasta file (--fasta)",
                 -verbose => 0,
                 -exitval => 2 } );
}

my $ostream     = IO::File->new();

# Manage input gff file
my $format = select_gff_format($opt_gfffile);
my $ref_in = Bio::Tools::GFF->new(-file => $opt_gfffile, -gff_version => $format);

# Manage output fasta file
my ($fasta_in,$path,$ext) = fileparse($opt_fastafile,qr/\.[^.]*/);
my $fasta_out = $fasta_in."_rt".$ext;
open(my $fh_fasta, '>', $fasta_out) or die "Could not open file '$fasta_out' $!";
$fasta_out = Bio::SeqIO->new(-fh => $fh_fasta , -format => 'Fasta');

# Manage Output
my $gffout;
if ($outfile) {
  open(my $fh, '>', $outfile) or die "Could not open file '$outfile' $!";
  $gffout= Bio::Tools::GFF->new(-fh => $fh, -gff_version => 3 );
}
else{
  $gffout = Bio::Tools::GFF->new(-fh => \*STDOUT, -gff_version => 3);
}


#### rt fasta 
my $seqio = Bio::SeqIO->new(-file => $opt_fastafile, -format => "fasta");
while(my $seqObj = $seqio->next_seq) {
    $seqObj->revcom;
    $fasta_out->write_seq($seqObj);
}

#### read fasta again for DB
my $nbFastaSeq=0;
my $db = Bio::DB::Fasta->new($opt_fastafile);
print ("Fasta file parsed\n");

# get all seq id from fasta and convert to hash
my @ids      = $db->get_all_primary_ids;
my %hash_id;
foreach my $id (@ids){
  my $length = $db->length($id);
  $hash_id{ lc( $id ) } = $length 
}

#time to calcul progression
my $startP=time;
my %info;
while (my $feature = $ref_in->next_feature() ) {
  my $seq_id = lc($feature->seq_id);
  if( exists_keys( \%hash_id, ( $seq_id ) ) ){

    # Flip location properly
    my $length_seq = $hash_id{$seq_id};
    my $start = $feature->start(); 
    my $end = $feature->end();    
    $feature->end($length_seq-$start+1); # set new end
    $feature->start($length_seq-$end+1); # set new start

    # Flip strand
    my $strand = $feature->strand;
    if ( ($strand == -1) or ($strand eq "-") ) {
      $strand = "+";
    }
    elsif ( ($strand == 1) or ($strand eq "+") ) {
      $strand = "-";
    }

    $gffout->write_feature($feature); 
    # to count number of sequence with annotation
    $info{"flip"}{$seq_id}++;
    $info{"all"}{$seq_id}++;
  }
  else{
    $info{"intact"}{$seq_id}++;
    $info{"all"}{$seq_id}++;
    $gffout->write_feature($feature);
  }
}

# Check fasta seq not used in gff
foreach my $seq_id (keys %hash_id){
  if(! exists_keys( \%info, ( "all", $seq_id ) ) ){
    $info{"error"}{$seq_id}++;
  }
}

my $nb_intact = 0;
$nb_intact =  keys %{$info{"intact"}};
my $nb_flip = 0;
$nb_flip   =  keys %{$info{"flip"}};
my $nb_error = 0;
$nb_error  =  keys %{$info{"error"}};


print "Annotations on $nb_flip sequences have been reverse complemented.\n";
print "Annotations on $nb_intact sequences have been kept intact (Sequences absent from the fasta file but present in the gff.\n";
print "$nb_error sequences from the fasta file were absent from the gff.\n";

my $end_run = time();
my $run_time = $end_run - $start_run;
print "Job done in $run_time seconds\n";

__END__

=head1 NAME

agat_sq_reverse_complement.pl

=head1 DESCRIPTION

This script will reverse complement the annotation of all annotation from the gff that are hold by sequences described in the fasta file.
The match between sequence name in the fasta file and the 1st column of the gff3 file is case sensitive.

=head1 SYNOPSIS

    agat_sq_reverse_complement.pl --gff <gff_file.gff> --fasta <fasta_file.fa> [-o <output file>]
    agat_sq_reverse_complement.pl --help

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
