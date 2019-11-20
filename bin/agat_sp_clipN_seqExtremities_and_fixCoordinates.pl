#!/usr/bin/env perl

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use Bio::SeqIO ;
use Bio::DB::Fasta;
use Bio::Tools::GFF;
use File::Basename;
use AGAT::Omniscient;

my $header = get_agat_header();
my $start_run = time();
my $opt_gfffile;
my $opt_fastafile;
my $opt_output_fasta;
my $opt_output_gff;
my $opt_help;
my $width = 60; # line length printed

# OPTION MANAGMENT
my @copyARGV=@ARGV;
if ( !GetOptions( 'g|gff=s' => \$opt_gfffile,
                  'f|fa|fasta=s' => \$opt_fastafile,
                  'of=s'      => \$opt_output_fasta,
                  'og=s'      => \$opt_output_gff,
                  'h|help!'         => \$opt_help ) )
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

if ( (! (defined($opt_gfffile)) ) or (! (defined($opt_fastafile)) ) ){
    pod2usage( {
           -message => "$header\nAt least 2 parametes are mandatory:\nInput reference gff file (-g);  Input reference fasta file (-f)\n\n".
           "Output is optional. Look at the help documentation to know more.\n",
           -verbose => 0,
           -exitval => 2 } );
}


my $ostream;
if ($opt_output_fasta) {
  $opt_output_fasta=~ s/.fasta//g;
  $opt_output_fasta=~ s/.fa//g;
  open(my $fh, '>', $opt_output_fasta.".fa") or die "Could not open file '$opt_output_fasta' $!";
  $ostream= Bio::SeqIO->new(-fh => $fh, -format => 'Fasta' );
}
else{
  $ostream = Bio::SeqIO->new(-fh => \*STDOUT, -format => 'Fasta');
}

my $gffout;
if ($opt_output_gff) {
  my($opt_output_gff, $dirs, $suffix) = fileparse($opt_output_gff, (".gff",".gff1",".gff2",".gff3",".gtf",".gtf1",".gtf2",".gtf3",".txt")); #remove extension
  open(my $fh, '>', $opt_output_gff.".gff3") or die "Could not open file '$opt_output_gff' $!";
  $gffout= Bio::Tools::GFF->new(-fh => $fh, -gff_version => 3 );
  }
else{
  $gffout = Bio::Tools::GFF->new(-fh => \*STDOUT, -gff_version => 3);
}

##### MAIN ####
#### read gff file and save info in memory
######################
### Parse GFF input #
print "Reading file $opt_gfffile\n";
my ($hash_omniscient, $hash_mRNAGeneLink) = slurp_gff3_file_JD({ input => $opt_gfffile
                                                              });
print "Parsing Finished\n";
### END Parse GFF input #
#########################

my $hash_l1_grouped = group_l1IDs_from_omniscient($hash_omniscient);

#### read fasta
my $nbFastaSeq=0;
my $db = Bio::DB::Fasta->new($opt_fastafile);
my @ids      = $db->get_all_primary_ids;
my %allIDs; # save ID in lower case to avoid cast problems
foreach my $id (@ids ){$allIDs{lc($id)}=$id;}


my $cpt_Nleft=0;
my $cpt_Nright=0;
my $cpt_Nboth=0;

foreach my $seq_id (@ids ){
  my $seqObject     = $db->get_Seq_by_id($seq_id);
  my $seq = $seqObject->seq;

  my @letters = split (//,$seq);

  ################
  # look at N at the beginning of the sequence
  my $nb_N_start = 0;
  foreach my $letter (@letters){

    if ( lc($letter) eq 'n'){
      $nb_N_start++
    }
    else{
      last;
    }
  }

  #start by N, let's remove them
  if ($nb_N_start != 0){
    $seq = substr $seq, $nb_N_start;
    shift_annotation($hash_omniscient, $hash_l1_grouped->{$seq_id}, $nb_N_start);
  }


  ####################
  # look at N at the end of the sequence
  my $nb_N_end = 0;
  foreach my $letter (reverse (@letters ) ){
    if ( lc($letter) eq 'n'){
      $nb_N_end++
    }
    else{
      last;
    }
  }

  ##############
  # CLIP Ns
  #start by N, let's remove them
  if ($nb_N_end != 0){
    $seq = substr $seq, 0, -$nb_N_end; # -0 will remove nothing at the end
  }

  if ($nb_N_end != 0 or $nb_N_start != 0){
    #create sequence object
    my $header = $db->header($seq_id);
    my $seqObj  = Bio::Seq->new( '-format' => 'fasta' , -seq => $seq, -id => $header );

    # print sequence object
    $ostream->write_seq($seqObj);
  }
  else{
    $ostream->write_seq($seqObject); #original object not modified
  }

  #####################
  #Handle counter to resume information
  if($nb_N_end != 0){
    $cpt_Nright++;
  }
  if($nb_N_start != 0){
    $cpt_Nleft++;
  }
    if($nb_N_end != 0 and $nb_N_start != 0){
    $cpt_Nboth++;
  }
}

# print annotation whith shifter location
print_omniscient($hash_omniscient, $gffout); #print gene modified

print "We found $cpt_Nleft sequence(s) starting with N\n";
print "We found $cpt_Nright sequence(s) ending with N\n";
print "We found $cpt_Nboth sequence(s) having N both extremities\n";

my $end_run = time();
my $run_time = $end_run - $start_run;
print "Job done in $run_time seconds\n";


########################################################################################

sub shift_annotation{
  my ($hash_omniscient, $list_id_l1, $nb_N_start) =@_;

  #Handle annotation
  foreach my $primary_tag_l1 (keys %{$hash_omniscient->{'level1'}}){ # primary_tag_l1 = gene or repeat etc...

    foreach my $id_tag_key_level1_raw (@$list_id_l1){
      my $id_tag_key_level1 = lc($id_tag_key_level1_raw);
      if(exists ($hash_omniscient->{'level1'}{$primary_tag_l1}{$id_tag_key_level1})){

        my $feature_level1 = $hash_omniscient->{'level1'}{$primary_tag_l1}{$id_tag_key_level1};
        # Shift position
        $feature_level1->start($feature_level1->start-$nb_N_start);
        $feature_level1->end($feature_level1->end-$nb_N_start);

        #################
        # == LEVEL 2 == #
        #################
        foreach my $primary_tag_l2 ( keys %{$hash_omniscient->{'level2'}}){ # primary_tag_l2 = mrna or mirna or ncrna or trna etc...

          if ( exists_keys( $hash_omniscient, ('level2', $primary_tag_l2, $id_tag_key_level1) ) ){
            foreach my $feature_level2 ( @{$hash_omniscient->{'level2'}{$primary_tag_l2}{$id_tag_key_level1}}) {
              # Shift position
              $feature_level2->start($feature_level2->start-$nb_N_start);
              $feature_level2->end($feature_level2->end-$nb_N_start);


              #################
              # == LEVEL 3 == #
              #################
              my $level2_ID = lc($feature_level2->_tag_value('ID'));

              ############
              # THEN ALL THE REST
              foreach my $primary_tag_l3 (keys %{$hash_omniscient->{'level3'}}){ # primary_tag_l3 = cds or exon or start_codon or utr etc...
                if ( exists_keys( $hash_omniscient, ('level3', $primary_tag_l3, $level2_ID) ) ){
                  foreach my $feature_level3 (@{$hash_omniscient->{'level3'}{$primary_tag_l3}{$level2_ID}}) {

                    # Shift position
                    $feature_level3->start($feature_level3->start-$nb_N_start);
                    $feature_level3->end($feature_level3->end-$nb_N_start);
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

__END__

=head1 NAME

agat_sp_clipN_seqExtremities_and_fixCoordinates.pl

=head1 DESCRIPTION

The script aims to clip the N's extremities of the sequences.
The annotation from the sequence clipped are modified accrodingly to stay consistent

=head1 SYNOPSIS

    agat_sp_clipN_seqExtremities_and_fixCoordinates.pl -g infile.gff -f infile.fasta  [ -o outfile ]
    agat_sp_clipN_seqExtremities_and_fixCoordinates.pl --help

=head1 OPTIONS

=over 8

=item B<-g>, B<--gff> or B<-ref>

Input GTF/GFF file.

=item B<-f> or B<--fasta>

Input fasta file.

=item B<--of>

Output fixed fasta file.  If no output file is specified, the output will be
written to STDOUT.

=item B<--og>

Output fixed GFF file.  If no output file is specified, the output will be
written to STDOUT

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
