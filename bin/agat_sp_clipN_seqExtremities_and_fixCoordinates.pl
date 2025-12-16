#!/usr/bin/env perl

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use Bio::SeqIO ;
use Bio::DB::Fasta;
use File::Basename;
use AGAT::AGAT;

start_script();
my $header = get_agat_header();
# -------------------------------- LOAD OPTIONS --------------------------------
my $opt_gfffile;
my $opt_fastafile;
my $opt_output_fasta;
my $opt_output_gff;
my $opt_help;
my $width = 60; # line length printed

my @copyARGV=@ARGV;
# OPTION MANAGEMENT: partition @ARGV into shared vs script options via library
my ($shared_argv, $script_argv) = split_argv_shared_vs_script(\@ARGV);

# Parse script-specific options from its own list
my $script_parser = Getopt::Long::Parser->new;
$script_parser->configure('bundling','no_auto_abbrev');
if ( ! $script_parser->getoptionsfromarray(
        $script_argv,
        'g|gff=s'        => \$opt_gfffile,
        'f|fa|fasta=s'   => \$opt_fastafile,
        'of=s'           => \$opt_output_fasta,
        'og=s'           => \$opt_output_gff,
        'h|help!'        => \$opt_help,
    ) ) {
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

# Parse shared options
my ($shared_opts) = parse_shared_options($shared_argv);

# --- Manage config ---
initialize_agat({ config_file_in => ($shared_opts->{config}), input => $opt_gfffile, shared_opts => $shared_opts });
# -----------------------------------------------------------------------------------------------

######################
# Manage output file #
my $ostream;
if ($opt_output_fasta) {
  open(my $fh, '>', $opt_output_fasta) or die "Could not open file $opt_output_fasta $!";
  $ostream= Bio::SeqIO->new(-fh => $fh, -format => 'Fasta' );
}
else{
  $ostream = Bio::SeqIO->new(-fh => \*STDOUT, -format => 'Fasta');
}

my $gffout = prepare_gffout( $opt_output_gff );

##### MAIN ####
#### read gff file and save info in memory
######################
### Parse GFF input #
my ($hash_omniscient) = slurp_gff3_file_JD({ input => $opt_gfffile });
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
print_omniscient( {omniscient => $hash_omniscient, output => $gffout} );

dual_print1 "We found $cpt_Nleft sequence(s) starting with N\n";
dual_print1 "We found $cpt_Nright sequence(s) ending with N\n";
dual_print1 "We found $cpt_Nboth sequence(s) having N both extremities\n";

# --- final messages ---
end_script();


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

=item B<-f>, B<--fa> or B<--fasta>

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
