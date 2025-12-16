#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Pod::Usage;
use Getopt::Long;
use Bio::DB::Fasta;
use IO::File ;
use File::Basename;
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

if ((!defined($opt_gfffile) or !defined($opt_fastafile) ) ){
   pod2usage( { -message => "At least 2 parameters are mandatory:\n * A gff/gtf file (--gff)\n * A fasta file (--fasta)",
                 -verbose => 0,
                 -exitval => 2 } );
}

# Parse shared options (CPU, config, etc.)
my ($shared_opts) = parse_shared_options($shared_argv);

# --- Manage config ---
initialize_agat({ config_file_in => ( $shared_opts->{config} ), input => $opt_gfffile, shared_opts => $shared_opts });

# Manage input gff file
my $format = $CONFIG->{force_gff_input_version};
if(! $format ){ $format = select_gff_format($opt_gfffile); }
my $inputfh = open_maybe_gz($opt_gfffile);
my $ref_in = AGAT::BioperlGFF->new(-fh => $inputfh, -gff_version => $format);

# Manage output fasta file
my ($fasta_in,$path,$ext) = fileparse($opt_fastafile,qr/\.[^.]*/);
my $fasta_out = $fasta_in."_rt".$ext;
open(my $fh_fasta, '>', $fasta_out) or die "Could not open file '$fasta_out' $!";
$fasta_out = Bio::SeqIO->new(-fh => $fh_fasta , -format => 'Fasta');

# Manage Output
my $gffout = prepare_gffout( $outfile);

#### rt fasta
my $seqio = Bio::SeqIO->new(-file => $opt_fastafile, -format => "fasta");
while(my $seqObj = $seqio->next_seq) {
    $seqObj->revcom;
    $fasta_out->write_seq($seqObj);
}

#### read fasta again for DB
my $nbFastaSeq=0;
my $db = Bio::DB::Fasta->new($opt_fastafile);
dual_print1 "Fasta file parsed\n";

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
      $feature->strand($strand);
    }
    elsif ( ($strand == 1) or ($strand eq "+") ) {
      $strand = "-";
      $feature->strand($strand);
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

# --- final messages ---
dual_print1 "Annotations on $nb_flip sequences have been reverse complemented.\n";
dual_print1 "Annotations on $nb_intact sequences have been kept intact (Sequences absent from the fasta file but present in the gff.\n";
dual_print1 "$nb_error sequences from the fasta file were absent from the gff.\n";
end_script();

# ---------------------------- FUNCTIONS ----------------------------

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

=item B<-o> or B<--output>

STRING: Output file.  If no output file is specified, the output will be written to STDOUT.


=item B<--help> or B<-h>

BOOLEAN: Display this helpful text.

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
