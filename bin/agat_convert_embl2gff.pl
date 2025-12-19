#!/usr/bin/env perl

## TO DO => Deal With sequences. Write the DNA sequence of the "source" primary tag within the output gff3

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use Bio::SeqIO;
use AGAT::AGAT;

start_script();
my $header = get_agat_header();
# -------------------------------- LOAD OPTIONS --------------------------------
my $outfile;
my $embl;
my $emblmygff3;
my $primaryTags;
my $discard;
my $keep;
my $help;

# OPTION MANAGEMENT: partition @ARGV into shared vs script options via library
my ($shared_argv, $script_argv) = split_argv_shared_vs_script(\@ARGV);

# Parse script-specific options from its own list
my $script_parser = Getopt::Long::Parser->new;
$script_parser->configure('bundling','no_auto_abbrev');
if( ! $script_parser->getoptionsfromarray(
  $script_argv,
  "h|help"                     => \$help,
  "embl=s"                     => \$embl,
  "primary_tag|pt|t=s"         => \$primaryTags,
  "d|discard!"                 => \$discard,
  "k|keep!"                    => \$keep,
  "emblmygff3!"                => \$emblmygff3,
  "o|out|output=s"             => \$outfile,
) )
{
  pod2usage( { -message => "Failed to parse command line\n$header",
         -verbose => 1,
         -exitval => 1 } );
}

# Print Help and exit
if ($help) {
  pod2usage( { -verbose => 99,
         -exitval => 0,
         -message => "$header\n" } );
}

if ( ! (defined($embl)) ){
  pod2usage( {
       -message => "$header\nMissing the --embl argument",
       -verbose => 0,
       -exitval => 1 } );
}

## Parse shared options (e.g., config, cpu) from shared_argv
my ($shared_opts) = parse_shared_options($shared_argv);

# --- Manage config ---
initialize_agat({ config_file_in => ($shared_opts->{config}) , input => $embl, shared_opts => $shared_opts });
my $throw_fasta=$CONFIG->{"throw_fasta"};

##################
# MANAGE OPTION  #
if($discard and $keep){
  print "Cannot discard and keep the same primary tag. You have to choose if you want to discard it or to keep it.\n";
}

### If primaryTags given, parse them:

my @listprimaryTags;
if ($primaryTags){
  @listprimaryTags= split(/,/, $primaryTags);

  if($discard){
    print "We will not keep the following primary tag:\n";
    foreach my $tag (@listprimaryTags){
      print $tag,"\n";
    }
  }
  elsif($keep){ # Attribute we have to replace by a new name
    print "We will keep only the following primary tag:\n";
    foreach my $tag (@listprimaryTags){
      print $tag,"\n";
    }
  }
  else{print "You gave a list of primary tag wihtout telling me what you want I do with. Discard them or keep only them ?\n";}
}


##################
# MANAGE OUTPUT  #
my $gff_out = prepare_gffout( $outfile);

### Read embl input file.
my $embl_in = Bio::SeqIO->new(-file => $embl, -format => 'embl');


### MAIN ###
my $sequence="";
while( my $seq_obj = $embl_in->next_seq) {

  for my $feat_obj ($seq_obj->get_SeqFeatures) {
    my $skipit=undef;

    # In case we should discard some
    if($discard){

      foreach my $pTag (@listprimaryTags){
        if(lc($pTag) eq lc($feat_obj->primary_tag)){
          $skipit=1;last;
        }
      }
    }
    # In case we should keep only some
    elsif($keep){
      my $skipit=1;
      foreach my $pTag (@listprimaryTags){
        if(lc($pTag) eq lc($feat_obj->primary_tag)){
          $skipit=undef;last;
        }
      }

    }

    if(! $skipit){

      if($emblmygff3){
        # ------ Get seqId name when EMBLmyGFF3 file -----
        # get the second AC line
        my @arr = $seq_obj->get_secondary_accessions;
        my $index = 0;
        $index++ until $arr[$index] eq '*';
        splice(@arr, $index, 1);
        my $seqid_raw = $arr[0];
        my $seqid_clean = substr $seqid_raw, 1; # remove the _ at the beginning added by EMBLmyGFF3
        $feat_obj->seq_id($seqid_clean); # replace the default seq_id
      }

      $gff_out->write_feature($feat_obj);
      $sequence.= ">".$seq_obj->seq();

    }
  }
}

# Close the gff input FH opened by OmniscientI
$embl_in->close();

# if user want to keep the fasta file
if (! $throw_fasta){

  ### Read embl input file to cach the fasta sequences now.
  $embl_in = Bio::SeqIO->new(-file => $embl, -format => 'embl');

  # Print sequences
  _write_fasta($gff_out, $embl_in, $emblmygff3);

  # Close the gff input FH opened by OmniscientI
  $embl_in->close();
}

# --- final messages ---
end_script();

#################################### methods ####################################

# Catch sequences from embl file and write all of them at the end of the gff file
sub _write_fasta {
  my ($gffout, $embl_in, $emblmygff3) = @_;

  $gffout->_print("##FASTA\n");

  while( my $Bio_Seq_obj = $embl_in->next_seq ) {

    my $seq_header = ">".$Bio_Seq_obj->display_id();

    if($emblmygff3){
      # ------ Get seqId name when EMBLmyGFF3 file -----
      # get the second AC line
      my @arr = $Bio_Seq_obj->get_secondary_accessions;
      my $index = 0;
      $index++ until $arr[$index] eq '*';
      splice(@arr, $index, 1);
      my $seqid_raw = $arr[0];
      my $seqid_clean = substr $seqid_raw, 1; # remove the _ at the beginning added by EMBLmyGFF3
      $seq_header=">".$seqid_clean; # replace the default seq_id
    }


    if( $Bio_Seq_obj->desc ){
      $gffout->_print($seq_header." ".$Bio_Seq_obj->desc."\n");
    }
    else{
      $gffout->_print($seq_header."\n");
    }

    my $str = $Bio_Seq_obj->seq;
    my $nuc = 80;       # Number of nucleotides per line
    my $length = length($str);

    # Calculate the number of nucleotides which fit on whole lines
    my $whole = int($length / $nuc) * $nuc;

    # Print the whole lines
    my( $i );
    for ($i = 0; $i < $whole; $i += $nuc) {
        my $blocks = substr($str, $i, $nuc);
        $gffout->_print("$blocks\n") ;
    }
    # Print the last line
    if (my $last = substr($str, $i)) {
        $gffout->_print("$last\n") ;
    }
  }
}

__END__

=head1 NAME

agat_converter_embl2gff.pl

=head1 DESCRIPTION

The script takes an EMBL file as input, and will translate it in gff format.

=head1 SYNOPSIS

    agat_converter_embl2gff.pl --embl infile.embl [ -o outfile ]

=head1 OPTIONS

=over 8

=item B<--embl> <file>

Input EMBL file that will be read

=item B<--emblmygff3>
Means that the EMBL flat file comes from the EMBLmyGFF3 software.
This is an EMBL format dedicated for submission and contains particularity to deal with.
This parameter is needed to get a proper sequence id in the GFF3 from an embl made with EMBLmyGFF3.

=item B<--primary_tag>, B<--pt> or B<-t> <list>

List of "primary tag". Useful to discard or keep specific features.
Multiple tags must be coma-separated.

=item B<-d> or B<--discard>

Means that primary tags provided by the option "primary_tag" will be discarded.

=item B<-k> or B<--keep>

Means that only primary tags provided by the option "primary_tag" will be kept.

=item B<-o>, B<--out> or B<--output> <file>

Output GFF file to create. If no output file is specified, the output will be
written to STDOUT.

=item B<-h> or B<--help>

Display this helpful text.


=item B<-v> or B<--verbose> <int>

Verbosity, choice are 0,1,2,3,4. 0 is quiet, 1 is normal, 2,3,4 is more verbose. Default 1.

=back

=head1 FEEDBACK

For questions, suggestions, or general discussions about AGAT, please use the AGAT community forum:
https://github.com/NBISweden/AGAT/discussions

=head1 BUG REPORTING

Bug reports should be submitted through the AGAT GitHub issue tracker:
https://github.com/NBISweden/AGAT/issues

=cut

AUTHOR - Jacques Dainat
