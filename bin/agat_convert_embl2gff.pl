#!/usr/bin/env perl

## TO DO => Deal With sequences. Write the DNA sequence of the "source" primary tag within the output gff3

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use Bio::Tools::GFF;
use Bio::SeqIO;
use AGAT::Omniscient;

my $header = get_agat_header();
my $outfile = undef;
my $embl = undef;
my $primaryTags = undef;
my $discard = undef;
my $keep = undef;
my $help;

if( !GetOptions(
    "help" => \$help,
    "embl=s" => \$embl,
    "ptag|t=s" => \$primaryTags,
    "d|s" => \$discard,
    "k" => \$keep,
    "outfile|output|o|out|gff=s" => \$outfile))
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
my $gff_out;
if ($outfile) {
open(my $fh, '>', $outfile) or die "Could not open file '$outfile' $!";
  $gff_out= Bio::Tools::GFF->new(-fh => $fh, -gff_version => 3);
}
else{
  $gff_out = Bio::Tools::GFF->new(-fh => \*STDOUT, -gff_version => 3);
}

### Read gb input file.
my $embl_in = Bio::SeqIO->new(-file => $embl, -format => 'embl');


### MAIN ###

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
        $gff_out->write_feature($feat_obj);
    }
  }
}

__END__

=head1 NAME

gaas_converter_embl2gff.pl

=head1 DESCRIPTION

The script takes an EMBL file as input, and will translate it in gff format.

=head1 SYNOPSIS

    gaas_converter_embl2gff.pl --embl infile.embl [ -o outfile ]

=head1 OPTIONS

=over 8

=item B<--embl>

Input EMBL file that will be read

=item B<-primary_tag>, B<--pt>, B<-t>

List of "primary tag". Useful to discard or keep specific features.
Multiple tags must be coma-separated.

=item B<-d>

Means that primary tags provided by the option "prinary_tag" will be discarded.

=item B<-d>

Means that only primary tags provided by the option "prinary_tag" will be kept.

=item B<-o> , B<--output> , B<--out> , B<--outfile> or B<--gff>

Output GFF file. If no output file is specified, the output will be
written to STDOUT.

=item B<-h> or B<--help>

Display this helpful text.

=back

=head1 FEEDBACK

=head2 Did you find a bug?

Do not hesitate to report bugs to help us keep track of the bugs and their
resolution. Please use the GitHub issue tracking system available at this
address:

            https://github.com/NBISweden/GAAS/issues

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
https://github.com/NBISweden/GAAS/blob/master/CONTRIBUTING.md

=cut

AUTHOR - Jacques Dainat
