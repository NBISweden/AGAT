#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Pod::Usage;
use Getopt::Long;
use IO::File ;
use Bio::Tools::GFF;
use AGAT::Omniscient;

my $header = get_agat_header();
my $start_run = time();
my $inputFile=undef;
my $outfile=undef;
my $outformat=undef;
my $opt_help = 0;

Getopt::Long::Configure ('bundling');
if ( !GetOptions ('file|input|gff|i=s' => \$inputFile,
      'of=i' => \$outformat,
      'o|output=s' => \$outfile,
      'h|help!'         => \$opt_help )  )
{
    pod2usage( { -message => 'Failed to parse command line',
                 -verbose => 1,
                 -exitval => 1 } );
}

if ($opt_help) {
    pod2usage( { -verbose => 99,
                 -exitval => 0,
                 -message => "$header\n" } );
}

if ((!defined($inputFile)) ){
   pod2usage( { -message => "$header\nAt least 1 parameter is mandatory: -i",
                 -verbose => 0,
                 -exitval => 1 } );
}


# Manage input fasta file
my $format = select_gff_format($inputFile);
my $ref_in = Bio::Tools::GFF->new(-file => $inputFile, -gff_version => $format);


# Manage Output
if(! $outformat){
  $outformat=$format;
}

my $gffout;
if ($outfile) {
  open(my $fh, '>', $outfile) or die "Could not open file '$outfile' $!";
  $gffout= Bio::Tools::GFF->new(-fh => $fh, -gff_version => $outformat );

}
else{
  $gffout = Bio::Tools::GFF->new(-fh => \*STDOUT, -gff_version => $outformat);
}

#time to calcul progression
my $startP=time;
my $nbLine=`wc -l < $inputFile`;
$nbLine =~ s/ //g;
chomp $nbLine;
print "$nbLine line to process...\n";

my $line_cpt=0;
my %hash_IDs;
my %featCount;
my %mapID;
while (my $feature = $ref_in->next_feature() ) {
  $line_cpt++;

  _uniq_ID ($feature, \%hash_IDs, \%featCount, \%mapID);

  if($feature->has_tag('Parent')){
    my $parent = lc($feature->_tag_value('Parent'));
    if(! exists($mapID{$parent})){
      print "How is it possible ? This parent hasn't been seen before\n";
    }
     create_or_replace_tag($feature,'Parent', $mapID{$parent});
  }

  $gffout->write_feature($feature);

  #####################
  #Display progression
  if ((30 - (time - $startP)) < 0) {
    my $done = ($line_cpt*100)/$nbLine;
    $done = sprintf ('%.0f', $done);
        print "\rProgression : $done % processed.\n";
    $startP= time;
  }
}

##Last round
my $end_run = time();
my $run_time = $end_run - $start_run;
print "Job done in $run_time seconds\n";



sub _uniq_ID{
  my ($feature, $hash_IDs, $miscCount, $mapID) = @_;


  my  $key=lc($feature->primary_tag);
  $miscCount->{$key}++;
  my $id = $key."-".$miscCount->{$key};

  while( exists_keys($hash_IDs, ($id) ) ){  #loop until we found an uniq tag
    $miscCount->{$key}++;
    $id = $key."-".$miscCount->{$key};
  }

  #push the new ID
  $hash_IDs->{$id}++;
  my $original_id = undef;
  if($feature->has_tag('ID')){
   $original_id = lc($feature->_tag_value('ID'));
  }
  else{
    $original_id = lc($id);
  }
  $mapID->{$original_id} = $id;

  # modify the feature ID with the correct one chosen
  create_or_replace_tag($feature,'ID', $id); #modify ID to replace by parent value
}

__END__

=head1 NAME

agat_sq_manage_ID.pl

=head1 DESCRIPTION

The script changes IDs to give uniq one and reflect the change in Parent attribute
of impacted features.

=head1 SYNOPSIS

    agat_sq_manage_ID.pl --gff <input file> [-o <output file>]
    agat_sq_manage_ID.pl --help

=head1 OPTIONS

=over 8

=item B<--gff>, B<--file> or B<--input>

STRING: Input GTF/GFF file.

=item B<--of>

Output format, if no ouput format is given, the same as the input one detected will be used. Otherwise you can force to have a gff version 1 or 2 or 3 by giving the corresponding number.

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
