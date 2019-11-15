#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Pod::Usage;
use Getopt::Long;
use IO::File ;
use Bio::Tools::GFF;

my $header = get_agat_header();
my $start_run = time();
my $inputFile;
my $outfile;
my $opt_help = 0;

Getopt::Long::Configure ('bundling');
if ( !GetOptions ('i|file|input|gff=s' => \$inputFile,
      'o|output=s' => \$outfile,
      'h|help!'         => \$opt_help )  )
{
    pod2usage( { -message => 'Failed to parse command line',
                 -verbose => 1,
                 -exitval => 1 } );
}

# Print Help and exit
if ($opt_help) {
    pod2usage( { -verbose => 2,
                 -exitval => 2,
                 -message => "$header\n" } );
}

if ((!defined($inputFile)) ){
   pod2usage( { -message => 'at least 1 parameter is mandatory: -i',
                 -verbose => 0,
                 -exitval => 2 } );
}

my $ostream     = IO::File->new();

# Manage input fasta file
my $ref_in = Bio::Tools::GFF->new(-file => $inputFile, -gff_version => 3);

# Manage Output
my $gffout;
if ($outfile) {
  $outfile=~ s/.gff//g;
open(my $fh, '>', $outfile.".gff") or die "Could not open file '$outfile' $!";
  $gffout= Bio::Tools::GFF->new(-fh => $fh, -gff_version => 3 );
}
else{
  $gffout = Bio::Tools::GFF->new(-fh => \*STDOUT, -gff_version => 3);
}


#time to calcul progression
my $startP=time;
my $nbLine=`wc -l < $inputFile`;
$nbLine =~ s/ //g;
chomp $nbLine;
print "$nbLine line to process...\n";
my $line_cpt=0;

my $count=0;
my %check; # keep track of signature seen

while (my $feature = $ref_in->next_feature() ) {
  $line_cpt++;

  my $position=lc($feature->seq_id)."".lc($feature->primary_tag)."".$feature->start()."".$feature->end(); #uniq position

  if(exists ($check{$position} ) ){
    $count++;
    next;
  }
  else{
    $gffout->write_feature($feature);
    $check{$position}++;
  }

  #Display progression
  if ((30 - (time - $startP)) < 0) {
    my $done = ($line_cpt*100)/$nbLine;
    $done = sprintf ('%.0f', $done);
        print "Progression : $done % processed.\n";
    $startP= time;
  }
}

if($count > 0){
  print "$count entries removed !\n";
}
else{print "No entry removed !\n";}
my $end_run = time();
my $run_time = $end_run - $start_run;
print "Job done in $run_time seconds\n";

__END__

=head1 NAME

gff3_remove_redundant_entries.pl -
remove redundant entries: same seq_id,primary_tag,start,stop.

=head1 SYNOPSIS

    gff3_remove_redundant_entries.pl -i <input file> [-o <output file>]
    gff3_remove_redundant_entries.pl --help

=head1 OPTIONS

=over 8

=item B<-i>, B<--gff>, B<--file> or B<--input>

STRING: Input gff file that will be read.


=item B<-o> or B<--output>

STRING: Output file.  If no output file is specified, the output will be written to STDOUT. The result is in tabulate format.

=item B<--help> or B<-h>

Display this helpful text.

=back

=cut
