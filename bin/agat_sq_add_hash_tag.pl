#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Pod::Usage;
use Getopt::Long;
use IO::File ;
use AGAT::AGAT;

my $header = get_agat_header();
my $config = get_agat_config();
my $start_run = time();
my $inputFile=undef;
my $outfile=undef;
my $opt_help = 0;
my $interval=1;

Getopt::Long::Configure ('bundling');
if ( !GetOptions ('file|input|gff=s' => \$inputFile,
      'i|interval=i' => \$interval,
      'o|output=s' => \$outfile,
      'h|help!'         => \$opt_help )  )
{
    pod2usage( { -message => 'Failed to parse command line',
                 -verbose => 1,
                 -exitval => 1 } );
}

# Print Help and exit
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

if (( $interval > 2 or $interval < 1) ){
   pod2usage( { -message => 'interval must be 1 or 2. Have a look to the help to know more',
                 -verbose => 1,
                 -exitval => 1 } );
}

# Manage input gff file
my $format = $config->{gff_output_version};
if(! $format ){ $format = select_gff_format($inputFile); }
my $ref_in = Bio::Tools::GFF->new(-file => $inputFile, -gff_version => $format);

# Manage Output
my $gffout = prepare_gffout($config, $outfile);
my $gffXtra=$gffout->{"_filehandle"}; #to add extra lines to gff!!

#time to calcul progression
my $startP=time;
my $nbLine=`wc -l < $inputFile`;
$nbLine =~ s/ //g;
chomp $nbLine;
print "$nbLine line to process...\n";

my $line_cpt=0;
my $count=0;
my $nextGroup=0;
my @bucket=();
my $before="";
my $actual="";
while (my $feature = $ref_in->next_feature() ) {
  $line_cpt++;

  #What do we follow

  if ($interval eq "1"){ #per sequence

    $actual=lc($feature->seq_id);
    if( ($actual ne $before) and ($actual ne "" and $before ne "") ) {
      _write_bucket(\@bucket, $gffout);
      $count++;
      $nextGroup=0;
      @bucket=();
    }
    push (@bucket,$feature);
    $before=lc($feature->seq_id);
  }


  if ($interval eq "2"){ #per feature group
    $actual=lc($feature->primary_tag);
    if ( ($actual ne $before) and ($before ne "") and ($actual eq "gene" or $actual eq "expressed_sequence_match" or $actual eq "match") ) {
      _write_bucket(\@bucket, $gffout);
      $count++;
      $nextGroup=0;
      @bucket=();
    }
    push (@bucket,$feature);
    $before=lc($feature->primary_tag);
  }

  #Display progression
  if ((30 - (time - $startP)) < 0) {
    my $done = ($line_cpt*100)/$nbLine;
    $done = sprintf ('%.0f', $done);
        print "\rProgression : $done % processed.\n";
    $startP= time;
  }
}

##Last round
 _write_bucket(\@bucket, $gffout);
$count++;

if($count > 0){
  print "$count line added !\n";
}
else{print "No line added !\n";}
my $end_run = time();
my $run_time = $end_run - $start_run;
print "Job done in $run_time seconds\n";


sub _write_bucket{
  my($bucket, $gffout)=@_;
  foreach my $feature (@$bucket){
    $gffout->write_feature($feature);
  }

  # Get the filehandle
  print $gffXtra "###\n";
}

__END__

=head1 NAME

agat_sq_add_hash_tag.pl

=head1 DESCRIPTION

The script aims to introduce hash tag (###) into the file. It allows for some tools
using gff3 to handle independantly file chucks separated by the ### signal. Can make
them more efficient.

=head1 SYNOPSIS

    agat_sq_add_hash_tag.pl -i <input file> [-o <output file>]
    agat_sq_add_hash_tag.pl --help

=head1 OPTIONS

=over 8

=item B<--gff>, B<--file> or B<--input>

STRING: Input GTF/GFF file.

=item B<-i> or B<--interval>

Integer: 1 or 2. 1 will add ### after each new sequence (column1 of the gff), while 2 will add the ### after each group of feature (gene).
By default the value is 1.

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
