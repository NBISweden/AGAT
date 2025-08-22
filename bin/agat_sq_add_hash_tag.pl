#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Pod::Usage;
use Getopt::Long::Descriptive;
use IO::File;
use AGAT::AGAT;

my $header = get_agat_header();
my $start_run = time();
my @copyARGV = @ARGV;
my ( $opt, $usage, $config ) = AGAT::AGAT::describe_script_options( $header,
    [ 'file|input|gff=s', 'Input GFF file', { required => 1 } ],
    [ 'i|interval=i', 'Interval (1 or 2)',
      { default => 1,
        callbacks => { valid => sub { $_[0] == 1 || $_[0] == 2 or die 'interval must be 1 or 2. Have a look to the help to know more' } } } ],
);

my $inputFile = $opt->file;
my $interval  = $opt->i;

my $log;
if ( my $log_name = $config->{log_path} ) {
    open( $log, '>', $log_name ) or die "Can not open $log_name for printing: $!";
    dual_print( $log, $header, 0 );
}

# Manage input gff file
my $format = $config->{force_gff_input_version};
if(! $format ){ $format = select_gff_format($inputFile); }
my $ref_in = AGAT::BioperlGFF->new(-file => $inputFile, -gff_version => $format);

# Manage Output
my $outfile = $config->{output};
my $gffout = prepare_gffout($config, $outfile);
my $gffXtra=$gffout->{"_filehandle"}; #to add extra lines to gff!!

#time to calcul progression
my $startP=time;
my $nbLine=`wc -l < $inputFile`;
$nbLine =~ s/ //g;
chomp $nbLine;
dual_print($log, "$nbLine line to process...\n", $config->{verbose});

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
        dual_print($log, "\rProgression : $done % processed.\n", $config->{verbose});
    $startP= time;
  }
}

##Last round
 _write_bucket(\@bucket, $gffout);
$count++;

if($count > 0){
  dual_print($log, "$count line added !\n", $config->{verbose});
}
else{dual_print($log, "No line added !\n", $config->{verbose});}
my $end_run = time();
my $run_time = $end_run - $start_run;
dual_print($log, "Job done in $run_time seconds\n", $config->{verbose});

close $log if $log;


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

=item B<-c> or B<--config>

String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any, 
otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

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
