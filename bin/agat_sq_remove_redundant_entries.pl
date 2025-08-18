#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Pod::Usage;
use Getopt::Long;
use IO::File ;
use AGAT::AGAT;

my $header = get_agat_header();
my $config;
my $start_run = time();
my $inputFile;
my $verbose;
my $outfile;
my $opt_help = 0;

my $common = parse_common_options() || {};
$config   = $common->{config};
$outfile  = $common->{output};
$verbose  = $common->{verbose};
$opt_help = $common->{help};

if ( !GetOptions ('i|file|input|gff=s' => \$inputFile,
                    )  )
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
                 -exitval => 2 } );
}

# --- Manage config ---
$config = get_agat_config({config_file_in => $config});

# Manage input gff file
my $format = $config->{force_gff_input_version};
if(! $format ){ $format = select_gff_format($inputFile); }
my $ref_in = AGAT::BioperlGFF->new(-file => $inputFile, -gff_version => $format);

# Manage Output
my $gffout = prepare_gffout($config, $outfile);

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

  my $parent="";
  if($feature->has_tag('Parent')){
    $parent = $feature->_tag_value('Parent');
  }
  my $id="";
  if($feature->has_tag('ID')){
    $id = $feature->_tag_value('ID');
  }
  if ($parent eq "" and $id eq "" ){next;}

  my $position=lc($feature->seq_id)."".lc($feature->primary_tag)."".$feature->start()."".$feature->end()."".$id."".$parent; #uniq position

  if(exists ($check{$position} ) ){
    $count++;
    if ($verbose){print "remove: ".$feature->gff_string."\n";}
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

agat_sq_remove_redundant_entries.pl

=head1 DESCRIPTION

The script remove redundant entries: same seq_id,primary_tag,start,stop,ID,Parent.
If ID and Parent attribute is not present, we do no remove the feature. If one of them
do not exists we use "" instead.

=head1 SYNOPSIS

    agat_sq_remove_redundant_entries.pl -i <input file> [-o <output file>]
    agat_sq_remove_redundant_entries.pl --help

=head1 OPTIONS

=over 8

=item B<-i>, B<--gff>, B<--file> or B<--input>

STRING: Input GTF/GFF file.

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
