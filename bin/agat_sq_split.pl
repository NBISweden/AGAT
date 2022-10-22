#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Pod::Usage;
use Getopt::Long;
use File::Basename;
use IO::File ;
use AGAT::AGAT;

my $header = get_agat_header();
my $config = get_agat_config();
my $start_run = time();
my $inputFile=undef;
my $outfolder=undef;
my $opt_help = 0;
my $interval=10;
my $feature_type="gene";

Getopt::Long::Configure ('bundling');
if ( !GetOptions ('file|input|gff=s' => \$inputFile,
      'ft|feature_type=s' => \$feature_type,
      'i|interval=i' => \$interval,
      'o|output=s' => \$outfolder,
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

if ( !(defined($inputFile)) or !(defined($outfolder)) ){
   pod2usage( { -message => "$header\nAt least 2 parameters are mandatory: -i inputFile and -o $outfolder",
                 -verbose => 0,
                 -exitval => 1 } );
}

# Manage input fasta file
my $format = $config->{gff_output_version};
if(! $format ){ $format = select_gff_format($inputFile); }
my $ref_in = Bio::Tools::GFF->new(-file => $inputFile, -gff_version => $format);

# Manage Output
if (-d $outfolder) {
  print "The output directory <$outfolder> already exists.\n";exit;
}
else{
  my ($path,$ext);
  ($outfolder,$path,$ext) = fileparse($outfolder,qr/\.[^.]*/);
  print "Creating the $outfolder folder\n";
  mkdir $outfolder;
}

print "I will split the file into files containing $interval group of feature. The top feature of the group of feature is currenlty defined by <$feature_type>.\n";

#time to calcul progression
my $startP=time;
my $nbLine=`wc -l < $inputFile`;
$nbLine =~ s/ //g;
chomp $nbLine;
print "$nbLine line to process...\n";
my $line_cpt=0;

my $count_feature=0;
my $count_file=1;
my ($file_name,$path,$ext) = fileparse($inputFile,qr/\.[^.]*/);

my $gffout = prepare_gffout($config, $outfolder."/".$file_name."_".$count_file.".gff");

while (my $feature = $ref_in->next_feature() ) {
  $line_cpt++;

  #What do we follow
  if($feature->primary_tag eq $feature_type){
    if($count_feature == $interval){
      close $gffout;
      $count_file++;
			$gffout = prepare_gffout($config,  $outfolder."/".$file_name."_".$count_file.".gff");
      $count_feature=0;
    }
    $count_feature++;
  }
  $gffout->write_feature($feature);

  #Display progression
  if ((30 - (time - $startP)) < 0) {
    my $done = ($line_cpt*100)/$nbLine;
    $done = sprintf ('%.0f', $done);
        print "\rProgression : $done % processed.\n";
    $startP= time;
  }
}
close $gffout;

my $end_run = time();
my $run_time = $end_run - $start_run;
print "Job done in $run_time seconds\n";

__END__

=head1 NAME

agat_sq_split.pl

=head1 DESCRIPTION

split gff3 file into several files.
By default we create files containing 1000 genes and all sub-features associated.
GFF3 input file must be sequential.

=head1 SYNOPSIS

    agat_sq_split.pl -i <input file> -o <output file>
    agat_sq_split.pl --help

=head1 OPTIONS

=over 8

=item B<--gff>, B<--file> or B<--input>

STRING: Input GTF/GFF file.

=item B<-i> or B<--interval>
Integer.  Number of group of feature to include in each file. 1000 by default.

=item B<--ft> or B<--feature_type>
The top feature of the feature group. By default "gene".

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
