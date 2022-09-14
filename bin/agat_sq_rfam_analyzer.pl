#!/usr/bin/env perl

use Carp;
use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use IO::File ;
use Bio::SeqIO;
use Bio::Tools::GFF;
use AGAT::Omniscient;

my $header = get_agat_header();
my $config = get_agat_config();
my $start_run = time();
my @inputFile;
my $outputFile;
my $genome;
my $opt_help = 0;

Getopt::Long::Configure ('bundling');
if ( !GetOptions ('i|file|input|gff=s' => \@inputFile,
      'o|output=s' => \$outputFile,
      'g|genome=s' => \$genome,
      'h|help!'         => \$opt_help )  )
{
    pod2usage( { -message => "Failed to parse command line",
                 -verbose => 1,
                 -exitval => 1 } );
}

if ($opt_help) {
    pod2usage( { -verbose => 99,
                 -exitval => 0,
                 -message => "$header \n" } );
}

if (! @inputFile ){
   pod2usage( { -message => "$header\nAt least 1 input file is mandatory",
                 -verbose => 0,
                 -exitval => 1 } );
}

# Manage Output
my $ostream = prepare_fileout($outputFile);

#check genome size
my $genomeSize=undef;
  if($genome){
    if( $genome =~ /^[0-9]+$/){ #check if it's a number
      $genomeSize=$genome;
    }
    elsif($genome){
      my $seqio = Bio::SeqIO->new(-file => $genome, '-format' => 'Fasta');
      while(my $seq = $seqio->next_seq) {
          my $string = $seq->seq;
          $genomeSize += length($string);
        }
    }
  printf("%-45s%d%s", "Total sequence length", $genomeSize,"\n");
  }

#time to calcul progression
my $type_count;
my $type_bp;
my %check; #track the repeat already annotated to not. Allow to skip already read repeats

foreach my $file (@inputFile){
# Manage input fasta file
  print "Reading $file\n";
	my $format = $config->{gff_output_version};
	if(! $format ){ $format = select_gff_format($file); }
  my $ref_in = Bio::Tools::GFF->new(-file => $file, -gff_version => $format);

  my $startP=time;
  my $nbLine=`wc -l < $file`;
  $nbLine =~ s/ //g;
  chomp $nbLine;
  print "$nbLine line to process...\n";
  my $line_cpt=0;

  local $| = 1; # Or use IO::Handle; STDOUT->autoflush; Use to print progression bar
  while (my $feature = $ref_in->next_feature() ) {
    $line_cpt++;
    my $type = lc($feature->primary_tag);
    ## repeatMasker or repeatRunner
    if (($type eq 'ncrna') or ($type eq 'nc_rna')){

      my $ID = lc($feature->_tag_value('ID'));
      my $position=$feature->seq_id."".$feature->start()."".$feature->end(); #uniq position
      my $genus=$feature->_tag_value('rfam-id');

      if(exists ($check{$ID}) ) {
        if(! exists($check{$ID}{$position} ) ){
          $type_bp->{$genus} += ($feature->end()-$feature->start())+1;
          $check{$ID}{$position}++;
        }
      }
      else{
       $type_count->{$genus}++;
       $type_bp->{$genus} += ($feature->end()-$feature->start())+1;
       $check{$ID}{$position}++;
      }
    }

    #Display progression
    if ((30 - (time - $startP)) < 0) {
      my $done = ($line_cpt*100)/$nbLine;
      $done = sprintf ('%.0f', $done);
          print "\rProgress : $done %";
      $startP= time;
    }
  }
  print "\rProgress : 100 %\n";
}

my $totalNumber=0;
my $totalSize=0;

if(defined($genomeSize)){
print $ostream "ncRNA type\tNumber\tSize total (kb)\tSize mean (bp)\t% of the genome\t/!\\Results are rounding to two decimal places \n";
  foreach my $gnx (sort {$a cmp $b} keys(%$type_count)) {
    my $Sitotal=sprintf("%0.2f",($type_bp->{$gnx}/1000));
    my $SizeMean=sprintf("%0.2f",($type_bp->{$gnx}/$type_count->{$gnx}));
    my $xGenome=sprintf("%0.2f",($type_bp->{$gnx}/$genomeSize)*100);
    print $ostream $gnx,"\t",$type_count->{$gnx},"\t",$Sitotal,"\t",$SizeMean,"\t",$xGenome,"\n";

    $totalNumber += $type_count->{$gnx};
    $totalSize += $type_bp->{$gnx};

  }
}
else{
  print $ostream "ncRNA type\tNumber\tSize total (kb)\tSize mean (bp)\t/!\\Results are rounding to two decimal places \n";
  foreach my $gnx (sort {$a cmp $b} keys(%$type_count)) {
    my $Sitotal=sprintf("%0.2f",($type_bp->{$gnx}/1000));
    my $SizeMean=sprintf("%0.2f",($type_bp->{$gnx}/$type_count->{$gnx}));
    print $ostream $gnx,"\t",$type_count->{$gnx},"\t",$Sitotal,"\t",$SizeMean,"\n";

    $totalNumber += $type_count->{$gnx};
    $totalSize += $type_bp->{$gnx};

  }
}

if($totalSize){
  my $goodTotalSize=sprintf("%0.2f",($totalSize/1000));
  my $goodTotalSizeMean=sprintf("%0.2f",($totalSize/$totalNumber));

  if(defined($genomeSize)){
    my $goodxGenome=sprintf("%0.2f",($totalSize/$genomeSize)*100);
    print $ostream "Total\t",$totalNumber,"\t",$goodTotalSize,"\t",$goodTotalSizeMean,"\t",$goodxGenome,"\n";
  }
  else{
    print $ostream "Total\t",$totalNumber,"\t",$goodTotalSize,"\t",$goodTotalSizeMean,"\n";
  }
}
else{
  print $ostream "None found\n";
}

  my $end_run = time();
  my $run_time = $end_run - $start_run;
  print "Job done in $run_time seconds\n";

__END__

=head1 NAME

agat_sq_rfam_analyzer.pl

=head1 DESCRIPTION

The script allows to generate a tabulated format report of rfam-id annotated from a gff file
containing rfam results (type of the 3rd column must be ncRNA or nc_RNA - not case sensitive. And the 9th column must contain the rfam-id attribute).
    e.g:
ScG6Pog_82  Rfam  ncRNA 737595  737663  20.7  + 0 ID=RF00134_ScG6Pog_82_737595;Name=RF00134_ScG6Pog_82_737595;evalue=0.45;gc-content=0.28;model_end=1;model_start=1;rfam-acc=RF00134;rfam-id=snoZ196
ScG6Pog_82  Rfam  ncRNA 305023  305103  20.8  + 0 ID=RF00227_ScG6Pog_82_305023;Name=RF00227_ScG6Pog_82_305023;evalue=0.35;gc-content=0.31;model_end=1;model_start=1;rfam-acc=RF00227;rfam-id=FIE3

=head1 SYNOPSIS

    agat_sq_rfam_analyzer.pl -i <input file> [-g <integer or fasta> -o <output file>]
    agat_sq_rfam_analyzer.pl --help

=head1 OPTIONS

=over 8

=item B<-i>, B<--gff>, B<--file> or B<--input>

STRING: Input GTF/GFF file(s). Several files can be processed at once: -i file1 -i file2

=item B<-g>, B<--genome>

That input is design to know the genome size in order to calculate the percentage of the genome represented by each kind of rfam-id.
You can provide an INTEGER or the genome in fasta format. If you provide the fasta, the genome size will be calculated on the fly.

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
