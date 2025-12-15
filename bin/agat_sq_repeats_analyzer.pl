#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Pod::Usage;
use Getopt::Long;
use IO::File ;
use Bio::SeqIO;
use AGAT::AGAT;

start_script();
my $header = get_agat_header();
# ---------------------------- OPTIONS ----------------------------
my @inputFile;
my $outputFile;
my $genome;
my $opt_help = 0;

# OPTION MANAGEMENT: split shared vs script-specific argv
my ($shared_argv, $script_argv) = split_argv_shared_vs_script(\@ARGV);

# Parse script options
my $script_parser = Getopt::Long::Parser->new;
$script_parser->configure('bundling','no_auto_abbrev');
if ( ! $script_parser->getoptionsfromarray(
  $script_argv,
  'i|file|input|gff=s' => \@inputFile,
  'o|output=s'         => \$outputFile,
  'g|genome=s'         => \$genome,
  'h|help!'            => \$opt_help )  )
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

if (! @inputFile ){
   pod2usage( { -message => "$header\nAt least 1 input file is mandatory",
                 -verbose => 0,
                 -exitval => 1 } );
}

# Parse shared options
my ($shared_opts) = parse_shared_options($shared_argv);

# --- Manage config ---
initialize_agat({ config_file_in => $shared_opts->{config}, input => $inputFile[0], shared_opts => $shared_opts });

# ----------------------------

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
  dual_print1 sprintf("%-45s%d%s", "Total sequence length", $genomeSize,"\n");
  }

#time to calcul progression
my $type_count;
my $type_bp;
my %check; #track the repeat already annotated to not. Allow to skip already read repeats

foreach my $file (@inputFile){
# Manage input fasta file
  dual_print1 "Reading $file\n";
	my $format = $CONFIG->{force_gff_input_version};
	if(! $format ){ $format = select_gff_format($file); }
  my $ref_in = AGAT::BioperlGFF->new(-file => $file, -gff_version => $format);

  my $startP=time;
  my $nbLine=`wc -l < $file`;
  $nbLine =~ s/ //g;
  chomp $nbLine;
  dual_print1 "$nbLine line to process...\n";
  my $line_cpt=0;

  local $| = 1; # Or use IO::Handle; STDOUT->autoflush; Use to print progression bar
  while (my $feature = $ref_in->next_feature() ) {
    $line_cpt++;
    my $type = lc($feature->primary_tag);
    ## repeatMasker or repeatRunner
    if (($type eq 'match') or ($type eq 'protein_match')){

      my $position=$feature->seq_id."".$feature->start()."".$feature->end(); #uniq position
      if(exists ($check{$position} ) ){next;}
      else{

       my $nameAtt=$feature->_tag_value('Name');
       my $genus=(split ":", (split /\|/, (split /\s+/,$nameAtt)[0])[-1])[-1];
       $type_count->{$genus}++;
       $type_bp->{$genus} += ($feature->end()-$feature->start())+1;
       $check{$position}++;
      }
    }

    #Display progression
    if ((30 - (time - $startP)) < 0) {
      my $done = ($line_cpt*100)/$nbLine;
      $done = sprintf ('%.0f', $done);
          dual_print1 "\rProgress : $done %";
      $startP= time;
    }
  }
  dual_print1 "\rProgress : 100 %\n";
}

my $totalNumber=0;
my $totalSize=0;

if(defined($genomeSize)){
print $ostream "Repeat type\tNumber\tSize total (kb)\tSize mean (bp)\t% of the genome\t/!\\Results are rounding to two decimal places \n";
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
  print $ostream "Repeat type\tNumber\tSize total (kb)\tSize mean (bp)\t/!\\Results are rounding to two decimal places \n";
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

# --- final messages ---
end_script();

__END__

=head1 NAME

agat_sq_repeats_analyzer.pl

=head1 DESCRIPTION

The script allows to generate a tabulated format report of repeats annotated
from a gff file containing repeats (feature type must be match or protein_match).

=head1 SYNOPSIS

    agat_sq_repeats_analyzer.pl -i <input file> [-g <integer or fasta> -o <output file>]
    agat_sq_repeats_analyzer.pl --help

=head1 OPTIONS

=over 8

=item B<-i>, B<--gff>, B<--file> or B<--input>

STRING: Input GTF/GFF file(s). Several files can be processed at once: -i file1 -i file2

=item B<-g>, B<--genome>

That input is design to know the genome size in order to calculate the percentage of the genome represented by each kind of repeats.
You can provide an INTEGER or the genome in fasta format. If you provide the fasta, the genome size will be calculated on the fly.

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
