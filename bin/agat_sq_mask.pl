#!/usr/bin/env perl

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use Bio::SeqIO ;
use IO::File ;
use AGAT::AGAT;

my $header = get_agat_header();
my $config;
my $start_run = time();
my $opt_HardMask;
my $opt_SoftMask;
my $opt_gfffile;
my $opt_fastafile;
my $opt_output;
my $opt_help = 0;

# Character for hardMask
my $hardMaskChar;
my $width = 60; # line length printed

# OPTION MANAGMENT
my $common = parse_common_options() || {};
$config     = $common->{config};
$opt_output = $common->{output};
$opt_help   = $common->{help};

if ( !GetOptions( 'g|gff=s'         => \$opt_gfffile,
                  'f|fa|fasta=s'    => \$opt_fastafile,
                  'hm:s'            => \$opt_HardMask,
                  'sm'              => \$opt_SoftMask,
                  ) )
{ 
    pod2usage( { -message => 'Failed to parse command line',
                 -verbose => 1,
                 -exitval => 1 } );
}

if ($opt_help) {
    pod2usage( {  -verbose => 99,
                  -exitval => 0,
                  -message => "$header\n" } );
}

if ( (! (defined($opt_gfffile)) ) || (! (defined($opt_fastafile)) ) || ( (! defined($opt_HardMask) && (! defined($opt_SoftMask))) ) ){
    pod2usage( {
           -message => "$header\nAt least 3 parametes are mandatory:\nInput reference gff file (-g);  Input reference fasta file (-f); Mask type (--hm for hard mask or --sm for soft mask)\n\n".
           "Ouptut is optional. Look at the help documentation to know more.\n",
           -verbose => 0,
           -exitval => 1 } );
}

# --- Manage config ---
$config = get_agat_config({config_file_in => $config});

if (defined ($opt_HardMask) && defined ($opt_SoftMask)){
  print "It is not possible to HardMask and SoftMask at the same time. Choose only one the options and try again !\n"; exit();
}

my $ostream = prepare_fileout($opt_output);

if (defined( $opt_HardMask)){
  print "You choose to Hard Mask the genome.\n";
	if (! $opt_HardMask){
	  $hardMaskChar = "n";
	}
	elsif(length($opt_HardMask) == 1){
	  $hardMaskChar = $opt_HardMask;
	}
	else{print "$opt_HardMask cannot be used to Mask. A character is mandatory.\n";exit;}
	print "Charcater uses for Mask: $hardMaskChar\n";
}
if (defined( $opt_HardMask)){
  print "You choose to Soft Mask the genome.\n";
}
##### MAIN ####

#### read gff file and save info in memory
my %gff; my $nbLineRead=0;

# Manage input gff file
my $format = $config->{force_gff_input_version};
if(! $format ){ $format = select_gff_format($opt_gfffile); }
my $gff_in = AGAT::BioperlGFF->new(-file => $opt_gfffile, -gff_version => $format);


print( "Reading features from $opt_gfffile...\n");
  while (my $feature = $gff_in->next_feature()) {
    my $seqname=$feature->seq_id();
    my $start=$feature->start();
    my $end=$feature->end();
   	push @{$gff{uc $seqname}},"$start $end";
    $nbLineRead++;
   }
$gff_in->close();
print "$nbLineRead lines read\n";

#### read fasta
my $nbFastaSeq=0;
my $nucl_masked=0;
my $inFasta  = Bio::SeqIO->new(-file => "$opt_fastafile" , '-format' => 'Fasta');

while ($_=$inFasta->next_seq()) {
    my $seqname = $_->id;
    my $sequence = $_->seq;

    foreach (@{$gff{uc $seqname}}) {
	    my ($start,$end) = split;
      if ($opt_SoftMask){
        my $strinTolo = substr($sequence,$start-1,$end+1-$start);
        substr($sequence,$start-1,$end+1-$start) = lc $strinTolo;
      }
      else{
	     substr($sequence,$start-1,$end+1-$start) = $hardMaskChar x ($end+1-$start);
      }
      $nucl_masked += ($end-$start+1);
    }

    print $ostream ">$seqname\n";
    for (my $i=0; $i<length $sequence; $i += $width) { print $ostream substr($sequence,$i,$width)."\n" }
    $nbFastaSeq++;
}
$inFasta->close();
print "$nbFastaSeq fasta sequences read.\n";
print "$nucl_masked nucleotides masked.\n";
my $end_run = time();
my $run_time = $end_run - $start_run;
print "Job done in $run_time seconds\n";
__END__

=head1 NAME

agat_sq_mask.pl

=head1 DESCRIPTION

The script masks (hard or soft) GFF-denoted segments out of a FASTA format file.
It needs 3 input parameters: a gff3 file, a fasta file, and a Mask method.
The result is written to the specified output file, or to STDOUT.

=head1 SYNOPSIS

    agat_sq_mask.pl -g infile.gff -f infile.fasta  (-hm or -sm) [ -o outfile ]
    agat_sq_mask.pl --help

=head1 OPTIONS

=over 8

=item B<-g>, B<--gff> or B<-ref>

Input GTF/GFF file.

=item B<-f>, B<--fa> or B<--fasta>

Input fasta file that will be masked

=item B<-sm>

SoftMask option =>Sequences masked will be in lowercase

=item B<-hm>

HardMask option => Sequences masked will be replaced by a character. By default the character used is 'n'. But you are allowed to speceify any character of your choice. To use 'z' instead of 'n' type: -hm z

=item B<-o> or B<--output>

Output GFF file.  If no output file is specified, the output will be
written to STDOUT.

=item B<-c> or B<--config>

String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any, 
otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

=item B<-h> or B<--help>

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
