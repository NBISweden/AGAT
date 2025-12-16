#!/usr/bin/env perl

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use Bio::SeqIO ;
use IO::File ;
use AGAT::AGAT;

start_script();
my $header = get_agat_header();
# ---------------------------- OPTIONS ----------------------------
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
my ($shared_argv, $script_argv) = split_argv_shared_vs_script(\@ARGV);
my $script_parser = Getopt::Long::Parser->new;
$script_parser->configure('bundling','no_auto_abbrev');
if ( !$script_parser->getoptionsfromarray( $script_argv,
          'g|gff=s'         => \$opt_gfffile,
          'f|fa|fasta=s'    => \$opt_fastafile,
          'hm:s'            => \$opt_HardMask,
          'sm'              => \$opt_SoftMask,
          'o|output=s'      => \$opt_output,
          'h|help!'         => \$opt_help ) )
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

my ($shared_opts) = parse_shared_options($shared_argv);
initialize_agat({ config_file_in => $shared_opts->{config}, input => $opt_gfffile, shared_opts => $shared_opts });

# ------------------------------------------------------------------------------

if (defined ($opt_HardMask) && defined ($opt_SoftMask)){
  dual_print1 "It is not possible to HardMask and SoftMask at the same time. Choose only one option and try again !\n"; exit();
}

my $ostream = prepare_fileout($opt_output);

if (defined( $opt_HardMask)){
  dual_print1 "You choose to Hard Mask the genome.\n";
	if (! $opt_HardMask){
	  $hardMaskChar = "n";
	}
	elsif(length($opt_HardMask) == 1){
	  $hardMaskChar = $opt_HardMask;
	}
  else{dual_print1 "$opt_HardMask cannot be used to Mask. A character is mandatory.\n";exit;}
  dual_print1 "Character used for Mask: $hardMaskChar\n";
}
if (defined( $opt_SoftMask)){
  dual_print1 "You choose to Soft Mask the genome.\n";
}
##### MAIN ####

#### read gff file and save info in memory
my %gff; my $nbLineRead=0;

# Manage input gff file
my $format = $CONFIG->{force_gff_input_version};
if(! $format ){ $format = select_gff_format($opt_gfffile); }
my $gff_in = AGAT::BioperlGFF->new(-file => $opt_gfffile, -gff_version => $format);


dual_print1 "Reading features from $opt_gfffile...\n";
  while (my $feature = $gff_in->next_feature()) {
    my $seqname=$feature->seq_id();
    my $start=$feature->start();
    my $end=$feature->end();
   	push @{$gff{uc $seqname}},"$start $end";
    $nbLineRead++;
   }
$gff_in->close();
dual_print1 "$nbLineRead lines read\n";

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
dual_print1 "$nbFastaSeq fasta sequences read.\n";
dual_print1 "$nucl_masked nucleotides masked.\n";

# --- final messages ---
end_script();

# ---------------------------- FUNCTIONS ----------------------------
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

=item B<-h> or B<--help>

Display this helpful text.

=back

=head1 SHARED OPTIONS

Shared options are defined in the AGAT configuration file and can be overridden via the command line for this script only.
Common shared options are listed below; for the full list, please refer to the AGAT agat_config.yaml.
Note: For _sq_ scripts, only the following options are supported: verbose, output_format, gff_output_version, gtf_output_version, progress_bar, and tabix.

=over 8

=item B<--config>

String - Path to a custom AGAT configuration file.  
By default, AGAT uses `agat_config.yaml` from the working directory if present, otherwise the default file shipped with AGAT
(available locally via `agat config --expose`).

=item B<-v> or B<--verbose>

Integer - Verbosity, choice are 0,1,2,3,4. 0 is quiet, 1 is normal, 2,3,4 is more verbose. Default 1.

=back

=head1 FEEDBACK

For questions, suggestions, or general discussions about AGAT, please use the AGAT community forum:
https://github.com/NBISweden/AGAT/discussions

=head1 BUG REPORTING

Bug reports should be submitted through the AGAT GitHub issue tracker:
https://github.com/NBISweden/AGAT/issues

=cut

AUTHOR - Jacques Dainat
