#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use warnings;
use Pod::Usage;
use Getopt::Long;
use IO::File ;
use Bio::SeqIO;
use AGAT::AGAT;

start_script();
my $header = get_agat_header();
# ---------------------------- OPTIONS ----------------------------
my $start_run = time();
my @inputFile;
my $outputFile;
my $genome;
my $inflate;
my $opt_help = 0;

# OPTION MANAGEMENT: partition @ARGV into shared vs script options via library
my ($shared_argv, $script_argv) = split_argv_shared_vs_script(\@ARGV);

# Parse script-specific options from its own list
my $script_parser = Getopt::Long::Parser->new;
$script_parser->configure('bundling','no_auto_abbrev');
if ( ! $script_parser->getoptionsfromarray(
  $script_argv,
  'i|file|input|gff=s' => \@inputFile,
  'o|out|output=s'     => \$outputFile,
  'inflate!'        => \$inflate,
  'g|genome=s'      => \$genome,
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

if (! @inputFile ){
   pod2usage( {  -message => "$header\nAt least 1 input file is mandatory",
                 -verbose => 0,
                 -exitval => 1 } );
}

# Parse shared options (CPU, config, etc.)
my ($shared_opts) = parse_shared_options($shared_argv);

# --- Manage config ---
initialize_agat({ config_file_in => ( $shared_opts->{config} ), input => $inputFile[0], shared_opts => $shared_opts });

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
# Manage input gff file
  dual_print1 "Reading $file\n";
	my $format = $CONFIG->{force_gff_input_version};
	if(! $format ){ $format = select_gff_format($file); }
  my $inputfh = open_maybe_gz($file);
  my $ref_in = AGAT::BioperlGFF->new(-fh => $inputfh, -gff_version => $format);

  # set counter for progression bar
  set_progression_counter( $file);
  my $line_cpt=0;

  # parse gff
  while (my $feature = $ref_in->next_feature() ) {
    $line_cpt++;
    my $type = lc($feature->primary_tag);

		my $nb_parent = 1;

		# count number of parent if option activated
		if ($inflate){
			if($feature->has_tag('Parent')){
				my @parentList = $feature->get_tag_values('Parent');
				$nb_parent = scalar @parentList;
			}
		}

		# count several time is several parents
		for (my $i=0; $i<$nb_parent; $i++){
    	$type_count->{$type}++;
    	$type_bp->{$type} += ($feature->end()-$feature->start())+1;
		}

    #Display progression
    update_progression_counter($line_cpt);
  }
  dual_print1 "\rProgress : 100 %\n";
end_script();
}

my $totalNumber=0;
my $totalSize=0;

if(defined($genomeSize)){
print $ostream "Type (3rd column)\tNumber\tSize total (kb)\tSize mean (bp)\t% of the genome\t/!\\Results are rounding to two decimal places \n";
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
  print $ostream "Type (3rd column)\tNumber\tSize total (kb)\tSize mean (bp)\t/!\\Results are rounding to two decimal places \n";
  foreach my $gnx (sort {$a cmp $b} keys(%$type_count)) {
    my $Sitotal=sprintf("%0.2f",($type_bp->{$gnx}/1000));
    my $SizeMean=sprintf("%0.2f",($type_bp->{$gnx}/$type_count->{$gnx}));
    print $ostream $gnx,"\t",$type_count->{$gnx},"\t",$Sitotal,"\t",$SizeMean,"\n";

    $totalNumber += $type_count->{$gnx};
    $totalSize += $type_bp->{$gnx};

  }
}

my $goodTotalSize=sprintf("%0.2f",($totalSize/1000));
my $goodTotalSizeMean=sprintf("%0.2f",($totalSize/$totalNumber));

if(defined($genomeSize)){
  my $goodxGenome=sprintf("%0.2f",($totalSize/$genomeSize)*100);
  print $ostream "Total\t",$totalNumber,"\t",$goodTotalSize,"\t",$goodTotalSizeMean,"\t",$goodxGenome,"\n";
}
else{
  print $ostream "Total\t",$totalNumber,"\t",$goodTotalSize,"\t",$goodTotalSizeMean,"\n";
}

# --- final messages ---
end_script();

__END__

=head1 NAME

agat_sq_stat_basic.pl

=head1 DESCRIPTION

The script aims to provide basic statistics of a gtf/gff file.

=head1 SYNOPSIS

    agat_sq_stat_basic.pl -i <input file> [-g <integer or fasta> -o <output file>]
    agat_sq_stat_basic.pl --help

=head1 OPTIONS

=over 8

=item B<-i>, B<--gff>, B<--file> or B<--input> <file>

STRING: Input GTF/GFF file. Several files can be processed at once: -i file1 -i file2

=item B<-g>, B<--genome> <integer or fasta>

That input is design to know the genome size in order to calculate the percentage of the genome represented by each kind of feature type.
You can provide an INTEGER or the genome in fasta format. If you provide the fasta, the genome size will be calculated on the fly.

=item B<--inflate>

Inflate the statistics taking into account feature with multi-parents.
Indeed to avoid redundant information, some gff factorize identical features.
e.g: one exon used in two different isoform will be defined only once, and will have multiple parent.
By default the script count such feature only once. Using the inflate option allows
to count the feature and its size as many time there are parents.

=item  B<-o>, B<--out> or B<--output> <file>

Output file to create (default GFF3 - see config to modify output format).
If no output file is specified, the output will be written to STDOUT.

=item B<--help> or B<-h>

Display this helpful text.

=back

=head1 SHARED OPTIONS

Shared options are defined in the AGAT configuration file and can be overridden via the command line for this script only.
Common shared options are listed below; for the full list, please refer to the AGAT agat_config.yaml.
Note: For _sq_ scripts, only the following options are supported: verbose, output_format, gff_output_version, gtf_output_version, progress_bar, and tabix.

=over 8

=item B<--config> <file>

Path to a custom AGAT configuration file.  
By default, AGAT uses `agat_config.yaml` from the working directory if present, otherwise the default file shipped with AGAT
(available locally via `agat config --expose`).

=item B<-v> or B<--verbose> <int>

Verbosity, choice are 0,1,2,3,4. 0 is quiet, 1 is normal, 2,3,4 is more verbose. Default 1.

=back

=head1 FEEDBACK

For questions, suggestions, or general discussions about AGAT, please use the AGAT community forum:
https://github.com/NBISweden/AGAT/discussions

=head1 BUG REPORTING

Bug reports should be submitted through the AGAT GitHub issue tracker:
https://github.com/NBISweden/AGAT/issues

=cut

AUTHOR - Jacques Dainat
