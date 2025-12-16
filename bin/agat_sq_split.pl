#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Pod::Usage;
use Getopt::Long;
use File::Basename;
use IO::File ;
use AGAT::AGAT;

start_script();
my $header = get_agat_header();
# ---------------------------- OPTIONS ----------------------------
my $start_run = time();
my $inputFile=undef;
my $outfolder=undef;
my $opt_help = 0;
my $interval=10;
my $feature_type="gene";

# OPTION MANAGEMENT: partition @ARGV into shared vs script options via library
my ($shared_argv, $script_argv) = split_argv_shared_vs_script(\@ARGV);

# Parse script-specific options from its own list
my $script_parser = Getopt::Long::Parser->new;
$script_parser->configure('bundling','no_auto_abbrev');
if ( ! $script_parser->getoptionsfromarray(
  $script_argv,
  'file|input|gff=s' => \$inputFile,
  'ft|feature_type=s'        => \$feature_type,
  'i|interval=i'             => \$interval,
  'o|output=s'               => \$outfolder,
  'h|help!'                  => \$opt_help )  )
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

# Parse shared options (CPU, config, etc.)
my ($shared_opts) = parse_shared_options($shared_argv);

# --- Manage config ---
initialize_agat({ config_file_in => ( $shared_opts->{config} ), input => $inputFile, shared_opts => $shared_opts });

# Manage input gff file
my $format = $CONFIG->{force_gff_input_version};
if(! $format ){ $format = select_gff_format($inputFile); }
my $inputfh = open_maybe_gz($inputFile);
my $ref_in = AGAT::BioperlGFF->new(-fh => $inputfh, -gff_version => $format);

# Manage Output
if (-d $outfolder) {
  die "The output directory <$outfolder> already exists.\n";
}
else{
  my ($path,$ext);
  ($outfolder,$path,$ext) = fileparse($outfolder,qr/\.[^.]*/);
  dual_print1 "Creating the $outfolder folder\n";
  mkdir $outfolder;
}

dual_print1 "I will split the file into files containing $interval group of feature. The top feature of the group of feature is currenlty defined by <$feature_type>.\n";

#time to calcul progression
set_progression_counter($inputFile);
my $line_cpt=0;

my $count_feature=0;
my $count_file=1;
my ($file_name,$path,$ext) = fileparse($inputFile,qr/\.[^.]*/);

my $gffout = prepare_gffout( $outfolder."/".$file_name."_".$count_file.".gff");

# parse gff
while (my $feature = $ref_in->next_feature() ) {
  $line_cpt++;

  #What do we follow
  if($feature->primary_tag eq $feature_type){
    if($count_feature == $interval){
      close $gffout;
      $count_file++;
			$gffout = prepare_gffout(  $outfolder."/".$file_name."_".$count_file.".gff");
      $count_feature=0;
    }
    $count_feature++;
  }
  $gffout->write_feature($feature);

  #Display progression
  update_progression_counter($line_cpt);
}

# print fasta in asked and any
write_fasta($gffout, $ref_in);

# --- final messages ---
end_script();

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
