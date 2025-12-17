#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Pod::Usage;
use Getopt::Long;
use IO::File ;
use AGAT::AGAT;

start_script();
my $header = get_agat_header();
# -------------------------------- LOAD OPTIONS --------------------------------
my $inputFile=undef;
my $outfile=undef;
my $opt_help = 0;
my $interval=1;

# OPTION MANAGEMENT: partition @ARGV into shared vs script options via library
my ($shared_argv, $script_argv) = split_argv_shared_vs_script(\@ARGV);

# Parse script-specific options from its own list
my $script_parser = Getopt::Long::Parser->new;
$script_parser->configure('bundling','no_auto_abbrev');
if ( ! $script_parser->getoptionsfromarray(
  $script_argv,
  'file|input|gff=s' => \$inputFile,
  'i|interval=i'     => \$interval,
  'o|output=s'       => \$outfile,
  'h|help!'          => \$opt_help )  )
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

# Parse shared options (CPU, config, etc.)
my ($shared_opts) = parse_shared_options($shared_argv);

# --- Manage config ---
initialize_agat({ config_file_in => ( $shared_opts->{config} ), input => $inputFile, shared_opts => $shared_opts });

# ------------------------------------------------------------------------------

# Manage input gff file
my $format = $CONFIG->{force_gff_input_version};
if(! $format ){ $format = select_gff_format($inputFile); }
my $inputfh = open_maybe_gz($inputFile);
my $ref_in = AGAT::BioperlGFF->new(-fh => $inputfh, -gff_version => $format, -keep_fasta => !$CONFIG->{throw_fasta} );

# Manage Output
my $gffout = prepare_gffout( $outfile);
my $gffXtra=$gffout->{"_filehandle"}; #to add extra lines to gff!!

# set progression bar
set_progression_counter( $inputFile);
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
  update_progression_counter($line_cpt);
}

##Last round
 _write_bucket(\@bucket, $gffout);

# print fasta in asked and any
write_fasta($gffout, $ref_in);

# --- final counting messages ---
$count++;
if($count > 0){
  dual_print1 "$count line added !\n";
}
else{ dual_print1 "No line added !\n"; }

# --- final messages ---
end_script();

# ---------------------------- FUNCTIONS ----------------------------


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
