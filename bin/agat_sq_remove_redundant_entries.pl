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
# ---------------------------- OPTIONS ----------------------------
my $inputFile;
my $outfile;
my $opt_help = 0;

# OPTION MANAGEMENT: partition @ARGV into shared vs script options via library
my ($shared_argv, $script_argv) = split_argv_shared_vs_script(\@ARGV);

# Parse script-specific options from its own list
my $script_parser = Getopt::Long::Parser->new;
$script_parser->configure('bundling','no_auto_abbrev');
if ( ! $script_parser->getoptionsfromarray(
    $script_argv,
    'i|file|input|gff=s' => \$inputFile,
    'o|output=s'         => \$outfile,
    'h|help!'            => \$opt_help )  )
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

# Parse shared options (CPU, config, etc.)
my ($shared_opts) = parse_shared_options($shared_argv);

# --- Manage config ---
initialize_agat({ config_file_in => $shared_opts->{config}, input => $inputFile, shared_opts => $shared_opts });

# Manage input gff file
my $format = $CONFIG->{force_gff_input_version};
if(! $format ){ $format = select_gff_format($inputFile); }
my $inputfh = open_maybe_gz($inputFile);
my $ref_in = AGAT::BioperlGFF->new(-fh => $inputfh, -gff_version => $format);

# Manage Output
my $gffout = prepare_gffout( $outfile );

# set progression bar
set_progression_counter( $inputFile);
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
    dual_print2 "remove: ".$feature->gff_string."\n";
    next;
  }
  else{
    $gffout->write_feature($feature);
    $check{$position}++;
  }

  #Display progression
  update_progression_counter($line_cpt);
}

# print fasta in asked and any
write_fasta($gffout, $ref_in);

if($count > 0){
  dual_print1 "$count entries removed !\n";
}
else{ dual_print1 "No entry removed !\n"; }

# --- final messages ---
end_script();

# ---------------------------- FUNCTIONS ----------------------------

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
