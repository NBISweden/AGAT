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
my $inputFile=undef;
my $outfile=undef;
my $opt_help = 0;

# OPTION MANAGEMENT: split shared vs script-specific argv
my ($shared_argv, $script_argv) = split_argv_shared_vs_script(\@ARGV);

# Parse script options
my $script_parser = Getopt::Long::Parser->new;
$script_parser->configure('bundling','no_auto_abbrev');
if ( ! $script_parser->getoptionsfromarray(
  $script_argv,
  'file|input|gff|i=s' => \$inputFile,
  'o|out|output=s'     => \$outfile,
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

if ((!defined($inputFile)) ){
   pod2usage( { -message => "$header\nAt least 1 parameter is mandatory: -i",
                 -verbose => 0,
                 -exitval => 1 } );
}

# Parse shared options
my ($shared_opts) = parse_shared_options($shared_argv);

# --- Manage config ---
initialize_agat({ config_file_in => $shared_opts->{config}, input => $inputFile, shared_opts => $shared_opts });

# ----------------------------

# Manage input fasta file
my $format = $CONFIG->{force_gff_input_version};
if(! $format ){ $format = select_gff_format($inputFile); }
my $inputfh = open_maybe_gz($inputFile);
my $ref_in = AGAT::BioperlGFF->new(-fh => $inputfh, -gff_version => $format);

my $gffout = prepare_gffout( $outfile );

# set progression bar
set_progression_counter( $inputFile);
my $line_cpt=0;

my %hash_IDs;
my %featCount;
my %mapID;
while (my $feature = $ref_in->next_feature() ) {
  $line_cpt++;

  _uniq_ID ($feature, \%hash_IDs, \%featCount, \%mapID);

  if($feature->has_tag('Parent')){
    my $parent = lc($feature->_tag_value('Parent'));
    if(! exists($mapID{$parent})){
      dual_print1 "How is it possible ? This parent hasn't been seen before\n";
    }
     create_or_replace_tag($feature,'Parent', $mapID{$parent});
  }

  $gffout->write_feature($feature);

  #####################
  #Display progression
  update_progression_counter($line_cpt);
}

# print fasta in asked and any
write_fasta($gffout, $ref_in);

# --- final messages ---
end_script();

      # ---------------------------- FUNCTIONS ----------------------------

sub _uniq_ID{
  my ($feature, $hash_IDs, $miscCount, $mapID) = @_;


  my  $key=lc($feature->primary_tag);
  $miscCount->{$key}++;
  my $id = $key."-".$miscCount->{$key};

  while( exists_keys($hash_IDs, ($id) ) ){  #loop until we found an uniq tag
    $miscCount->{$key}++;
    $id = $key."-".$miscCount->{$key};
  }

  #push the new ID
  $hash_IDs->{$id}++;
  my $original_id = undef;
  if($feature->has_tag('ID')){
   $original_id = lc($feature->_tag_value('ID'));
  }
  else{
    $original_id = lc($id);
  }
  $mapID->{$original_id} = $id;

  # modify the feature ID with the correct one chosen
  create_or_replace_tag($feature,'ID', $id); #modify ID to replace by parent value
}

__END__

=head1 NAME

agat_sq_manage_ID.pl

=head1 DESCRIPTION

The script changes IDs to give uniq one and reflect the change in Parent attribute
of impacted features.

=head1 SYNOPSIS

    agat_sq_manage_ID.pl --gff <input file> [-o <output file>]
    agat_sq_manage_ID.pl --help

=head1 OPTIONS

=over 8

=item B<--gff>, B<--file> or B<--input>

STRING: Input GTF/GFF file.

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
