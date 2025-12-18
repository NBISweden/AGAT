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
my $primaryTag=undef;
my $opt_help = 0;
my $locus_tag="locus";
my $locus_cpt=1;
my $tag_in=undef;

# OPTION MANAGEMENT: partition @ARGV into shared vs script options via library
my ($shared_argv, $script_argv) = split_argv_shared_vs_script(\@ARGV);

# Parse script-specific options from its own list
my $script_parser = Getopt::Long::Parser->new;
$script_parser->configure('bundling','no_auto_abbrev');
if ( ! $script_parser->getoptionsfromarray(
    $script_argv,
    'file|input|gff=s'  => \$inputFile,
    'to|lo=s'           => \$locus_tag,
    'ti|li=s'           => \$tag_in,
    'p|type|l=s'        => \$primaryTag,
    'o|out|output=s'    => \$outfile,
    'h|help!'           => \$opt_help )  )
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

# Parse shared options (CPU, config, etc.)
my ($shared_opts) = parse_shared_options($shared_argv);

# --- Manage config ---
initialize_agat({ config_file_in => ( $shared_opts->{config} ), input => $inputFile, shared_opts => $shared_opts });

# ------------------------------------------------------------------------------

# Manage input gff file
my $format = $CONFIG->{force_gff_input_version};
if(! $format ){ $format = select_gff_format($inputFile); }
my $inputfh = open_maybe_gz($inputFile);
my $ref_in = AGAT::BioperlGFF->new(-fh => $inputfh, -gff_version => $format);

# Manage Output
my $gffout = prepare_gffout( $outfile);

#define the locus tag
if(! $locus_tag){
  $locus_tag="locus_tag";
}

# Manage $primaryTag
my @ptagList;
my $hash_level1 = $LEVELS->{'level1'};

if(! $primaryTag){
  dual_print1 "We will work on attributes from all Level1 features.\n";
  push(@ptagList, "all");
}
else{
   @ptagList= split(/,/, $primaryTag);
   foreach my $tag (@ptagList){
      if ( exists_keys ( $hash_level1, ( lc($tag) ) ) ){
        dual_print1 "We will work on attributes from <$tag> feature.\n";
      }
      else{
        dual_print1 "<$tag> feature is not a level1 feature. Current accepted value are:\n";
        foreach my $key ( keys %{$hash_level1}){
          dual_print1 $key." ";
        }
        dual_print1 "\n"; exit;
      }
   }
}

# set counter for progression bar
set_progression_counter( $inputFile);
my $line_cpt=0;

my $locus=undef;
while (my $feature = $ref_in->next_feature() ) {
  $line_cpt++;

  my $ptag = lc($feature->primary_tag());

  if ( exists_keys( $hash_level1, ( $ptag ) ) and  ( $hash_level1->{$ptag} ne "topfeature" ) ){

    # initialize locus_tag
    if ( grep( /^$ptag/, @ptagList ) or   grep( /^all/, @ptagList ) ) {

      # if locus_tag has to be the value of an existing attribute.
      if( $tag_in){
        if( $feature->has_tag($tag_in)){
          $locus = $feature->_tag_value($tag_in);
        }
        else{
          dual_print1 "No attribute $tag_in for the following feature:\n".$feature->gff_string()."\n";
          $locus = $locus_tag.$locus_cpt;$locus_cpt++;
          dual_print1 "We will use the created locus_tag value: $locus instead to name the locus!\n";
        }
      }
      else{
        $locus = $locus_tag.$locus_cpt;$locus_cpt++;
      }
    }
    else{
      $locus=undef;
    }
    # if level1 and part of those to provide locus_tag
    if($locus){
      create_or_replace_tag($feature,$locus_tag, $locus);
    }
  }
  else{
  # if not level 1 and we have to spread locus_tag to sub feature.
    if($locus){
      create_or_replace_tag($feature,$locus_tag, $locus);
    }
  }

  $gffout->write_feature($feature);

  #Display progression
  update_progression_counter($line_cpt);
}

# print fasta in asked and any
write_fasta($gffout, $ref_in);

# --- final messages ---
end_script();

# ---------------------------- FUNCTIONS ----------------------------

__END__

=head1 NAME

agat_sq_add_locus_tag.pl

=head1 DESCRIPTION

Add a shared locus tag per record. A record is all features linked by each other
by parent/children relationship (e.g Gene,mRNA,exon, CDS).

=head1 SYNOPSIS

    agat_sq_add_locus_tag.pl --gff <input file> [-o <output file>]
    agat_sq_add_locus_tag.pl --help

=head1 OPTIONS

=over 8

=item B<--gff>, B<--file> or B<--input> <file>

Input GTF/GFF file.

=item B<-p>,  B<--type> or  B<-l> <string>

Primary tag option, case insensitive, list. Allow to specied the Level1 feature types that will be handled.
By default all feature Level1 are taken into account.

=item B<--lo> or B<--to> <string>

Locus tag output, by defaut it will be called locus_tag, but using this option you can specied the name of this attribute.

=item B<--li> or B<--ti> <string>

Tag input, by default the value of the locus tag attribute will be locusX where X is an incremented number.
You can use the values of an existing attribute instead e.g the ID value: --li ID.

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
