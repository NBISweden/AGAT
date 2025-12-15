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
my $outformat=undef;
my $opt_help = 0;

# OPTION MANAGEMENT: split shared vs script-specific argv
my ($shared_argv, $script_argv) = split_argv_shared_vs_script(\@ARGV);

# Parse script options
my $script_parser = Getopt::Long::Parser->new;
$script_parser->configure('bundling','no_auto_abbrev');
if ( ! $script_parser->getoptionsfromarray(
  $script_argv,
  'file|input|gff|i=s' => \$inputFile,
  'of=i'               => \$outformat,
  'o|output=s'         => \$outfile,
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
my $ref_in = AGAT::BioperlGFF->new(-file => $inputFile, -gff_version => $format);

my $gffout = prepare_gffout( $outfile );

#time to calcul progression
my $startP=time;
my $nbLine=`wc -l < $inputFile`;
$nbLine =~ s/ //g;
chomp $nbLine;
dual_print1 "$nbLine line to process...\n";

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
  if ((30 - (time - $startP)) < 0) {
    my $done = ($line_cpt*100)/$nbLine;
    $done = sprintf ('%.0f', $done);
        dual_print1 "\rProgression : $done % processed.\n";
    $startP= time;
  }
}

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

=item B<--of>

Output format, if no ouput format is given, the same as the input one detected will be used. Otherwise you can force to have a gff version 1 or 2 or 3 by giving the corresponding number.

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
