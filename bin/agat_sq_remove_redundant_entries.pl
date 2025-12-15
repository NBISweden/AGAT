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
my $ref_in = AGAT::BioperlGFF->new(-file => $inputFile, -gff_version => $format);

# Manage Output
my $gffout = prepare_gffout( $outfile );

#time to calcul progression
my $startP=time;
my $nbLine=`wc -l < $inputFile`;
$nbLine =~ s/ //g;
chomp $nbLine;
dual_print1 "$nbLine line to process...\n";
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
  if ((30 - (time - $startP)) < 0) {
    my $done = ($line_cpt*100)/$nbLine;
    $done = sprintf ('%.0f', $done);
        dual_print1 "Progression : $done % processed.\n";
    $startP= time;
  }
}

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
