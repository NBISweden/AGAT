#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Pod::Usage;
use IO::File ;
use AGAT::AGAT;

my $header = get_agat_header();
my ( $opt, $usage, $config ) = AGAT::AGAT::describe_script_options(
    $header,
    [ 'file|input|gff|i=s', 'Input GTF/GFF file', { required => 1 } ],
    # [ 'of=i', 'Output format' ], # currently ineffective (as of https://github.com/NBISweden/AGAT/commit/a14978012da04e62f83ea0d5bb7d6861de361ea5)
);

my $inputFile = $opt->file;
my $outfile   = $config->{output};

my $log;
if ( my $log_name = $config->{log_path} ) {
    open( $log, '>', $log_name )
      or die "Can not open $log_name for printing: $!";
    dual_print( $log, $header,  3 );
}

my $start_run = time();

# Manage input fasta file
my $format = $config->{force_gff_input_version};
if ( !$format ) { $format = select_gff_format($inputFile); }
my $ref_in = AGAT::BioperlGFF->new( -file => $inputFile, -gff_version => $format );

my $gffout = prepare_gffout( $config, $outfile );

#time to calcul progression
my $startP=time;
my $nbLine=`wc -l < $inputFile`;
$nbLine =~ s/ //g;
chomp $nbLine;
dual_print( $log, "$nbLine line to process...\n");

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
      warn "How is it possible ? This parent hasn't been seen before\n" if $config->{verbose};
    }
     create_or_replace_tag($feature,'Parent', $mapID{$parent});
  }

  $gffout->write_feature($feature);

  #####################
  #Display progression
  if ((30 - (time - $startP)) < 0) {
    my $done = ($line_cpt*100)/$nbLine;
    $done = sprintf ('%.0f', $done);
        dual_print( $log, "\rProgression : $done % processed.\n");
    $startP= time;
  }
}

##Last round
my $end_run = time();
my $run_time = $end_run - $start_run;
dual_print( $log, "Job done in $run_time seconds\n");



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

# currently ineffective (as of https://github.com/NBISweden/AGAT/commit/a14978012da04e62f83ea0d5bb7d6861de361ea5)
# =item B<--of>
# 
# Output format, if no ouput format is given, the same as the input one detected will be used. Otherwise you can force to have a gff version 1 or 2 or 3 by giving the corresponding number.


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
