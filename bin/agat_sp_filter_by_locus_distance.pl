#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long;
use Pod::Usage;
use List::MoreUtils qw(uniq);
use AGAT::AGAT;

start_script();
my $header = get_agat_header();
# -----------------------------------------------------------------------------------------------
my $outfile = undef;
my $gff = undef;
my $add_flag=undef;
my $opt_dist=500;
my $opt_help= 0;

# ---------------------------- OPTIONS ----------------------------
# Partition @ARGV into shared vs script options
my ($shared_argv, $script_argv) = split_argv_shared_vs_script(\@ARGV);

# Parse script-specific options
my $script_parser = Getopt::Long::Parser->new;
$script_parser->configure('bundling','no_auto_abbrev');
if ( ! $script_parser->getoptionsfromarray(
  $script_argv,
  'h|help!'                 => \$opt_help,
  'gff=s'                   => \$gff,
  'add_flag|af!'            => \$add_flag,
  'd|dist=i'                => \$opt_dist,
  'output|out|o=s'          => \$outfile ) )
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

if ( ! defined($gff) ){
    pod2usage( {
           -message => "$header\nAt least 1 parameter is mandatory:\nInput reference gff file (--gff)\n\n",
           -verbose => 0,
           -exitval => 1 } );
}

# Parse shared options and initialize AGAT
my ($shared_opts) = parse_shared_options($shared_argv);
initialize_agat({ config_file_in => ( $shared_opts->{config} ), input => $gff, shared_opts => $shared_opts });

# -------------------------------------------------------------------------------

######################
# Manage output file #
my $gffout = prepare_gffout( $outfile );

                #####################
                #     MAIN          #
                #####################

######################
### Parse GFF input #
my ($omniscient) = slurp_gff3_file_JD({ input => $gff });

#counters
my $geneCounter_skip=0;
my $geneCounter_ok=0;
my $total=0;
my @gene_id_ok;

my $sortBySeq = gather_and_sort_l1_location_by_seq_id($omniscient);

#get top feature first
my $top_features = get_feature_type_by_agat_value($omniscient, 'level1', 'topfeature');

foreach my $locusID ( sort keys %{$sortBySeq}){ # tag_l1 = gene or repeat etc...

  # print top features per sequence.
  foreach my $type_top_feature (keys %{$top_features}){
    if (exists_keys( $sortBySeq, ($locusID, $type_top_feature) ) ){
      foreach my $location ( @{$sortBySeq->{$locusID}{$type_top_feature}} ){
        push @gene_id_ok, lc($location->[0]); # print feature
      }
      delete  $sortBySeq->{$locusID}{$type_top_feature};
    }
  }

  foreach my $tag_l1 ( sort keys %{$sortBySeq->{$locusID}} ) {

    # Go through location from left to right ### !! if more than 2
    if (scalar @{$sortBySeq->{$locusID}{$tag_l1}} > 1){
      while ( @{$sortBySeq->{$locusID}{$tag_l1}} ){
        $total++;

        #location A
        my $location = shift  @{$sortBySeq->{$locusID}{$tag_l1}};# This location will be updated on the fly
        my $id_l1 = $location->[0];

        #it was the last location we can keep it and get out of the while loop
        if (scalar @{$sortBySeq->{$locusID}{$tag_l1}} eq 0){
          push @gene_id_ok, lc($id_l1);last;
        }

        my $continue = 1;
        my $overlap = 0;
        #loop to look at potential set of overlaping genes otherwise go through only once
        while ( $continue ){

          # Next location
          my $location2 = @{$sortBySeq->{$locusID}{$tag_l1}}[0];
          my $id2_l1 = $location2->[0];
          my $dist = $location2->[1] - $location->[2] + 1;
          dual_print2 "distance $id_l1 - id2_l1 = $dist\n";

          ############################
          #deal with overlap
          if ( ($location->[1] <= $location2->[2]) and ($location->[2] >= $location2->[1])){

                if( ! $overlap){

                  foreach my $tag_level1 (keys %{$omniscient->{'level1'}}){
                    if (exists_keys($omniscient, ('level1', $tag_level1, lc($id_l1) ) ) ){
                      my $level1_feature = $omniscient->{'level1'}{$tag_level1}{lc($id_l1)};
                      add_info($level1_feature, 'O');
                    }
                  }
                }

                $overlap=1;

                foreach my $tag_level1 (keys %{$omniscient->{'level1'}}){
                  if (exists_keys($omniscient, ('level1', $tag_level1, lc($id2_l1) ) ) ){
                    my $level1_feature = $omniscient->{'level1'}{$tag_level1}{lc($id2_l1)};
                    add_info($level1_feature, 'O');
                  }
                }

                if($location2->[2] < $location->[2]){
                  my $tothrow = shift  @{$sortBySeq->{$locusID}{$tag_l1}};# Throw location B. We still need to use location A to check the left extremity of the next locus
                                                                            #location A  -------------------------                                --------------------------
                                                                            #location B                ---------
                  $total++;
                  next;
                }
                else{
                                                                            # We need to use the location B to check the left extremity of the next locus
                                                                            #location A  -------------------------                                --------------------------
                                                                            #location B                ------------------------
                  last;
                }

          }
          #
          ############################
          $continue = 0;

          # locus distance is under minimum distance
          if( $dist < $opt_dist)  {

            foreach my $tag_level1 (keys %{$omniscient->{'level1'}}){
              if (exists_keys($omniscient, ('level1', $tag_level1, lc($id_l1) ) ) ){
                my $level1_feature = $omniscient->{'level1'}{$tag_level1}{lc($id_l1)};
                add_info($level1_feature, 'R'.$dist);
              }
            }

            foreach my $tag_level1 (keys %{$omniscient->{'level1'}}){
              if (exists_keys($omniscient, ('level1', $tag_level1, lc($id2_l1) ) ) ){
                my $level1_feature = $omniscient->{'level1'}{$tag_level1}{lc($id2_l1)};
                add_info($level1_feature, 'L'.$dist);
              }
            }
          }

          # distance with next is ok but we have to check what was the result with the previous locus
          else{
           foreach my $tag_level1 (keys %{$omniscient->{'level1'}}){
              if (exists_keys($omniscient, ('level1', $tag_level1, lc($id_l1) ) ) ){
                my $level1_feature = $omniscient->{'level1'}{$tag_level1}{lc($id_l1)};
                if(! $level1_feature->has_tag('low_dist')){
                  $geneCounter_ok ++;
                  push @gene_id_ok, lc($id_l1);
                }
              }
            }
          }
        }
      }
    }
    # if we have only 1 locations
    else{
      #push @gene_id_ok, $sortBySeq->{$locusID}{$tag_l1}->[0];
      my $location = shift  @{$sortBySeq->{$locusID}{$tag_l1}};
      my $id_l1 = $location->[0];
      push @gene_id_ok, lc($id_l1);
    }
  }
}

if($add_flag){
  print_omniscient( {omniscient => $omniscient, output => $gffout} );
}
else{
  print_omniscient_from_level1_id_list( {omniscient => $omniscient, level_id_list =>\@gene_id_ok, output => $gffout} );
}

#END messages
my $string_to_print = "Results:\n".
  "Total number investigated: $total\n".
  "Number of skipped loci: $geneCounter_skip\n".
  "Number of loci with distance to the surrounding loci over $opt_dist: $geneCounter_ok \n";
dual_print1 $string_to_print;

end_script();
#######################################################################################################################
        ####################
         #     METHODS    #
          ################
           ##############
            ############
             ##########
              ########
               ######
                ####
                 ##


sub add_info{
  my ($feature, $value)=@_;

  if($feature->has_tag('low_dist')){
    $feature->add_tag_value('low_dist', $value);
    dual_print2 $feature->_tag_value('ID')." add $value\n";
  }
  else{
    create_or_replace_tag($feature, 'low_dist', $value);
    $geneCounter_skip++;
    dual_print2 $feature->_tag_value('ID')." create $value\n";
  }

}
__END__

=head1 NAME

agat_sp_filter_by_locus_distance.pl

=head1 DESCRIPTION

The script aims to remove or flag loci that are too close to each other.
Close loci are important to remove when training abinitio tools in order
to train intergenic region properly. Indeed if intergenic region
(surrouneded part of a locus) contain part of another locus,
the training on intergenic part will be biased.

=head1 SYNOPSIS

    agat_sp_filter_by_locus_distance.pl -gff infile.gff [ -o outfile ]
    agat_sp_filter_by_locus_distance.pl --help

=head1 OPTIONS

=over 8

=item B<-gff> <file>

Input GTF/GFF file.

=item B<--dist> or B<-d> <int>

The minimum inter-loci distance to allow. No default (will not apply
filter by default).

=item B<--af> or B<--add_flag>

Instead of filter the result into two output files, write only one and add the flag <low_dist> in the gff.(tag = Lvalue or tag = Rvalue  where L is left and R right and the value is the distance with accordingle the left or right locus)

=item B<-o>, B<--out> or B<--output> <file>

Output file to create (default GFF3 - see config to modify output format).
If no output file is specified, the output will be written to STDOUT.

=item B<-h> or B<--help>

Display this helpful text.

=back

=head1 SHARED OPTIONS

Shared options are defined in the AGAT configuration file and can be overridden via the command line for this script only.
Common shared options are listed below; for the full list, please refer to the AGAT agat_config.yaml.

=over 8

=item B<--config> <file>

Path to a custom AGAT configuration file.  
By default, AGAT uses `agat_config.yaml` from the working directory if present, otherwise the default file shipped with AGAT
(available locally via `agat config --expose`).

=item B<--cpu>, B<--core>, B<--job> or B<--thread> <int>

Number of parallel processes to use for file input parsing (via forking).

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
