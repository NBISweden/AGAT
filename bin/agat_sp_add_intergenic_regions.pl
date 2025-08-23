#!/usr/bin/env perl

use strict;
use warnings;
use POSIX qw(strftime);
use List::MoreUtils  qw(natatime);;
use Carp;
use Getopt::Long::Descriptive;
use Pod::Usage;
use Clone 'clone';
use AGAT::AGAT;

my $header = get_agat_header();
my ( $opt, $usage, $config ) = AGAT::AGAT::describe_script_options(
    $header,
    [ 'gff|f|ref=s', 'Input GTF/GFF file', { required => 1 } ],
);

my $opt_file   = $opt->gff;
my $opt_output = $opt->out;
my $intergenicID = 1;

my $log;
if ( my $log_name = $config->{log_path} ) {
    open( $log, '>', $log_name ) or die "Can not open $log_name for printing: $!";
    dual_print( $log, $header, 0 );
}

my $gffout = prepare_gffout($config, $opt_output);

#                         #######################
# >>>>>>>>>>>>>>>>>>>>>>>>#        MAIN         #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#                         #######################

  ######################
  ### Parse GFF input #
  my ($omniscient, $hash_mRNAGeneLink) = slurp_gff3_file_JD({ input => $opt_file,
                                                              config => $config });
  dual_print($log, "Parsing Finished\n", $config->{verbose});
  ### END Parse GFF input #
  #########################

if(! exists_keys($omniscient,('level1', "gene") ) ){
  dual_print($log, "No gene feature found in $opt_file, intergenic regions cannot be determinded!\n", $config->{verbose});
  exit 0;
}

# Gather all Level1 features
my $sortBySeq = gather_and_sort_l1_by_seq_id_for_l1type($omniscient, 'gene');

# --------------------------- COLLECT GENE LOCATIONS -----------------------
dual_print($log, "Now colleting the gene locations\n", $config->{verbose});
my $flattened_locations = {};
foreach my $locusID ( sort keys %{$sortBySeq}){ # tag_l1 = gene or repeat etc...
  # check if gene  exits for this sequence
  if( exists_keys($sortBySeq,($locusID, "gene") ) ){

    # Go through location from left to right ### !! if not empty
    while ( my $gene_feature = shift @{$sortBySeq->{$locusID}{"gene"}} ){

      # Define location l1
      my $current_location_l1 = [$gene_feature->start(), $gene_feature->end()];

      # save appending a list
      push @{$flattened_locations->{$locusID}}, $current_location_l1 ;
    }
  }
}

# --------------------------- FIX OVERLAPPING LOCATIONS -----------------------
# Will merge locations that overlap

dual_print($log, "Now flattening the locations\n", $config->{verbose});
foreach my $locusID (  keys %{$flattened_locations} ){

  my @newlocations;
  my $previous_location = undef;
  foreach my $location ( sort {$a->[0] <=> $b->[0]} @{$flattened_locations->{$locusID}} ){

    # first round
    if (! $previous_location){
        push @newlocations, $location;
        $previous_location = $location;
    }
    # Not first round
    else{
      #  location A  -------------------------
      #  location B           -------------------------
      if ( ($previous_location->[0] <= $location->[1]) and ($previous_location->[1] >= $location->[0])){

        #  location A  -------------------------
        #  location B           -------------------------
        if($previous_location->[1] <= $location->[1]){
          # take back last location pushed in the array
          $previous_location = pop @newlocations;
          $previous_location  = [$previous_location->[0], $location->[1]];
          # push back into the array the previous location that has been modified
          push @newlocations, $previous_location ;
        }
      }
      else{
        push @newlocations, $location ;
        $previous_location = $location;
      }
    }
  }
  @{$flattened_locations->{$locusID}} = @newlocations ;
}

# --------------------------- NOW creating intergenic location -----------------------
dual_print($log, "Now creating intergenic regions\n", $config->{verbose});
my $intergenic_added=0;
# Go through location from left to right ### !! if not empty
foreach my $locusID ( sort keys %{$flattened_locations}){ # tag_l1 = gene or repeat etc...
  my $round = 0;
  my $previous_location = undef;
  my $last_location = undef;

  while ( my $location = shift @{$flattened_locations->{$locusID}} ){
    $round++;
     
    $last_location = 1 if (scalar @{$flattened_locations->{$locusID}} == 0);
    if($round == 1 or $last_location){
      $previous_location=$location;
      next;
    }

    my $uniq_ID = "intergenic_region_".$intergenicID;
    my $intergenic_feature = Bio::SeqFeature::Generic->new( -seq_id => $locusID, 
                                                            -source_tag => "AGAT",
                                                            -primary_tag => 'intergenic_region' , 
                                                            -start => $previous_location->[1]+1,  
                                                            -end => $location->[0]-1, 
                                                            -frame => ".", 
                                                            -strand => "." , 
                                                            -tag => { 'ID' => $uniq_ID  }) ;
    # add feature in omniscient
    $omniscient->{"level1"}{"intergenic_region"}{lc($uniq_ID)}=$intergenic_feature; 
    
    # update local variables
    $intergenic_added++;
    $intergenicID++;
    $previous_location=$location;
  }     
}

# print result
print_omniscient( {omniscient => $omniscient, output => $gffout} );

dual_print($log, "$intergenic_added intergenic_region added!\nBye Bye\n", $config->{verbose});
      #########################
      ######### END ###########
      #########################

#######################################################################################################################

        ####################
         #     methods    #
          ################
           ##############
            ############
             ##########
              ########
               ######
                ####
                 ##

__END__

=head1 NAME

agat_sp_add_intergenic_regions.pl

=head1 DESCRIPTION

The script aims to add intergenic features (intergenic_region) to gtf/gff file.
The intergenic regions are deduced from gene features (feature type gene from the 3rd column).

=head1 SYNOPSIS

    agat_sp_add_intergenic_regions.pl --gff infile --out outFile
    agat_sp_add_intergenic_regions.pl --help

=head1 OPTIONS

=over 8

=item B<--gff>, B<-f> or B<--ref>

Input GTF/GFF file.

=item  B<--out>, B<--output> or B<-o>

Output file (default GFF3 - see config to modify output format).

=item B<-c> or B<--config>

String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any, 
otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
The --config option gives yo the possibility to use your own AGAT config file (located elsewhere or named differently).

=item B<-v> or B<--verbose>

Add verbosity

=item B<-h> or B<--help>

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
