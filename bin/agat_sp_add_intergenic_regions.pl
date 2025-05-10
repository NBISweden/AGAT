#!/usr/bin/env perl

use strict;
use warnings;
use POSIX qw(strftime);
use List::MoreUtils  qw(natatime);;
use Carp;
use Getopt::Long;
use Pod::Usage;
use Clone 'clone';
use AGAT::AGAT;

my $header = get_agat_header();
my $config;
my $threads;
my $verbose;
my $intergenicID = 1;
my $opt_file;
my $opt_output=undef;
my $opt_help = 0;

my @copyARGV=@ARGV;
if ( !GetOptions( 'f|gff|ref=s'         => \$opt_file,
                  'o|out|output=s'      => \$opt_output,
                  'c|config=s'          => \$config,
                  'v|verbose!'          => \$verbose,
                  'h|help!'             => \$opt_help ) )
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

if ( ! defined( $opt_file) ) {
    pod2usage( {
           -message => "$header\nMust specify at least 1 parameters:\nReference data GFF/GTF file (--gff)\n",
           -verbose => 0,
           -exitval => 1 } );
}

# --- Manage config ---
initialize_agat({ config_file_in => $config, input => $opt_file });
$CONFIG->{threads} = $threads if defined($threads);

# # START Manage Option #
# #######################
my $gffout = prepare_gffout( $opt_output );

#                         #######################
# >>>>>>>>>>>>>>>>>>>>>>>>#        MAIN         #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#                         #######################

  ######################
  ### Parse GFF input #
  my ($omniscient, $hash_mRNAGeneLink) = slurp_gff3_file_JD({ input => $opt_file });
  ### END Parse GFF input #
  #########################

if(! exists_keys($omniscient,('level1', "gene") ) ){
  print "No gene feature found in $opt_file, intergenic regions cannot be determinded!", exit 0;
}

# Gather all Level1 features
my $sortBySeq = gather_and_sort_l1_by_seq_id_for_l1type($omniscient, 'gene');

# --------------------------- COLLECT GENE LOCATIONS -----------------------
print "Now colleting the gene locations\n" if ($verbose);
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

print "Now flattening the locations\n" if ($verbose);
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
print "Now creating intergenic regions\n" if ($verbose);
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

print "$intergenic_added intergenic_region added!\nBye Bye\n";
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
