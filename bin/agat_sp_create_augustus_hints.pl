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
my $opt_file;
my $opt_output=undef;
my $opt_help = 0;

my @copyARGV=@ARGV;
if ( !GetOptions( 'f|gff|ref|reffile=s' => \$opt_file,
                  'o|out|output=s'      => \$opt_output,
                  'c|config=s'          => \$config,
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
$config = get_agat_config({config_file_in => $config});

# #######################
# # START Manage Option #
# #######################
my $gffout = prepare_gffout($config, $opt_output);

# #####################################
# # END Manage OPTION
# #####################################

#                         #######################
# >>>>>>>>>>>>>>>>>>>>>>>>#        MAIN         #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#                         #######################

#PART 1
###################################
# Read input gff3 files one by one and save value in hash of list


######################
### Parse GFF input #
my ($omniscient, $hash_mRNAGeneLink) = slurp_gff3_file_JD({ input => $opt_file,
                                                                  config => $config });
print("Parsing Finished\n\n");
### END Parse GFF input #
#########################

# ---------------------------------------------------------
# ------------------- Introns and Splices -------------------
# ---------------------------------------------------------

# counter for creating IDs
my $spliceID = 1;
my $intronID = 1;

# counter for tracking created features
my $splice_added=0;
my $intron_added=0;

# get nb of each feature in omniscient;
foreach my $tag_l2 (sort keys %{$omniscient->{'level2'}}){
  foreach my $id_l1 (sort keys %{$omniscient->{'level2'}{$tag_l2}}){
    my $one_f2 = $omniscient->{'level2'}{$tag_l2}{$id_l1}[0];

    #######################
    #get feature1 and info
    my $feature_l1=undef;
    my $tag_l1;
    foreach my $tag_level1 (sort keys %{$omniscient->{'level1'}}){
      if (exists ($omniscient->{'level1'}{$tag_level1}{$id_l1})){
        $feature_l1=$omniscient->{'level1'}{$tag_level1}{$id_l1};
        $tag_l1=$tag_level1;
        last;
      }
    }
    if(! $feature_l1){print "Problem ! We didnt retrieve the level1 feature with id $id_l1\n";exit;}

    #####
    # get all level2
    my $All_l2_single=1;
    foreach my $feature_l2 ( @{$omniscient->{'level2'}{$tag_l2}{$id_l1}} ){

      my @introns=();
      my $feature_example;

      ######
      #get all level3
      my $id_l2=lc($feature_l2->_tag_value('ID'));

        if(exists_keys($omniscient, ('level3','cds',$id_l2) ) ){

        my $counterL3=-1;
        #Initialize intron to 0 to avoid error during printing results
        my $indexLast = $#{$omniscient->{'level3'}{'cds'}{$id_l2}};

        my @sortedList = sort {$a->start <=> $b->start} @{$omniscient->{'level3'}{'cds'}{$id_l2}};

          foreach my $feature_l3 ( @sortedList ){

            #count number feature of tag_l3 type
            $counterL3++;

            ################
            #Manage Introns#
            # from the second intron to the last (from index 1 to last index of the table sortedList)
            # We go inside this loop only if we have more than 1 feature.
            if($counterL3 > 0 and $counterL3 <= $indexLast){
              my $intronStart = $sortedList[$counterL3-1]->end+1;
              my $intronEnd = $sortedList[$counterL3]->start-1;
              push @introns, ($intronStart, $intronEnd);
              $feature_example=clone($sortedList[$counterL3]);
            }
          }# END FOREACH L3
        }

      #Now add introns features
      if(@introns){
        my $it = natatime 2, @introns;
        while (my @tuple = $it->()) {

          $splice_added++;
          $intron_added++;

          # --- Create introns ---
          my $ID='intron_added-'.$intronID; $intronID++;
          my $intron_feature = clean_clone( { omniscient => $omniscient,
                                feature => $feature_example,
                                new_id => $ID,
                                new_primary_tag => "intron"
                                } );
          $intron_feature->start($tuple[0]);
          $intron_feature->end($tuple[1]);
          push (@{$omniscient->{"level3"}{'intron'}{lc($id_l2)}}, $intron_feature);

          # --- Create Splice sites ---
          # Create splice template feature and assume + strand
          my $left_ID='splice5_added-'.$spliceID;
          my $left_fetaure_type = 'dss';
          my $right_ID='splice3_added-'.$spliceID;
          my $right_fetaure_type = 'ass';

          # Flip if minus strand
          if ( ($feature_example->strand() eq "-" ) or ($feature_example->strand() eq "-1" ) ) {            
            # flip ID              
            my $tmp = $left_ID;
            $left_ID = $right_ID;
            $right_ID = $tmp;
            # flip feature type
            my $tmp_ft = $left_fetaure_type;
            $left_fetaure_type = $right_fetaure_type;
            $right_fetaure_type = $tmp_ft;
          }

          # Deal with left splice
          my $left_feature = clean_clone( { omniscient => $omniscient,
                                feature => $feature_example,
                                new_id => $left_ID,
                                new_primary_tag => $left_fetaure_type
                                } );
          $left_feature->start($tuple[0]);
          $left_feature->end($tuple[0]+1);

          # Deal with right splice
          my $right_feature = clean_clone( { omniscient => $omniscient,
                                feature => $feature_example,
                                new_id => $right_ID,
                                new_primary_tag => $right_fetaure_type
                                } );
          $right_feature->start($tuple[1]-1);
          $right_feature->end($tuple[1]);
          
          $spliceID++; # for donnor or acceptor
          push (@{$omniscient->{"level3"}{lc($left_fetaure_type)}{lc($id_l2)}}, $left_feature);
          push (@{$omniscient->{"level3"}{lc($right_fetaure_type)}{lc($id_l2)}}, $right_feature);
        }
      }
    }
  }
}

# ---------------------------------------------------------
#    ------------------- Intergenic -------------------
# ---------------------------------------------------------
if(! exists_keys($omniscient,('level1', "gene") ) ){
  print "No gene feature found in $opt_file, intergenic regions cannot be determinded!";
} else {

  # counter for creating intergenic feature IDs
  my $intergenicID = 1;

  # collect level1
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

  print "$intergenic_added intergenic_region added!\n";
}



# print results
print_omniscient( {omniscient => $omniscient, output => $gffout} );

print "$splice_added five_prime_cis_splice_site and $splice_added three_prime_cis_splice_site added!\n";
print "$intron_added introns added\n";
print "Bye Bye\n";
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

agat_sp_add_splice_sites.pl

=head1 DESCRIPTION

The script aims to create hints that can be used to train Augustus.
There are 16 types of hints accepted by Augustus: start, stop, tss, tts, ass, dss, exonpart, exon, intronpart, intron, CDSpart, CDS, UTRpart, UTR, irpart, nonexonpart. 
AGAT will prepare the following hints: start, stop, tss, tts, ass, dss, exon, intron, CDS, UTR, irpart

start: translation start (start codon), specifies an interval that contains the start codon. The interval can be larger than 3bp, in which case every ATG in the interval gets a bonus. The highest bonus is given to ATGs in the middle of the interval, the bonus fades off towards the ends.
stop: translation end (stop codon), see 'start'
tss: transcription start site, see 'start'
tts: transcription termination site, see 'start'
ass: acceptor (3') splice site, the last intron position, for only approximately known ass an interval can be specified
dss: donor (5') splice site, the first intron position, for only approximately known dss an interval can be specified
exonpart: part of an exon in the biological sense. The bonus applies only to exons that contain the interval from the hint. Just overlapping means no bonus at all. The malus applies to every base of an exon. Therefore the malus for an exon is exponential in the length of an exon: malus=exonpartmalus^length. Therefore the malus should be close to 1, e.g. 0.99.
exon: exon in the biological sense. Only exons that exactly match the hint get a bonus. Exception: The exons that contain the start codon and stop codon. This malus applies to a complete exon independent of its length.
intronpart: introns both between coding and non-coding exons. The bonus applies to every intronic base in the interval of the hint.
intron: An intron gets the bonus if and only if it is exactly as in the hint.
CDSpart: part of the coding part of an exon. (CDS = coding sequence)
CDS: coding part of an exon with exact boundaries. For internal exons of a multi exon gene this is identical to the biological boundaries of the exon. For the first and the last coding exon the boundaries are the boundaries of the coding sequence (start, stop).
UTR: exact boundaries of a UTR exon or the untranslated part of a partially coding exon.
UTRpart: The hint interval must be included in the UTR part of an exon.
irpart: The bonus applies to every base of the intergenic region. If UTR prediction is turned on (--UTR=on) then UTR is considered genic. If you choose against the usual meaning the bonus of irparts to be much smaller than 1 in the configuration file you can force AUGUSTUS to not predict an intergenic region in the specified interval. This is useful if you want to tell AUGUSTUS that two distant exons belong to the same gene, when AUGUSTUS tends to split that gene into smaller genes.
nonexonpart: intergenic region or intron. The bonus applies to very non-exon base that overlaps with the interval from the hint. It is geometric in the length of that overlap, so choose it close to 1.0. This is useful as a weak kind of masking, e.g. when it is unlikely that a retroposed gene contains a coding region but you do not want to completely forbid exons.
genicpart: everything that is not intergenic region, i.e. intron or exon or UTR if applicable. The bonus applies to every genic base that overlaps with the interval from the hint. This can be used in particular to make Augustus predict one gene between positions a and b if a and b are experimentally confirmed to be part of the same gene, e.g. through ESTs from the same clone. alias: nonirpart

=head1 SYNOPSIS

    agat_sp_add_splice_sites.pl --gff infile --out outFile
    agat_sp_add_splice_sites.pl --help

=head1 OPTIONS

=over 8

=item B<--gff>, B<-f>, B<--ref> or B<-reffile>

Input GTF/GFF file.

=item  B<--out>, B<--output> or B<-o>

Output file (default GFF3 - see config to modify output format).

=item B<-c> or B<--config>

String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any, 
otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
The --config option gives yo the possibility to use your own AGAT config file (located elsewhere or named differently).

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
