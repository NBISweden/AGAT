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
my $opt;
eval {
    ( $opt ) = describe_options(
        "$header\n\n%c %o",
        [ 'gff|f|ref|reffile=s', 'Input GTF/GFF file', { required => 1 } ],
        [ 'help|h',              'Show this help', { shortcircuit => 1 } ],
        { getopt_conf => ['pass_through'] },
    );
    1;
} or pod2usage( { -message => $@, -exitstatus => 1, -verbose => 1 } );
pod2usage( { -verbose => 99, -exitstatus => 0, -message => "$header\n" } )
  if $opt->help;

my $opt_file = $opt->gff;
my $config   = resolve_common_options( \@ARGV );
my $intronID = 1;

# #######################
# # START Manage Option #
# #######################
my $gffout = prepare_gffout( $config, $config->{output} );

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
  my ($hash_omniscient, $hash_mRNAGeneLink) = slurp_gff3_file_JD({ input => $opt_file,
                                                                   config => $config });
  print("Parsing Finished\n\n") if $config->{verbose};
  ### END Parse GFF input #
  #########################


my $intron_added=0;

  ######################
  ### Parse GFF input #
  # get nb of each feature in omniscient;
  foreach my $tag_l2 (sort keys %{$hash_omniscient->{'level2'}}){
    foreach my $id_l1 (sort keys %{$hash_omniscient->{'level2'}{$tag_l2}}){
      my $one_f2 = $hash_omniscient->{'level2'}{$tag_l2}{$id_l1}[0];

      #######################
      #get feature1 and info
      my $feature_l1=undef;
      my $tag_l1;
      foreach my $tag_level1 (sort keys %{$hash_omniscient->{'level1'}}){
        if (exists ($hash_omniscient->{'level1'}{$tag_level1}{$id_l1})){
          $feature_l1=$hash_omniscient->{'level1'}{$tag_level1}{$id_l1};
          $tag_l1=$tag_level1;
          last;
        }
      }
      if(! $feature_l1){print "Problem ! We didnt retrieve the level1 feature with id $id_l1\n";exit;}

      #####
      # get all level2
      my $All_l2_single=1;
      foreach my $feature_l2 ( @{$hash_omniscient->{'level2'}{$tag_l2}{$id_l1}} ){

        my @introns=();
        my $feature_example;

        ######
        #get all level3
        my $id_l2=lc($feature_l2->_tag_value('ID'));

          if(exists_keys($hash_omniscient, ('level3','exon',$id_l2) ) ){

          my $counterL3=-1;
          #Initialize intron to 0 to avoid error during printing results
          my $indexLast = $#{$hash_omniscient->{'level3'}{'exon'}{$id_l2}};

          my @sortedList = sort {$a->start <=> $b->start} @{$hash_omniscient->{'level3'}{'exon'}{$id_l2}};

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
						$intron_added++;
            my $intron_feature = clone($feature_example);
            $intron_feature->primary_tag('intron');
            my $ID='intron_added-'.$intronID;
            $intronID++;
            create_or_replace_tag($intron_feature,'ID', $ID); #modify ID to replace by parent value
            $intron_feature->start($tuple[0]);
            $intron_feature->end($tuple[1]);
            push (@{$hash_omniscient->{"level3"}{'intron'}{lc($id_l2)}}, $intron_feature);
          }
        }
      }
    }
  }

print_omniscient( {omniscient => $hash_omniscient, output => $gffout} );

print "$intron_added introns added\nBye Bye\n" if $config->{verbose};
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

agat_sp_add_introns.pl

=head1 DESCRIPTION

The script aims to add intron features to gtf/gff file without intron features.

=head1 SYNOPSIS

    agat_sp_add_introns.pl --gff infile --out outFile --verbose 0
    agat_sp_add_introns.pl --help

=head1 OPTIONS

=over 8

=item B<--gff>, B<-f>, B<--ref> or B<-reffile>

Input GTF/GFF file.

=item  B<--out>, B<--output> or B<-o>

Output GFF3 file.

=item B<-c> or B<--config>

String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any, 
otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

=item B<-v> or B<--verbose>

Integer - Verbosity level (0-4). Command line overrides the config file. [Default 1]

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
