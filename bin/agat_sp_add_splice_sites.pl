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
my $cpu;
my $spliceID = 1;
my $opt_file;
my $opt_output=undef;
my $opt_help = 0;

my @copyARGV=@ARGV;
if ( !GetOptions( 'f|gff|ref|reffile=s' => \$opt_file,
                  'o|out|output=s'      => \$opt_output,
                  'c|config=s'          => \$config,
                    'thread|threads|cpu|cpus|core|cores|job|jobs=i' => \$cpu,
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
$CONFIG->{cpu} = $cpu if defined($cpu);

# #######################
# # START Manage Option #
# #######################
my $gffout = prepare_gffout( $opt_output );

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
  my ($hash_omniscient) = slurp_gff3_file_JD({ input => $opt_file });
  ### END Parse GFF input #
  #########################


my $splice_added=0;

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

          if(exists_keys($hash_omniscient, ('level3','cds',$id_l2) ) ){

          my $counterL3=-1;
          #Initialize intron to 0 to avoid error during printing results
          my $indexLast = $#{$hash_omniscient->{'level3'}{'cds'}{$id_l2}};

          my @sortedList = sort {$a->start <=> $b->start} @{$hash_omniscient->{'level3'}{'cds'}{$id_l2}};

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
            
            # Clean feature used as template
            $feature_example->frame(".");

            # Create splice template feature and assume + strand
            my $left_feature = clone($feature_example); 
            my $left_ID='splice5_added-'.$spliceID;
            my $left_fetaure_type = 'five_prime_cis_splice_site';
            my $right_feature = clone($feature_example);
            my $right_ID='splice3_added-'.$spliceID;
            my $right_fetaure_type = 'three_prime_cis_splice_site';

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
            $left_feature->primary_tag($left_fetaure_type);
            create_or_replace_tag($left_feature,'ID', $left_ID); #modify ID to replace by parent value
            $left_feature->start($tuple[0]);
            $left_feature->end($tuple[0]+1);

            # Deal with right splice
            $right_feature->primary_tag($right_fetaure_type);
            create_or_replace_tag($right_feature,'ID', $right_ID); #modify ID to replace by parent value
            $right_feature->start($tuple[1]-1);
            $right_feature->end($tuple[1]);
            
            $spliceID++; # for donnor or acceptor
            push (@{$hash_omniscient->{"level3"}{lc($left_fetaure_type)}{lc($id_l2)}}, $left_feature);
            push (@{$hash_omniscient->{"level3"}{lc($right_fetaure_type)}{lc($id_l2)}}, $right_feature);
          }
        }
      }
    }
  }

print_omniscient( {omniscient => $hash_omniscient, output => $gffout} );

print "$splice_added five_prime_cis_splice_site and $splice_added three_prime_cis_splice_site added!\nBye Bye\n";
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

The script aims to add splice sites features (five_prime_cis_splice_site and three_prime_cis_splice_site) to gtf/gff file.
The splice sites are deduced from CDS features.

=head1 SYNOPSIS

    agat_sp_add_splice_sites.pl --gff infile --out outFile
    agat_sp_add_splice_sites.pl --help

=head1 OPTIONS

=over 8

=item B<--gff>, B<-f>, B<--ref> or B<-reffile>

Input GTF/GFF file.

=item  B<--out>, B<--output> or B<-o>

Output file (default GFF3 - see config to modify output format).

=item B<-thread>, B<threads>, B<cpu>, B<cpus>, B<core>, B<cores>, B<job> or B<jobs>

Integer - Number of parallel processes to use for file input parsing (via forking).

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
