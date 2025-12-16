#!/usr/bin/env perl

use strict;
use warnings;
use List::MoreUtils  qw(natatime);;
use Carp;
use Getopt::Long;
use Pod::Usage;
use Clone 'clone';
use AGAT::AGAT;

start_script();
my $header = get_agat_header();
# -------------------------------- LOAD OPTIONS --------------------------------
my $intronID = 1;
my $opt_file;
my $opt_output = undef;
my $opt_help = 0;
my @copyARGV = @ARGV;

# OPTION MANAGEMENT: partition @ARGV into shared vs script options via library
my ($shared_argv, $script_argv) = split_argv_shared_vs_script(\@ARGV);

# Parse script-specific options from its own list
my $script_parser = Getopt::Long::Parser->new;
$script_parser->configure('bundling', 'no_auto_abbrev');
if ( ! $script_parser->getoptionsfromarray(
    $script_argv,
    'f|gff|ref|reffile=s' => \$opt_file,
    'o|out|output=s'      => \$opt_output,
    'h|help!'             => \$opt_help,
  ) ) {
  pod2usage({
    -message => 'Failed to parse command line',
    -verbose => 1,
    -exitval => 1
  });
}

# Print Help and exit
if ($opt_help) {
  pod2usage({ -verbose => 99,
        -exitval => 0,
        -message => "$header\n" });
}

if ( ! defined( $opt_file) ) {
  pod2usage({
    -message => "$header\nMust specify at least 1 parameters:\nReference data gff3 file (--gff)\n",
    -verbose => 0,
    -exitval => 1
  });
}

# Parse shared options
my ($shared_opts) = parse_shared_options($shared_argv);

# --- Manage config ---
initialize_agat({ config_file_in => ($shared_opts->{config}), input => $opt_file, shared_opts => $shared_opts });
# -----------------------------------------------------------------------------------------------

# #######################
# START Manage Option #

my $gffout = prepare_gffout( $opt_output );

# END Manage OPTION
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

dual_print1 "$intron_added introns added\nBye Bye\n";
      #########################
      ######### END ###########
      #########################

# --- final messages ---
end_script();


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

    agat_sp_add_introns.pl --gff infile --out outFile
    agat_sp_add_introns.pl --help

=head1 OPTIONS

=over 8

=item B<--gff>, B<-f>, B<--ref> or B<-reffile>

Input GTF/GFF file.

=item  B<--out>, B<--output> or B<-o>

Output GFF3 file.

=item B<--help> or B<-h>

Display this helpful text.

=back

=head1 SHARED OPTIONS

Shared options are defined in the AGAT configuration file and can be overridden via the command line for this script only.
Common shared options are listed below; for the full list, please refer to the AGAT agat_config.yaml.

=over 8

=item B<--config>

String - Path to a custom AGAT configuration file.  
By default, AGAT uses `agat_config.yaml` from the working directory if present, otherwise the default file shipped with AGAT
(available locally via `agat config --expose`).

=item B<--cpu>, B<--core>, B<--job> or B<--thread>

Integer - Number of parallel processes to use for file input parsing (via forking).

=item B<-v> or B<--verbose>

Integer - Verbosity, choice are 0,1,2,3,4. 0 is quiet, 1 is normal, 2,3,4 is more verbose. Default 1.

=back

=head1 FEEDBACK

For questions, suggestions, or general discussions about AGAT, please use the AGAT community forum:
https://github.com/NBISweden/AGAT/discussions

=head1 BUG REPORTING

Bug reports should be submitted through the AGAT GitHub issue tracker:
https://github.com/NBISweden/AGAT/issues

=cut

AUTHOR - Jacques Dainat
