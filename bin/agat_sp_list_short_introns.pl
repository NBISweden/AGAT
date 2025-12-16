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
start_script();
# ---------------------------- OPTIONS ----------------------------
my $opt_file;
my $INTRON_LENGTH = 10;
my $opt_output=undef;
my $opt_help = 0;

#############################
# >>>>>>>>>>>>> OPTIONS <<<<<<<<<<<<
#############################
my ($shared_argv, $script_argv) = split_argv_shared_vs_script(\@ARGV);

my $parser = Getopt::Long::Parser->new();
if ( !$parser->getoptionsfromarray(
    $script_argv,
    'f|gff|ref|reffile=s' => \$opt_file,
    'o|out|output=s'      => \$opt_output,
    'size|s=i'            => \$INTRON_LENGTH,
    'h|help!'             => \$opt_help,
  ) )
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
       -message => "$header\nMust specify at least 1 parameters:\nReference data gff3 file (--gff)\n",
       -verbose => 0,
       -exitval => 1 } );
}

#############################
# >>>>>>> Manage config <<<<<<<
#############################
my $shared_opts = parse_shared_options($shared_argv);
initialize_agat({ config_file_in => ( $shared_opts->{config} ), input => $opt_file, shared_opts => $shared_opts });

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    PARAMS    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

my $fh = prepare_fileout($opt_output);

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>     MAIN     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#PART 1
###################################
# Read input gff3 files one by one and save value in hash of list


  ######################
  ### Parse GFF input #
  my ($hash_omniscient) = slurp_gff3_file_JD({ input => $opt_file });
  ### END Parse GFF input #
  #########################

print $fh "List introns inferior to $INTRON_LENGTH nucleotides:\n\n";
print $fh "Seq_id\tGene_name\tintron_start\tintron_size\n";

my $total_intron = 0;
my %total_gene;
my %result;

  ######################
  ### Parse GFF input #
  # get nb of each feature in omniscient;
  foreach my $tag_l2 (keys %{$hash_omniscient->{'level2'}}){
    foreach my $id_l1 (keys %{$hash_omniscient->{'level2'}{$tag_l2}}){
      my $one_f2 = $hash_omniscient->{'level2'}{$tag_l2}{$id_l1}[0];

      #######################
      #get feature1 and info
      my $feature_l1=undef;
      my $tag_l1;
      foreach my $tag_level1 (keys %{$hash_omniscient->{'level1'}}){
        if (exists ($hash_omniscient->{'level1'}{$tag_level1}{$id_l1})){
          $feature_l1=$hash_omniscient->{'level1'}{$tag_level1}{$id_l1};
          $tag_l1=$tag_level1;
          last;
        }
      }
      if(! $feature_l1){ die "Problem ! We didnt retrieve the level1 feature with id $id_l1\n"; }

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
              my $intron_size = ($intronEnd - $intronStart + 1);
              if ($intron_size < $INTRON_LENGTH){
                my $seqid = $feature_l1->seq_id();

                $total_intron++;
                $total_gene{$id_l1}++;
                $result{$seqid}{$total_intron} = "$seqid\t$id_l1\t$intronStart\t$intron_size\n";

              }
            }
          }# END FOREACH L3
        }
      }
    }
  }
foreach my $seqid (keys %result){
  foreach my $cpt (keys %{ $result{$seqid} }){
    print $fh $result{$seqid}{$cpt};
  }
}

my $gene_number = keys %total_gene;
print $fh "\n$total_intron introns found for $gene_number uniq genes\n" if ($opt_output);
dual_print1 "\n$total_intron introns found for $gene_number uniq genes\n";

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

agat_sp_list_short_introns.pl

=head1 DESCRIPTION

The script aims to list all the introns inferior to a certain size.
Introns are calculated on the fly from exons. (intron feature will not be used).

=head1 SYNOPSIS

    agat_sp_list_short_introns.pl --gff infile [ --out outFile ]
    agat_sp_list_short_introns.pl --help

=head1 OPTIONS

=over 8

=item B<--gff>, B<-f>, B<--ref> or B<-reffile>

Input GTF/GFF file.

=item B<--size> or B<-s>

Minimum intron size accepted in nucleotide. All introns under this size will be reported.
Default value = 10.

=item  B<--out>, B<--output> or B<-o>

Output gff3 file where the gene incriminated will be write.


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
