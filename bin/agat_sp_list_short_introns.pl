#!/usr/bin/env perl

use strict;
use warnings;
use POSIX qw(strftime);
use List::MoreUtils  qw(natatime);;
use Carp;
use Clone 'clone';
use AGAT::AGAT;

my $header = get_agat_header();
my @copyARGV = @ARGV;
my ( $opt, $usage, $config ) = AGAT::AGAT::describe_script_options(
    $header,
    [ 'gff|f|ref|reffile=s', 'Reference data gff3 file', { required => 1 } ],
    [ 'size|s=i', 'Intron size threshold', { default => 10, callbacks => { positive => sub { shift() > 0 or die 'Intron size threshold must be positive' } } } ],
);

my $opt_file      = $opt->gff;
my $INTRON_LENGTH = $opt->size;
my $opt_output    = $config->{output};
my $verbose       = $config->{verbose};

my $log;
if ( my $log_name = $config->{log_path} ) {
    open( $log, '>', $log_name ) or die "Can not open $log_name for printing: $!";
    dual_print( $log, $header, 0 );
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    PARAMS    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

my $fh = prepare_fileout($opt_output);

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>     MAIN     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#PART 1
###################################
# Read input gff3 files one by one and save value in hash of list


  ######################
  ### Parse GFF input #
  my ($hash_omniscient, $hash_mRNAGeneLink) = slurp_gff3_file_JD({ input => $opt_file,
                                                                   config => $config
                                                              });
  dual_print( $log, "Parsing Finished\n\n", $verbose );
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
      if(! $feature_l1){
        my $msg = "Problem ! We didnt retrieve the level1 feature with id $id_l1\n";
        dual_print( $log, $msg, 0 );
        warn $msg if $verbose;
        exit;
      }

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
print $fh "\n$total_intron introns found for $gene_number uniq genes\n";
close $log if $log;
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
