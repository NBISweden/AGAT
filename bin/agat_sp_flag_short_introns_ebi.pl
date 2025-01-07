#!/usr/bin/env perl

use strict;
use warnings;
use POSIX qw(strftime);
use File::Basename;
use Carp;
use Getopt::Long;
use IO::File;
use Pod::Usage;
use AGAT::AGAT;

my $header = get_agat_header();
my $config;
my $opt_file;
my $opt_output=undef;
my $verbose=undef;
my $Xsize=10;
my $opt_help = 0;

my @copyARGV=@ARGV;
if ( !GetOptions( 'f|gff|ref|reffile=s' => \$opt_file,
                  'o|out|output=s'      => \$opt_output,
                  'v|verbose!'          => \$verbose,
                  'i|intron_size=i'     => \$Xsize,
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

if ( ! defined($opt_file) ) {
    pod2usage( {
           -message => "$header\nMust specify at least 1 parameters:\nReference data gff3 file (--gff)\n",
           -verbose => 0,
           -exitval => 1 } );
}

# --- Manage config ---
$config = get_agat_config({config_file_in => $config});

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    PARAMS    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

my $ostreamReport_file;
if (defined($opt_output) ) {
  my ($filename,$path,$ext) = fileparse($opt_output,qr/\.[^.]*/);
  $ostreamReport_file = $path.$filename."_report.txt";
}

my $gffout = prepare_gffout($config, $opt_output);
my $ostreamReport = prepare_fileout($ostreamReport_file);


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    EXTRA     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# Print info
my $string1 = strftime "%m/%d/%Y at %Hh%Mm%Ss", localtime;
$string1 .= "\n\nusage: $0 @copyARGV\n\n";

print $ostreamReport $string1;
if($opt_output){print $string1;}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>     MAIN     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

######################
### Parse GFF input #
print "Reading ".$opt_file,"\n";
my ($hash_omniscient, $hash_mRNAGeneLink) = slurp_gff3_file_JD({ input => $opt_file,
                                                                 config => $config
                                                              });
print("Parsing Finished\n\n");
### END Parse GFF input #
#########################

my $nb_cases=0;
my $tag = "pseudo";
######################
### Parse GFF input #
foreach my $tag_l1 (keys %{$hash_omniscient->{'level1'}}){
  foreach my $id_l1 (keys %{$hash_omniscient->{'level1'}{$tag_l1}}){
    my $shortest_intron=10000000000;
    foreach my $tag_l2 (keys %{$hash_omniscient->{'level2'}}){
      if (exists_keys($hash_omniscient,('level2',$tag_l2,$id_l1) ) ){
        foreach my $feature_l2 (@{$hash_omniscient->{'level2'}{$tag_l2}{$id_l1}}){
          my $level2_ID = lc($feature_l2->_tag_value('ID'));

          if ( exists_keys($hash_omniscient,('level3','cds',$level2_ID)) ){
            my $counterL3=-1;
            my $indexLast = $#{$hash_omniscient->{'level3'}{'cds'}{$level2_ID}};
            my @sortedList = sort {$a->start <=> $b->start} @{$hash_omniscient->{'level3'}{'cds'}{$level2_ID}};
            foreach my $feature_l3 ( @sortedList ){
              #count number feature of tag_l3 type
              $counterL3++;
              #Manage Introns## from the second intron to the last (from index 1 to last index of the table sortedList) ## We go inside this loop only if we have more than 1 feature.
              if($counterL3 > 0 and $counterL3 <= $indexLast){
                my $intronSize = $sortedList[$counterL3]->start - $sortedList[$counterL3-1]->end;
                $shortest_intron = $intronSize if($intronSize < $shortest_intron)
              }
            }
          }
        }
      }
    }
    print "Shortest intron for $id_l1:".$shortest_intron."\n" if($shortest_intron != 10000000000 and $verbose);
    if ($shortest_intron < $Xsize){
      print "flag the gene $id_l1\n";
      $nb_cases++;

      my $feature_l1 = $hash_omniscient->{'level1'}{$tag_l1}{$id_l1};
      $feature_l1->add_tag_value($tag, $shortest_intron);
      if($feature_l1->has_tag('product') ){
        $feature_l1->add_tag_value('note', $feature_l1->get_tag_values('product'));
        $feature_l1->remove_tag('product');
      }
      foreach my $tag_l2 (keys %{$hash_omniscient->{'level2'}}){
        if (exists_keys ($hash_omniscient, ('level2', $tag_l2, $id_l1) ) ) {
          foreach my $feature_l2 (@{$hash_omniscient->{'level2'}{$tag_l2}{$id_l1}}){
            my $level2_ID = lc($feature_l2->_tag_value('ID'));
            $feature_l2->add_tag_value($tag, $shortest_intron);
            if($feature_l2->has_tag('product') ){
              $feature_l2->add_tag_value('note', $feature_l2->get_tag_values('product'));
              $feature_l2->remove_tag('product');
            }

            foreach my $tag_l3 (keys %{$hash_omniscient->{'level3'}}){
              if ( exists_keys($hash_omniscient, ('level3', $tag_l3, $level2_ID) ) ){
                foreach my $feature_l3 (@{$hash_omniscient->{'level3'}{$tag_l3}{$level2_ID}}){
                  $feature_l3->add_tag_value($tag, $shortest_intron);
                  if($feature_l3->has_tag('product') ){
                    $feature_l3->add_tag_value('note', $feature_l3->get_tag_values('product'));
                    $feature_l3->remove_tag('product');
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

my $toprint = "We found $nb_cases cases where introns were < $Xsize, we flagged them with the attribute $tag. The value of this tag is size of the shortest intron found in this gene.\n";
print $ostreamReport $toprint;
if($opt_output){print $toprint;}

print_omniscient( {omniscient => $hash_omniscient, output => $gffout} );

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

agat_sp_flag_short_introns_ebi.pl

=head1 DESCRIPTION

By default, the script flags short introns within coding sequences (CDS) with the <pseudo> attribute and change the <product> attribute into a <note> attribute.
This is useful for avoiding ERROR messages when submitting data to the EBI.
(Typical EBI error message: ERROR: Intron usually expected to be at least 10 nt long. Please check the accuracy.)

=head1 SYNOPSIS

    agat_sp_flag_short_introns_ebi.pl --gff infile --out outfile
    agat_sp_flag_short_introns_ebi.pl --help

=head1 OPTIONS

=over 8

=item B<--gff>, B<-f>, B<--ref> or B<-reffile>

Input GTF/GFF file.

=item  B<--intron_size> or B<-i>

Minimum intron size, default 10. All genes with an intron < of this size will be
flagged with the pseudo attribute (the value will be the size of the smallest
intron found within the incriminated gene)

=item  B<--out>, B<--output> or B<-o>

Output gff3 file where the result will be printed.

=item B<-v>

Bolean. Verbose for debugging purpose.

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
