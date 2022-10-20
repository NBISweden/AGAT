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
my $config = get_agat_config();
my $opt_file;
my $opt_output=undef;
my $verbose=undef;
my $opt_help = 0;

my @copyARGV=@ARGV;
if ( !GetOptions( 'f|gff|ref=s' => \$opt_file,
                  'o|out|output=s' => \$opt_output,
                  'v|verbose!'      => \$verbose,
                  'h|help!'         => \$opt_help ) )
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

# #######################
# # START Manage Option #
# #######################
my $ostreamReport_filename;
if (defined($opt_output) ) {
  my ($filename,$path,$ext) = fileparse($opt_output,qr/\.[^.]*/);
  $ostreamReport_filename=$path.$filename."_report.txt";
}
my $gffout = prepare_gffout($config, $opt_output);
my $ostreamReport = prepare_fileout($ostreamReport_filename);

my $string1 = strftime "%m/%d/%Y at %Hh%Mm%Ss", localtime;
$string1 .= "\n\nusage: $0 @copyARGV\n\n";

print $ostreamReport $string1;
if($opt_output){print $string1;}

                                                      #######################
                                                      #        MAIN         #
#                     >>>>>>>>>>>>>>>>>>>>>>>>>       #######################       <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

######################
### Parse GFF input #
print "Reading ".$opt_file,"\n";
my ($hash_omniscient, $hash_mRNAGeneLink) = slurp_gff3_file_JD({ input => $opt_file,
                                                                 config => $config });
print("Parsing Finished\n\n");
### END Parse GFF input #
#########################

my $nb_cases_l1=0;
my $nb_cases_l2=0;
my $tag = "shortest_exon";
######################
### Parse GFF input #
foreach my $tag_l1 (keys %{$hash_omniscient->{'level1'}}){
  foreach my $id_l1 (keys %{$hash_omniscient->{'level1'}{$tag_l1}}){
    my $feature_l1 = $hash_omniscient->{'level1'}{$tag_l1}{$id_l1};
    my $shortest_exon=undef;
    foreach my $tag_l2 (keys %{$hash_omniscient->{'level2'}}){
      if (exists_keys($hash_omniscient,('level2',$tag_l2,$id_l1) ) ){
        my $shortest_exon_l2=undef;
        foreach my $feature_l2 (@{$hash_omniscient->{'level2'}{$tag_l2}{$id_l1}}){
          my $level2_ID = lc($feature_l2->_tag_value('ID'));

          if ( exists_keys($hash_omniscient,('level3','exon',$level2_ID)) ){

            foreach my $feature_l3 (  @{$hash_omniscient->{'level3'}{'exon'}{$level2_ID}} ){
              my $exonSize = $feature_l3->end - $feature_l3->start+1;
              if(! $shortest_exon_l2){
                $shortest_exon_l2 = $exonSize;
              }
              elsif( $exonSize < $shortest_exon_l2){
                $shortest_exon_l2 = $exonSize;
              }
            }
            $feature_l2->add_tag_value($tag, $shortest_exon_l2);
            $nb_cases_l2++;
            if(! $shortest_exon){
              $shortest_exon = $shortest_exon_l2;
            }
            elsif( $shortest_exon_l2 < $shortest_exon){
              $shortest_exon = $shortest_exon_l2;
            }
          }
        }
      }
    }
    if($shortest_exon){
      $feature_l1->add_tag_value($tag, $shortest_exon);
      $nb_cases_l1++;
    }
  }
}

my $toprint = "$nb_cases_l1 $tag flags/attributes added to level1 features and $nb_cases_l2 $tag flags/attributes added to level2 features. The value of the attribute is size of the shortest exon found.\n";
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

agat_sp_add_attribute_shortest_exon_size.pl

=head1 DESCRIPTION

The script add the attribute <shortest_exon> to each gene and rna, which will hold the size of the shortest exon in bp.

=head1 SYNOPSIS

    agat_sp_add_attribute_shortest_exon_size.pl --gff infile --out outfile
    agat_sp_add_attribute_shortest_exon_size.pl --help

=head1 OPTIONS

=over 8

=item B<--gff>, B<-f> or B<--ref>

STRING: Input GTF/GFF file.

=item  B<--out>, B<--output> or B<-o>

STRING: Output gff3 file where the result will be printed.

=item B<-v>

BOLEAN: Verbose for debugging purpose.

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
