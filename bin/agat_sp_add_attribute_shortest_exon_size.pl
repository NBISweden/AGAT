#!/usr/bin/env perl

use strict;
use warnings;
use POSIX qw(strftime);
use File::Basename;
use Carp;
use IO::File;
use Getopt::Long::Descriptive;
use Pod::Usage;
use AGAT::AGAT;

my $header = get_agat_header();
my @copyARGV = @ARGV;
my ( $opt, $usage, $config ) = AGAT::AGAT::describe_script_options( $header,
    [ 'gff|f|ref=s', 'Input GTF/GFF file', { required => 1 } ],
);

my $opt_file   = $opt->gff;
my $opt_output = $opt->out;

pod2usage( { -verbose => 99, -exitstatus => 0, -message => "$header\n" } ) if $opt->help;


# #######################
# # START Manage Option #
# #######################
my $ostreamReport_filename;
if (defined($opt_output) ) {
  my ($filename,$path,$ext) = fileparse($opt_output,qr/\.[^.]*/);
  $ostreamReport_filename=$path.$filename."_report.txt";
}
my $log;
if ( my $log_name = $config->{log_path} ) {
  open( $log, '>', $log_name ) or die "Can not open $log_name for printing: $!";
  dual_print( $log, $header,  3 );
}

my $gffout = prepare_gffout($config, $opt_output);
my $ostreamReport = prepare_fileout($ostreamReport_filename);

my $string1 = strftime "%m/%d/%Y at %Hh%Mm%Ss", localtime;
$string1 .= "\n\nusage: $0 @copyARGV\n\n";

print $ostreamReport $string1;
dual_print( $log, $string1) if $opt_output;

                                                      #######################
                                                      #        MAIN         #
#                     >>>>>>>>>>>>>>>>>>>>>>>>>       #######################       <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

######################
### Parse GFF input #
dual_print( $log, "Reading $opt_file\n");
my ($hash_omniscient, $hash_mRNAGeneLink) = slurp_gff3_file_JD({ input => $opt_file,
                                                                 config => $config });
dual_print( $log, "Parsing Finished\n\n");
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
dual_print( $log, $toprint) if $opt_output;
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

=item B<--out>, B<--output> or B<-o>

STRING: Output gff3 file where the result will be printed.

=item B<-v> or B<--verbose>

BOLEAN: Verbose for debugging purpose.

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
