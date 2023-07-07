#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use POSIX qw(strftime);
use Getopt::Long;
use Pod::Usage;
use AGAT::AGAT;

my $start_run = time();
my $header = get_agat_header();
my $config;
my $PROT_LENGTH = 100;
my $file_fasta=undef;
my $outfile = undef;
my $verbose = undef;
my $opt_test = undef;
my $gff = undef;
my $opt_help= 0;

my @copyARGV=@ARGV;
Getopt::Long::Configure ('bundling');
if ( !GetOptions(
    "help|h"   => \$opt_help,
    "g|gff=s"  => \$gff,
    't|test=s' => \$opt_test,
    "size|s=i" => \$PROT_LENGTH,
    "v!"       => \$verbose,
    "output|outfile|out|o=s" => \$outfile))

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

if ( ! (defined($gff)) ){
    pod2usage( {
           -message => "$header\nAt least 1 parameter is mandatory:\n Input reference gff file (--gff)\n\n",
           -verbose => 0,
           -exitval => 1 } );
}

######################
# Option check

if($opt_test){
  if($opt_test ne "<" and $opt_test ne ">" and $opt_test ne "<=" and $opt_test ne ">=" and $opt_test ne "=" and $opt_test ne "=="){
    print "The test to apply is Wrong: $opt_test.\nWe want something among this list: <,>,<=,>=,== or =.";exit;
  }
}
else{
  $opt_test = ">";
}

# To avoid > < character in output files
my $opt_test_to_print = $opt_test;
$opt_test_to_print =~ s/>/sup/ig;
$opt_test_to_print =~ s/</inf/ig;

######################
# Manage output file #
my $gffout_pass_file;
my $gffout_notpass_file;
if ($outfile) {
  $outfile=~ s/.gff//g;
  $gffout_pass_file = $outfile."_".$opt_test_to_print.$PROT_LENGTH.".gff";
  $gffout_notpass_file = $outfile."_NOT_".$opt_test_to_print.$PROT_LENGTH.".gff";
}

my $gffout_pass = prepare_gffout($config, $gffout_pass_file);
my $gffout_notpass = prepare_gffout($config, $gffout_notpass_file);

# print usage performed
my $stringPrint = strftime "%m/%d/%Y at %Hh%Mm%Ss", localtime;
$stringPrint = "Launched the ".$stringPrint."\nusage: $0 @copyARGV\n";
$stringPrint .= "We are filtering the gene with protein size $opt_test $PROT_LENGTH\n";
print $stringPrint;

                #####################
                #     MAIN          #
                #####################

######################
### Parse GFF input #
my ($hash_omniscient, $hash_mRNAGeneLink) = slurp_gff3_file_JD({ input => $gff,
                                                                 config => $config
                                                               });
print ("GFF3 file parsed\n");

my @good_gene_list;
my @bad_gene_list;
my $number_pass=0;
foreach my $primary_tag_l1 (keys %{$hash_omniscient->{'level1'}}){ # primary_tag_l1 = gene or repeat etc...
  foreach my $gene_id_l1 (keys %{$hash_omniscient->{'level1'}{$primary_tag_l1}}){
    my $gene_feature=$hash_omniscient->{'level1'}{$primary_tag_l1}{$gene_id_l1};
    print "Study gene $gene_id_l1\n" if($verbose);
		my $no_l2=1;# see if standalone or topfeature

    foreach my $primary_tag_l2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
      my $there_is_cds=undef;
      my $one_pass=undef;
      if ( exists_keys( $hash_omniscient, ('level2', $primary_tag_l2, $gene_id_l1) ) ){
				$no_l2 = undef;
        foreach my $level2_feature ( @{$hash_omniscient->{'level2'}{$primary_tag_l2}{$gene_id_l1}}) {

          # get level2 id
          my $id_level2 = lc($level2_feature->_tag_value('ID'));

          ##############################
          #If it's a mRNA = have CDS. #
					if( exists_keys( $hash_omniscient, ('level3', 'cds', $id_level2 ) ) ) {
            $there_is_cds="true";

            ##############
            # Manage CDS #
            my $cds_size=0;
            foreach my $cds (@{$hash_omniscient->{'level3'}{'cds'}{$id_level2}}){
              $cds_size+= ($cds->end() - $cds->start() + 1);
            }
            $cds_size = ($cds_size - 3) / 3; # Remove the stop codon and divide by 3 to get Amnino acid

            if(test_size($cds_size, $PROT_LENGTH, $opt_test) ){
              $one_pass="true";
            }
          }
        }
        if($there_is_cds){
          if($one_pass){
            push(@good_gene_list, $gene_id_l1);
            $number_pass++;
          }
          else{
            push(@bad_gene_list, $gene_id_l1);
          }
        }
        else{
          print "No cds for $gene_id_l1\n" if ($verbose);
          push(@good_gene_list, $gene_id_l1);
        }
      }
    }
		if($no_l2){ # case of l1 feature without child
			print "No child for $gene_id_l1\n" if ($verbose);
			push(@good_gene_list, $gene_id_l1);
		}
  }
}

#resume

my $number_notpass=$#bad_gene_list+1;
print_omniscient_from_level1_id_list( {omniscient => $hash_omniscient, level_id_list =>\@good_gene_list, output => $gffout_pass} );
print_omniscient_from_level1_id_list( {omniscient => $hash_omniscient, level_id_list =>\@bad_gene_list, output => $gffout_notpass} );

print "$number_pass genes passed the test.\n";
print "$number_notpass genes didn't pass the test.\n";

# END
my $end_run = time();
my $run_time = $end_run - $start_run;
print "Job done in $run_time seconds\n";
#######################################################################################################################
        ####################
         #     METHODS    #
          ################
           ##############
            ############
             ##########
              ########
               ######
                ####
                 ##

sub test_size{
  my ($size, $PROT_LENGTH, $operator) = @_;

  if ($operator eq ">"){
    if ($size > $PROT_LENGTH){
      return "true";
    }
  }
  if ($operator eq "<"){
    if ($size < $PROT_LENGTH){
      return "true";
    }
  }
  if ($operator eq "=" or $operator eq "=="){
    if ($size == $PROT_LENGTH){
      return "true";
    }
  }
  if ($operator eq "<="){
    if ($size <= $PROT_LENGTH){
       return "true";
    }
  }
  if ($operator eq ">="){
    if ($size >= $PROT_LENGTH){
      return "true";
    }
  }
return undef;
}

__END__

=head1 NAME

agat_sp_filter_by_ORF_size.pl

=head1 DESCRIPTION

The script reads a gff annotation file, and create two output files,
one contains the gene models with ORF passing the test, the other contains the rest.
By default the test is "> 100" that means all gene models that have ORF longer
than 100 Amino acids, will pass the test.

=head1 SYNOPSIS

    agat_sp_filter_by_ORF_size.pl --gff infile.gff [ -o outfile ]
    agat_sp_filter_by_ORF_size.pl -h

=head1 OPTIONS

=over 8

=item B<-g> or B<--gff>

Input GTF/GFF file.

=item B<-s> or B<--size>

ORF size to apply the test. Default 100.

=item B<-t> or B<--test>
Test to apply (> < = >= <=). If you us one of these two character >, <, please don't forget to quote you parameter liket that "<=". Else your terminal will complain.
By default it will be ">"

=item B<-v>

Verbose. Useful for debugging purpose. Bolean

=item B<-o> or B<--out> or B<--output> or B<--outfile>

Output GFF file.  If no output file is specified, the output will be
written to STDOUT.

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
