#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use POSIX qw(strftime);
use Pod::Usage;
use Bio::Tools::GFF;
use IO::File;
use AGAT::Omniscient;

my $header = get_agat_header();
my $opt_test="=";
my $opt_output= undef;
my $opt_nb = 0;
my $opt_gff = undef;
my $opt_verbose = undef;
my $opt_help;

# OPTION MANAGMENT
my @copyARGV=@ARGV;
if ( !GetOptions( 'f|ref|reffile|gff=s' => \$opt_gff,
                  't|test=s'            => \$opt_test,
                  "nb|number|n=i"       => \$opt_nb,
                  'o|output=s'          => \$opt_output,
                  'v|verbose!'          => \$opt_verbose,
                  'h|help!'             => \$opt_help ) )
{
    pod2usage( { -message => 'Failed to parse command line',
                 -verbose => 1,
                 -exitval => 1 } );
}

if ($opt_help) {
    pod2usage( { -verbose => 99,
                 -exitval => 0,
                 -message => "$header\n" } );
}

if ( ! $opt_gff ){
    pod2usage( {
           -message => "$header\nAt least 1 parameter is mandatory:\n1) Input reference gff file: --gff\n\n",
           -verbose => 0,
           -exitval => 2 } );
}

###############
# Manage Output

## FOR GFF FILE
my $gffout_ok ; my $gffout_notok ; my $ostreamReport ;

if ($opt_output) {
  my ($outfile,$path,$ext) = fileparse($opt_output,qr/\.[^.]*/);

  # set file names
  my $outfile_ok = $path.$outfile.$ext;
  my $outfile_notok = $path.$outfile."_remaining".$ext;
  my $outfile_report = $path.$outfile."_report.txt";

  # check existence
  if(-f $outfile_ok){  print "File $outfile_ok already exist.\n";exit;}
  if(-f $outfile_notok){  print "File $outfile_notok already exist.\n";exit;}
  if(-f $outfile_report){  print "File $outfile_report already exist.\n";exit;}

  # create fh
  open( my $fh, '>', $outfile_ok) or die "Could not open file $outfile_ok $!";
  $gffout_ok = Bio::Tools::GFF->new(-fh => $fh, -gff_version => 3 );
  open( my $fhnotok, '>', $outfile_notok) or die "Could not open file $outfile_notok $!";
  $gffout_notok = Bio::Tools::GFF->new(-fh => $fhnotok, -gff_version => 3 );
  open($ostreamReport, '>', $outfile_report) or die "Could not open file $outfile_report $!";
}
else{
  $gffout_ok = Bio::Tools::GFF->new(-fh => \*STDOUT, -gff_version => 3);
  $ostreamReport = Bio::Tools::GFF->new(-fh => \*STDOUT, -gff_version => 3);
}

#Manage test option
if($opt_test ne "<" and $opt_test ne ">" and $opt_test ne "<=" and $opt_test ne ">=" and $opt_test ne "="){
  print "The test to apply is Wrong: $opt_test.\nWe want something among this list: <,>,<=,>= or =.";exit;
}

# start with some interesting information
my $stringPrint = strftime "%m/%d/%Y at %Hh%Mm%Ss", localtime;
$stringPrint .= "\nusage: $0 @copyARGV\n";
$stringPrint .= "We will select genes that contain $opt_test $opt_nb introns.\n";

if ($opt_output){
  print $ostreamReport $stringPrint;
  print $stringPrint;
}
else{ print $stringPrint; }
                          #######################
# >>>>>>>>>>>>>>>>>>>>>>>>#        MAIN         #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                          #######################

######################
### Parse GFF input #
my ($hash_omniscient, $hash_mRNAGeneLink) =  slurp_gff3_file_JD({ input => $opt_gff,
                                                                  verbose => $opt_verbose
                                                                });
print("Parsing Finished\n");
### END Parse GFF input #
#########################
# sort by seq id
my $hash_sortBySeq = gather_and_sort_l1_by_seq_id($hash_omniscient);

my @listok;
my @list2;
#################
# == LEVEL 1 == #
#################
foreach my $seqid (sort { (($a =~ /(\d+)$/)[0] || 0) <=> (($b =~ /(\d+)$/)[0] || 0) } keys %{$hash_sortBySeq}){ # loop over all the feature level1

	foreach my $tag_l1 (sort {$a cmp $b} keys %{$hash_omniscient->{'level1'}}){
		foreach my $feature_l1 ( @{$hash_sortBySeq->{$seqid}{$tag_l1}} ){
			my $id_l1 = lc($feature_l1->_tag_value('ID'));
      my $success=undef;

	    #################
	    # == LEVEL 2 == #
	    #################
	    foreach my $tag_l2 (sort keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...

	      if ( exists_keys( $hash_omniscient, ('level2', $tag_l2, $id_l1) ) ){
	        my @list_fl2 = @{$hash_omniscient->{'level2'}{$tag_l2}{$id_l1}};
	        foreach my $feature_l2 ( @list_fl2 ) {

	          #################
	          # == LEVEL 3 == #
	          #################
	          my $id_l2 = lc($feature_l2->_tag_value('ID'));

	          if ( exists_keys( $hash_omniscient, ('level3', 'exon', $id_l2) ) ){
	            my $nb_exon = @{$hash_omniscient->{'level3'}{'exon'}{$id_l2}};
              my $nb_intron = $nb_exon-1;

              if( test_size( $nb_intron, $opt_test, $opt_nb ) ){
                push @listok, $id_l1;
                $success = 1;
                last;
              }
	          }
	        }
          if( $success ) {
            last ;
          }
          else{
            push @list2, $id_l1;
          }
	      }
	    }
    }
  }
}

# print ok
my $hash_ok = subsample_omniscient_from_level1_id_list($hash_omniscient, \@listok);
print_omniscient($hash_ok, $gffout_ok); #print gene modified in file
%{$hash_ok} = ();
# print remaining if an output is provided
if($opt_output){
  my $hash_remaining = subsample_omniscient_from_level1_id_list($hash_omniscient, \@list2);
  print_omniscient($hash_remaining, $gffout_notok); #print gene modified in file
  %{$hash_remaining} = ();
}

my $test_success = scalar @listok;
my $test_fail = scalar @list2;

$stringPrint = "$test_success genes selected with at least one RNA with $opt_test $opt_nb intron(s).\n";
$stringPrint .= "$test_fail remaining genes that not pass the test.\n";
if ($opt_output){
  print $ostreamReport $stringPrint;
  print $stringPrint;
} else{ print $stringPrint; }

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

sub test_size{
  my ($size, $operator, $nb_ref) = @_;

  if ($operator eq ">"){
    if ($size > $nb_ref){
      return "true";
    }
  }
  if ($operator eq "<"){
    if ($size < $nb_ref){
      return "true";
    }
  }
  if ($operator eq "=" or $operator eq "=="){
    if ($size == $nb_ref){
      return "true";
    }
  }
  if ($operator eq "<="){
    if ($size <= $nb_ref){
      return "true";
    }
  }
  if ($operator eq ">="){
    if ($size >= $nb_ref){
      return "true";
    }
  }
  return undef;
}

__END__

=head1 NAME

agat_sp_filter_gene_by_intron_numbers.pl

=head1 DESCRIPTION

The script aims to filter genes by intron numbers.
It will create two files. one with the genes passing the intron number filter,
the other one with the remaining genes.

Some examples:
Select intronless genes:
agat_sp_filter_gene_by_intron_numbers.pl --gff infile.gff -o result.gff
Select genes with more or equal 10 introns:
agat_sp_filter_gene_by_intron_numbers.pl --gff infile.gff --test ">=" --nb 10 [ --output outfile ]

=head1 SYNOPSIS

    agat_sp_filter_gene_by_intron_numbers.pl --gff infile.gff --test ">=" --nb 10 [ --output outfile ]
    agat_sp_filter_gene_by_intron_numbers.pl --help

=head1 OPTIONS

=over 8

=item B<-f>, B<--reffile>, B<--gff>  or B<-ref>

Input GFF3 file that will be read

=item B<-n>,  B<--nb> or B<--number>

Integer - Number of introns [Default 0]

=item B<-t> or B<--test>
Test to apply (>, <, =, >= or <=). If you use one of these two characters >, <,
please do not forget to quote your parameter like that "<=". Else your terminal will complain.
[Default "="]

=item B<-o> or B<--output>

Output GFF file.  If no output file is specified, the output will be
written to STDOUT.

=item B<-v>

Verbose option for debugging purpose.

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
