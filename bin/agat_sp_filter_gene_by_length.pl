#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use POSIX qw(strftime);
use Pod::Usage;
use IO::File;
use AGAT::AGAT;

my $header = get_agat_header();
my $config;
my $cpu;
my $opt_test=">";
my $opt_output= undef;
my $opt_size = 100;
my $opt_gff = undef;
my $opt_verbose = undef;
my $opt_help;

# OPTION MANAGMENT
my @copyARGV=@ARGV;
if ( !GetOptions( 'f|ref|reffile|gff=s' => \$opt_gff,
                  't|test=s'            => \$opt_test,
                  "s|size=i"            => \$opt_size,
                  'o|output=s'          => \$opt_output,
                  'v|verbose!'          => \$opt_verbose,
                  'c|config=s'          => \$config,
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

# --- Manage config ---
initialize_agat({ config_file_in => $config });
$CONFIG->{cpu} = $cpu if defined($cpu);

###############
# Manage Output

## FOR GFF FILE
my $gffout_ok_file ;
my $gffout_notok_file ;
my $ostreamReport_file ;

if ($opt_output) {
  my ($outfile,$path,$ext) = fileparse($opt_output,qr/\.[^.]*/);

  # set file names
  $gffout_ok_file = $path.$outfile.$ext;
  $gffout_notok_file = $path.$outfile."_remaining".$ext;
  $ostreamReport_file = $path.$outfile."_report.txt";
}

my $gffout_ok = prepare_gffout( $gffout_ok_file );
my $gffout_notok = prepare_gffout( $gffout_notok_file );
my $ostreamReport = prepare_fileout( $ostreamReport_file );

#Manage test option
if($opt_test ne "<" and $opt_test ne ">" and $opt_test ne "<=" and $opt_test ne ">=" and $opt_test ne "="){
  print "The test to apply is Wrong: $opt_test.\nWe want something among this list: <,>,<=,>= or =.";exit;
}

# start with some interesting information
my $stringPrint = strftime "%m/%d/%Y at %Hh%Mm%Ss", localtime;
$stringPrint .= "\nusage: $0 @copyARGV\n";
$stringPrint .= "We will select l1 feature (e.g. gene) that have length $opt_test $opt_size bp.\n";

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
my ($hash_omniscient, $hash_mRNAGeneLink) =  slurp_gff3_file_JD({ input => $opt_gff });
### END Parse GFF input #
#########################
# sort by seq id
my $hash_sortBySeq = gather_and_sort_l1_by_seq_id($hash_omniscient);

my @listok;
my @listNotOk;
#################
# == LEVEL 1 == #
#################
foreach my $seqid (sort { (($a =~ /(\d+)$/)[0] || 0) <=> (($b =~ /(\d+)$/)[0] || 0) } keys %{$hash_sortBySeq}){ # loop over all the feature level1

	foreach my $tag_l1 (sort {$a cmp $b} keys %{$hash_omniscient->{'level1'}}){
		foreach my $feature_l1 ( @{$hash_sortBySeq->{$seqid}{$tag_l1}} ){
			my $id_l1 = lc($feature_l1->_tag_value('ID'));
      my $gene_length=$feature_l1->end()-$feature_l1->start()+1;
      my $successl1 = test_size( $gene_length, $opt_test, $opt_size );
      my $longer_concat_exon=undef;
	    #################
	    # == LEVEL 2 == #
	    #################
	    foreach my $tag_l2 (sort keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
	      if ( exists_keys( $hash_omniscient, ('level2', $tag_l2, $id_l1) ) ){
	        foreach my $feature_l2 ( @{$hash_omniscient->{'level2'}{$tag_l2}{$id_l1}} ) {
            my $id_l2 = lc($feature_l2->_tag_value('ID'));
	          #################
	          # == LEVEL 3 == #
	          #################
            foreach my $tag_l3 (sort keys %{$hash_omniscient->{'level3'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
  	          if ( exists_keys( $hash_omniscient, ('level3', 'exon', $id_l2) ) ){
                my $local_size=0;
                foreach my $feature_l3 ( @{$hash_omniscient->{'level3'}{'exon'}{$id_l2}} ) {
                  $local_size += $feature_l3->end()-$feature_l3->start()+1;
                }
                if($longer_concat_exon and $longer_concat_exon<$local_size){
                  $longer_concat_exon = $local_size;
                }
                elsif(! $longer_concat_exon) {
                  $longer_concat_exon = $local_size;
                }
              }
	          }
	        }
	      }
	    }
      # case we had exon (we look at the longest mRNA)
      if($longer_concat_exon){
        print "$id_l1 does have exon(s). Longest concatenated exons: $longer_concat_exon\n" if $opt_verbose;
        if( test_size( $longer_concat_exon, $opt_test, $opt_size ) ){
          print "$id_l1 pass the test\n" if $opt_verbose;
          push @listok, $id_l1;
        }
        else{
          print "$id_l1 do not pass the test\n" if $opt_verbose;
          push @listNotOk, $id_l1;
        }
      }
      else{
        print "$id_l1 does not have any exon. $tag_l1 size: $gene_length\n" if $opt_verbose;
        # No exon, L1 pass test
        if($successl1){
          print "$id_l1 pass the test\n" if $opt_verbose;
          push @listok, $id_l1;
        }
        # No exon, L1 do not pass test
        else{
          print "$id_l1 do not pass the test\n" if $opt_verbose;
          push @listNotOk, $id_l1;
        }
      }
    }
  }
}

# print ok
my $hash_ok = subsample_omniscient_from_level1_id_list_delete($hash_omniscient, \@listok);
print_omniscient( {omniscient => $hash_ok, output => $gffout_ok} );
%{$hash_ok} = ();
# print remaining if an output is provided
if($opt_output){
  my $hash_remaining = subsample_omniscient_from_level1_id_list_delete($hash_omniscient, \@listNotOk);
  print_omniscient( {omniscient => $hash_remaining, output => $gffout_notok} );
  %{$hash_remaining} = ();
}

my $test_success = scalar @listok;
my $test_fail = scalar @listNotOk;

$stringPrint = "$test_success l1 feature (e.g. gene) selected with a length $opt_test $opt_size bp.\n";
$stringPrint .= "$test_fail remaining l1 feature (e.g. gene) do not pass the test.\n";
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

agat_sp_filter_gene_by_length.pl

=head1 DESCRIPTION

The script aims to filter level1 feature (e.g. gene, match, etc) by length.
It will create two files. one with the feature passing the length filter,
the other one with the remaining features.
If the level1 feature has exon features, the size is computed by concatenating
the exon together. If the level1 feature has several level2 features (e.g. mRNA)
we apply the test on the longest one (the longest concatenated exon set).

Some examples:
Select L1 feature shorter than 1000bp:
agat_sp_filter_gene_by_length.pl --gff infile.gff  --size 1000 --test "<" -o result.gff
Select genes longer than 200bp:
agat_sp_filter_gene_by_length.pl --gff infile.gff --size 200 --test ">" -o result.gff

=head1 SYNOPSIS

    agat_sp_filter_gene_by_length.pl --gff infile.gff --test ">=" --nb 10 [ --output outfile ]
    agat_sp_filter_gene_by_length.pl --help

=head1 OPTIONS

=over 8

=item B<-f>, B<--reffile>, B<--gff>  or B<-ref>

Input GFF3 file that will be read

=item B<-s> or B<--size>

Integer - Gene size in pb [Default 100]

=item B<-t> or B<--test>

Test to apply (>, <, =, >= or <=). If you use one of these two characters >, <,
please do not forget to quote your parameter like that "<=". Else your terminal will complain.
[Default "="]

=item B<-o> or B<--output>

Output GFF file.  If no output file is specified, the output will be
written to STDOUT.

=item B<-v>

Verbose option for debugging purpose.

=item B<-c> or B<--config>

String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any, 
otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

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
