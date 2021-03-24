#!/usr/bin/env perl

use strict;
use warnings;
use Sort::Naturally;
use Getopt::Long;
use File::Basename;
use POSIX qw(strftime);
use Pod::Usage;
use Bio::Tools::GFF;
use IO::File;
use AGAT::Omniscient;

my $header = get_agat_header();
my $opt_output = undef;
my $opt_coordinates = undef ;
my $opt_exclude_ov = undef ;
my $opt_gff = undef;
my $opt_verbose = undef;
my $opt_help;

# OPTION MANAGMENT
my @copyARGV=@ARGV;
if ( !GetOptions( 'i|input|gtf|gff=s'         => \$opt_gff,
                  "c|coordinates|tsv|r|ranges=s"=> \$opt_coordinates,
                  "e|exclude!"                => \$opt_exclude_ov,
                  'o|output=s'                => \$opt_output,
                  'v|verbose!'                => \$opt_verbose,
                  'h|help!'                   => \$opt_help ) )
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
           -message => "$header\nAt least 2 parameters are mandatory:\n1) Input reference gff file: --gff\n".
           "2) A file containing one range per line (see help for syntax): --coordinates\n\n",
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

# Manage ranges
my %range_hash;
open my $in_range, "<:encoding(utf8)", $opt_coordinates or die "$opt_coordinates: $!";
my $cpt_line=0;
while (my $line = <$in_range>) {
    $cpt_line++;
    chomp $line;
    $line =~ s/^\s+|\s+$//g; #removing leading and trailing white spaces
    my @array = split /\s/, $line;
    my $size_array = scalar @array;
    if ( $size_array >= 3){
      push @{$range_hash{lc($array[0])}}, [$array[1], $array[2]];
    }
    else{
      print "skip line $cpt_line (At least 3 values expected, only $size_array available): $line\n";
    }
}

my $nb_ranges = keys %range_hash;
# start with some interesting information
my $stringPrint = strftime "%m/%d/%Y at %Hh%Mm%Ss", localtime;
$stringPrint .= "\nusage: $0 @copyARGV\n";
$stringPrint .= "We will get features that are within the $nb_ranges selected ranges.\n";

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

my @listok;
my @listNotOk;
#################
# == LEVEL 1 == #
#################
# sort by seq id
my ( $hash_sortBySeq, $hash_sortBySeq_std, $hash_sortBySeq_topf ) = collect_l1_info_sorted_by_seqid_and_location($hash_omniscient);

#Read by seqId to sort properly the output by seq ID
foreach my $seqid ( sort { (($a =~ /(\d+)$/)[0] || 0) <=> (($b =~ /(\d+)$/)[0] || 0) } keys %{$hash_sortBySeq}){ # loop over all the feature level1
#################
# == LEVEL 1 == #
#################
  foreach my $locationid ( sort { ncmp ($a, $b) } keys %{$hash_sortBySeq->{$seqid} } ){

      my $tag_l1 = $hash_sortBySeq->{$seqid}{$locationid}{'tag'};
      my $id_l1 = $hash_sortBySeq->{$seqid}{$locationid}{'id'};
      my $feature_l1 = $hash_omniscient->{'level1'}{$tag_l1}{$id_l1};
      if ( test_overlap_with_ranges( $feature_l1, \%range_hash, $opt_exclude_ov ) ){
        push @listok, $id_l1;
      }
      else{
        push @listNotOk, $id_l1;
      }

  }
}

# print ok
my $hash_ok = subsample_omniscient_from_level1_id_list($hash_omniscient, \@listok);
print_omniscient($hash_ok, $gffout_ok); #print gene modified in file
%{$hash_ok} = ();
# print remaining if an output is provided
if($opt_output){
  my $hash_remaining = subsample_omniscient_from_level1_id_list($hash_omniscient, \@listNotOk);
  print_omniscient($hash_remaining, $gffout_notok); #print gene modified in file
  %{$hash_remaining} = ();
}

my $test_success = scalar @listok;
my $test_fail = scalar @listNotOk;

$stringPrint = "$test_success record(s) selected within the range(s).\n";
$stringPrint .= "$test_fail record(s) out of the range(s).\n";
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

sub test_overlap_with_ranges{
  my ($feature_l1, $range_hash, $opt_exclude_ov) = @_;
  my $result = undef;
  my $start = $feature_l1->start();
  my $end = $feature_l1->end();
  if (exists_keys($range_hash,( lc($feature_l1->seq_id) ) ) ){
    foreach my $range ( @{$range_hash{lc($feature_l1->seq_id)}} ){
      if(! $opt_exclude_ov){
        if(_overlap($range, [$start,$end])){
          print "feature [".$feature_l1->primary_tag." $start,$end] is included or overlap the range [@$range]\n" if $opt_verbose;
          return 1;
        }
      }
      else{
        if(_include($range, [$start,$end])){
          print "feature [".$feature_l1->primary_tag." $start,$end] is included in the range [@$range]\n" if $opt_verbose;
          return 1;
        }
      }
    }
  }
  return $result;
}

# feature must be in the range
sub _include{
	my($location1, $location2)=@_;
	my $overlap = undef;

	if (($location1->[0] <= $location2->[0]) and ($location1->[1] >= $location2->[1])){
		$overlap = 1;
	}
	return $overlap;
}

# feature must be in the range or overlaping the range
sub _overlap{
	my($location1, $location2)=@_;
	my $overlap = undef;

	if (($location1->[0] <= $location2->[1]) and ($location1->[1] >= $location2->[0])){
		$overlap = 1;
	}
	return $overlap;
}

__END__

=head1 NAME

agat_sp_filter_record_by_coordinates.pl

=head1 DESCRIPTION

The script aims to filter the records to keep only those contained within coordinates
defined in an input csv file.
A record can be a feature or a set of features with part-of relationships.
By default we keep records overlapping the coordinates. The --exclude parameter
allows to keep only record fully contained within the coordinates.

! With default paramater, an exon out of the coordinates can be kept if the gene
it is part of is overlaping the coordinates.

=head1 SYNOPSIS

    agat_sp_filter_record_by_coordinates.pl --gff infile.gff --tsv coordinates.tsv [ --output outfile ]
    agat_sp_filter_record_by_coordinates.pl --help

=head1 OPTIONS

=over 8

=item B<-i>, B<--input>, B<--gtf>  or B<--gff>

Input GTF/GFF file

=item B<-c>, B<--coordinates>, B<--tsv>, B<-r> or B<--ranges>

String - tsv file containing the coordinates.
Coordinates must be one per line.
Each line must contain 3 fields separated by a tabulation.
Field1 is the sequence id
Field2 is the start coordinate (included)
Field3 is the end coordinate (included)

=item B<-e> or B<--exclude>

Select only the features fully containined within the coordinates, exclude the overlapping
ones.

=item B<-o> or B<--output>

Output GFF file.  If no output file is specified, the output will be
written to STDOUT.

=item B<-v> or B<--verbose>

Verbosity.

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
