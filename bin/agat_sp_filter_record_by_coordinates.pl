#!/usr/bin/env perl

use strict;
use warnings;
use Sort::Naturally;
use Getopt::Long;
use File::Basename;
use POSIX qw(strftime);
use Pod::Usage;
use IO::File;
use AGAT::AGAT;

my $header = get_agat_header();
my $config;
my $cpu;
my $opt_output ;
my $opt_coordinates ;
my $opt_exclude_ov ;
my $opt_gff ;
my $opt_verbose ;
my $opt_help ;

# OPTION MANAGMENT
my @copyARGV=@ARGV;
if ( !GetOptions( 'i|input|gtf|gff=s'            => \$opt_gff,
                  "coordinates|tsv|r|ranges=s"   => \$opt_coordinates,
                  "e|exclude!"                   => \$opt_exclude_ov,
                  'o|output=s'                   => \$opt_output,
                  'v|verbose!'                   => \$opt_verbose,
                  'c|config=s'                   => \$config,
                    'thread|threads|cpu|cpus|core|cores|job|jobs=i' => \$cpu,
                  'h|help!'                      => \$opt_help ) )
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

# --- Manage config ---
initialize_agat({ config_file_in => $config, input => $opt_gff });
$CONFIG->{cpu} = $cpu if defined($cpu);

###############
# Manage Output

if (! $opt_output) {
  print "Default output name: filter_record_by_coordinates\n";
  $opt_output="filter_record_by_coordinates";
}

if (-d $opt_output){
  print "The output directory choosen already exists. Please give me another Name.\n";exit();
}
mkdir $opt_output;

## FOR GFF FILE
my $gffout_ok_file ;
my $gffout_notok_file ;
my $ostreamReport_file ;

if ($opt_output) {

  #my $outfile_ok = $path.$outfile.$ext;
  $gffout_notok_file = $opt_output."/remaining.gff3";
  $ostreamReport_file = $opt_output."/report.txt";
}

my $gffout_notok = prepare_gffout( $gffout_notok_file );
my $ostreamReport =  prepare_fileout( $ostreamReport_file );


# Manage ranges
my %range_hash;
open my $in_range, "<:encoding(utf8)", $opt_coordinates or die "$opt_coordinates: $!";
my $cpt_line=0;
my $nb_ranges = 0;
while (my $line = <$in_range>) {
    $cpt_line++;
    chomp $line;
    $line =~ s/^\s+|\s+$//g; #removing leading and trailing white spaces
    my @array = split /\s/, $line;
    my $size_array = scalar @array;
    if ( $size_array >= 3){
      push @{$range_hash{lc($array[0])}}, [$array[1], $array[2]];
      $nb_ranges++;
    }
    else{
      print "skip line $cpt_line (At least 3 values expected, only $size_array available): $line\n";
    }
}

# start with some interesting information
my $stringPrint = strftime "%m/%d/%Y at %Hh%Mm%Ss", localtime;
$stringPrint .= "\nusage: $0 @copyARGV\n";
$stringPrint .= "We will get features that are within the $nb_ranges selected ranges.\n";

if ($opt_output){
  print $ostreamReport $stringPrint;
  print $stringPrint;
}
else{ print $stringPrint; }

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>     MAIN     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

######################
### Parse GFF input #
my ($hash_omniscient) =  slurp_gff3_file_JD({ input => $opt_gff });
### END Parse GFF input #
#########################

my %hash_listok;
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

    my $list_ranges = test_overlap_with_ranges( $feature_l1, \%range_hash, $opt_exclude_ov );
    if ( @$list_ranges ){
      foreach my $range (@$list_ranges){
        push @{$hash_listok{$range}}, $id_l1;
      }
    }
    else{
      push @listNotOk, $id_l1;
    }
  }
}


# print ok
my $test_success=0;
foreach my $range ( sort { ncmp ($a, $b) } keys %hash_listok ){
  my $listok = $hash_listok{$range};
  $test_success += scalar @{ $listok };
  my $hash_ok = subsample_omniscient_from_level1_id_list_intact($hash_omniscient, $listok);

	$gffout_ok_file = "$opt_output/$range.gff3";
	my $gffout_ok = prepare_gffout( $gffout_ok_file);

  print_omniscient( {omniscient => $hash_ok, output => $gffout_ok} );
  %{$hash_ok} = (); #clean
}

# print remaining if an output is provided
if($opt_output){
  my $hash_remaining = subsample_omniscient_from_level1_id_list_intact($hash_omniscient, \@listNotOk);
  print_omniscient( {omniscient => $hash_remaining, output => $gffout_notok} );
  %{$hash_remaining} = ();
}
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

  my @list_ranges = (); #  An empty array.  @names = undef will leave the array with a single element which is undef.
  my $start = $feature_l1->start();
  my $end = $feature_l1->end();
  if (exists_keys($range_hash,( lc($feature_l1->seq_id) ) ) ){
    foreach my $range ( @{$range_hash{lc($feature_l1->seq_id)}} ){
      if(! $opt_exclude_ov){
        if(_overlap($range, [$start,$end])){
          print "feature [".$feature_l1->primary_tag." $start,$end] is included or overlap the range [@$range]\n" if $opt_verbose;
          my $range_string = $feature_l1->seq_id."_".$range->[0]."_".$range->[1];
          push @list_ranges, $range_string;
        }
      }
      else{
        if(_include($range, [$start,$end])){
          print "feature [".$feature_l1->primary_tag." $start,$end] is included in the range [@$range]\n" if $opt_verbose;
          my $range_string = $feature_l1->seq_id."_".$range->[0]."_".$range->[1];
          push @list_ranges, $range_string;
        }
      }
    }
  }
  return \@list_ranges;
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

=item B<--coordinates>, B<--tsv>, B<-r> or B<--ranges>

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

Output folder.

=item B<-v> or B<--verbose>

Verbosity.

=item B<-thread>, B<threads>, B<cpu>, B<cpus>, B<core>, B<cores>, B<job> or B<jobs>

Integer - Number of parallel processes to use for file input parsing (via forking).

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
