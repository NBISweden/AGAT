#!/usr/bin/env perl

use strict;
use warnings;
use Sort::Naturally;
use Getopt::Long;
use File::Basename;
use Pod::Usage;
use IO::File;
use AGAT::AGAT;

start_script();
my $header = get_agat_header();
# -----------------------------------------------------------------------------------------------
my $opt_output ;
my $opt_coordinates ;
my $opt_exclude_ov ;
my $opt_gff ;
my $opt_help ;

# OPTION MANAGMENT
my ($shared_argv, $script_argv) = split_argv_shared_vs_script(\@ARGV);

# Parse script-specific options
my $script_parser = Getopt::Long::Parser->new;
$script_parser->configure('bundling','no_auto_abbrev');
if ( ! $script_parser->getoptionsfromarray(
                  $script_argv,
                  'i|input|gtf|gff=s'            => \$opt_gff,
                  'coordinates|tsv|r|ranges=s'   => \$opt_coordinates,
                  'e|exclude!'                   => \$opt_exclude_ov,
                  'o|out|output=s'                   => \$opt_output,
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
my $shared_opts = parse_shared_options($shared_argv);
initialize_agat({ config_file_in => ( $shared_opts->{config} ), input => $opt_gff, shared_opts => $shared_opts });
# -----------------------------------------------------------------------------------------------

###############
# Manage Output

if (! $opt_output) {
  dual_print1 "Default output name: filter_record_by_coordinates\n";
  $opt_output="filter_record_by_coordinates";
}

if (-d $opt_output){
  die "The output directory choosen already exists. Please give me another Name.\n";
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
      dual_print1 "skip line $cpt_line (At least 3 values expected, only $size_array available): $line\n";
    }
}

# start with some interesting information
my $stringPrint = "We will get features that are within the $nb_ranges selected ranges.\n";

print $ostreamReport $stringPrint if ($opt_output);
dual_print1 $stringPrint;

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

my $stringPrint2 = "$test_success record(s) selected within the range(s).\n";
$stringPrint2 .= "$test_fail record(s) out of the range(s).\n";

print $ostreamReport $stringPrint2 if ($opt_output);
dual_print1 $stringPrint2;

# ----------------------------- END --------------------------------
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

sub test_overlap_with_ranges{
  my ($feature_l1, $range_hash, $opt_exclude_ov) = @_;

  my @list_ranges = (); #  An empty array.  @names = undef will leave the array with a single element which is undef.
  my $start = $feature_l1->start();
  my $end = $feature_l1->end();
  if (exists_keys($range_hash,( lc($feature_l1->seq_id) ) ) ){
    foreach my $range ( @{$range_hash{lc($feature_l1->seq_id)}} ){
      if(! $opt_exclude_ov){
        if(_overlap($range, [$start,$end])){
          dual_print2 "feature [".$feature_l1->primary_tag." $start,$end] is included or overlap the range [@$range]\n";
          my $range_string = $feature_l1->seq_id."_".$range->[0]."_".$range->[1];
          push @list_ranges, $range_string;
        }
      }
      else{
        if(_include($range, [$start,$end])){
          dual_print2 "feature [".$feature_l1->primary_tag." $start,$end] is included in the range [@$range]\n";
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

=item B<-i>, B<--input>, B<--gtf>  or B<--gff> <file>

Input GTF/GFF file

=item B<--coordinates>, B<--tsv>, B<-r> or B<--ranges> <file>

tsv file containing the coordinates.
Coordinates must be one per line.
Each line must contain 3 fields separated by a tabulation.
Field1 is the sequence id
Field2 is the start coordinate (included)
Field3 is the end coordinate (included)

=item B<-e> or B<--exclude>

Select only the features fully containined within the coordinates, exclude the overlapping
ones.

=item B<-o>, B<--out> or B<--output> <folder>

Output folder.

=item B<-h> or B<--help>

Display this helpful text.

=back

=head1 SHARED OPTIONS

Shared options are defined in the AGAT configuration file and can be overridden via the command line for this script only.
Common shared options are listed below; for the full list, please refer to the AGAT agat_config.yaml.

=over 8

=item B<--config> <file>

Path to a custom AGAT configuration file.  
By default, AGAT uses `agat_config.yaml` from the working directory if present, otherwise the default file shipped with AGAT
(available locally via `agat config --expose`).

=item B<--cpu>, B<--core>, B<--job> or B<--thread> <int>

Number of parallel processes to use for file input parsing (via forking).

=item B<-v> or B<--verbose> <int>

Verbosity, choice are 0,1,2,3,4. 0 is quiet, 1 is normal, 2,3,4 is more verbose. Default 1.

=back

=head1 FEEDBACK

For questions, suggestions, or general discussions about AGAT, please use the AGAT community forum:
https://github.com/NBISweden/AGAT/discussions

=head1 BUG REPORTING

Bug reports should be submitted through the AGAT GitHub issue tracker:
https://github.com/NBISweden/AGAT/issues

=cut

AUTHOR - Jacques Dainat
