#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Pod::Usage;
use IO::File;
use AGAT::AGAT;

start_script();
my $header = get_agat_header();
# -------------------------------- LOAD OPTIONS --------------------------------
my $primaryTag=undef;
my $opt_output= undef;
my $opt_attribute = undef;
my $opt_test = undef;
my $opt_gff = undef;
my $opt_help;

# ---------------------------- OPTIONS ----------------------------
# Partition @ARGV into shared vs script options
my ($shared_argv, $script_argv) = split_argv_shared_vs_script(\@ARGV);

# Parse script-specific options
my $script_parser = Getopt::Long::Parser->new;
$script_parser->configure('bundling','no_auto_abbrev');
if ( ! $script_parser->getoptionsfromarray(
  $script_argv,
  'f|ref|reffile|gff=s' => \$opt_gff,
  'p|type|l=s'          => \$primaryTag,
  'a|att|attribute=s'   => \$opt_attribute,
  'flip!'               => \$opt_test,
  'o|output=s'          => \$opt_output,
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

if ( ! $opt_gff or ! $opt_attribute ){
    pod2usage( {
           -message => "$header\nAt least 2 parameters are mandatory:\n1) Input reference gff file: -f\n".
           "2) An attribute tag: -a\n\n",
           -verbose => 0,
           -exitval => 2 } );
}

# Parse shared options and initialize AGAT
my ($shared_opts) = parse_shared_options($shared_argv);
initialize_agat({ config_file_in => ( $shared_opts->{config} ), input => $opt_gff, shared_opts => $shared_opts });

###############
# Manage Output

## FOR GFF FILE
my $gffout_ok_file ;
my $fhout_discarded_file ;
my $ostreamReport_file;

if ($opt_output) {
  my ($outfile,$path,$ext) = fileparse($opt_output,qr/\.[^.]*/);

  # set file names
  $gffout_ok_file = $path.$outfile.$ext;
  $fhout_discarded_file = $path.$outfile."_discarded.txt";
  $ostreamReport_file = $path.$outfile."_report.txt";
}

my $gffout_ok = prepare_gffout( $gffout_ok_file );
my $fhout_discarded = prepare_fileout($fhout_discarded_file);
my $ostreamReport = prepare_fileout($ostreamReport_file);


# Manage $primaryTag
my @ptagList;
my $print_feature_string;
if(! $primaryTag or $primaryTag eq "all"){
  $print_feature_string = "all features";
  push(@ptagList, "all");
}
elsif($primaryTag =~/^level[123]$/){
  $print_feature_string .= "$primaryTag features ";
  push(@ptagList, $primaryTag);
}
else{
   @ptagList= split(/,/, $primaryTag);
   foreach my $tag (@ptagList){
      if($tag =~/^level[123]$/){
        $print_feature_string .= "$primaryTag features ";
      }
      else{
        $print_feature_string .= "$tag feature ";
      }
   }
}

# Manage attributes if given
### If attributes given, parse them:
my @attListOk;
if ($opt_attribute){
  my @attList = split(/,/, $opt_attribute); # split at comma as separated value

  foreach my $attribute (@attList){
      push @attListOk, $attribute;
      dual_print1 "$attribute attribute will be processed.\n";

  }
  dual_print1 "\n";
}

# start with some interesting information
my $stringPrint = "\nWe will discard $print_feature_string that have the attribute $opt_attribute.\n";

print $ostreamReport $stringPrint if ($opt_output);
dual_print1 $stringPrint;

													#######################
# >>>>>>>>>>>>>>>>>>>>>>>>#        MAIN         #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
													#######################

my %all_cases = ('l1' => 0, 'l2' => 0, 'l3' => 0, 'all' => 0);
######################
### Parse GFF input #
my ($hash_omniscient) =  slurp_gff3_file_JD({ input => $opt_gff });;
### END Parse GFF input #
#########################
# sort by seq id
my $hash_sortBySeq = gather_and_sort_l1_by_seq_id($hash_omniscient);

my $removeit=undef;
#################
# == LEVEL 1 == #
#################
foreach my $seqid (sort { (($a =~ /(\d+)$/)[0] || 0) <=> (($b =~ /(\d+)$/)[0] || 0) } keys %{$hash_sortBySeq}){ # loop over all the feature level1

	foreach my $tag_l1 (sort {$a cmp $b} keys %{$hash_omniscient->{'level1'}}){
		foreach my $feature_l1 ( @{$hash_sortBySeq->{$seqid}{$tag_l1}} ){
			my $id_l1 = lc($feature_l1->_tag_value('ID'));

	    $removeit = check_feature($feature_l1, 'level1', \@ptagList, \@attListOk, $opt_test);
			# we can remove feature L1 now because we are looping over $hash_sortBySeq not $hash_omniscient
			if ($removeit){
	      my $cases = remove_l1_and_relatives($hash_omniscient, $feature_l1, $fhout_discarded);
				$all_cases{'l1'} += $cases->{'l1'};
				$all_cases{'l2'} += $cases->{'l2'};
				$all_cases{'l3'} += $cases->{'l3'};
				$all_cases{'all'} += $cases->{'all'};
				next;
	    }

	    #################
	    # == LEVEL 2 == #
	    #################
			my @list_l2_to_remove;
			my @list_l3_to_remove;
	    foreach my $tag_l2 (sort keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...

	      if ( exists_keys( $hash_omniscient, ('level2', $tag_l2, $id_l1) ) ){
	        my @list_fl2 = @{$hash_omniscient->{'level2'}{$tag_l2}{$id_l1}};
	        foreach my $feature_l2 ( sort { $a->start <=> $b->start } @list_fl2 ) {

	          $removeit = check_feature($feature_l2,'level2', \@ptagList, \@attListOk, $opt_test);
	          if ($removeit){
							push @list_l2_to_remove, [$feature_l2, $tag_l1, $id_l1, $fhout_discarded];
	            next;
	          }
	          #################
	          # == LEVEL 3 == #
	          #################
	          my $id_l2 = lc($feature_l2->_tag_value('ID'));

	          foreach my $tag_l3 (sort keys %{$hash_omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...
	            if ( exists_keys( $hash_omniscient, ('level3', $tag_l3, $id_l2) ) ){
	              my @list_fl3 = @{$hash_omniscient->{'level3'}{$tag_l3}{$id_l2}};
	              foreach my $feature_l3 ( sort { $a->start <=> $b->start } @list_fl3 ) {

	                $removeit = check_feature($feature_l3, 'level3', \@ptagList, \@attListOk, $opt_test);
	                if ($removeit){
										push @list_l3_to_remove, [$feature_l3, $tag_l1, $id_l1, $tag_l2, $id_l2, $fhout_discarded];
	                }
	              }
	            }
	          }
	        }
	      }
	    }
			# Should be removed after looping over them to avoid problems
			if (@list_l2_to_remove) {
				foreach my $infos (@list_l2_to_remove) {
					my $cases = remove_l2_and_relatives( $hash_omniscient, @$infos);
					$all_cases{'l1'} += $cases->{'l1'};
					$all_cases{'l2'} += $cases->{'l2'};
					$all_cases{'l3'} += $cases->{'l3'};
					$all_cases{'all'} += $cases->{'all'};
				}
			}
			if (@list_l3_to_remove) {
				foreach my $infos (@list_l3_to_remove) {
					my $cases = remove_l3_and_relatives( $hash_omniscient, @$infos);
					$all_cases{'l1'} += $cases->{'l1'};
					$all_cases{'l2'} += $cases->{'l2'};
					$all_cases{'l3'} += $cases->{'l3'};
					$all_cases{'all'} += $cases->{'all'};
				}
			}
	  }
	}
}

print_omniscient( {omniscient => $hash_omniscient, output => $gffout_ok} );

# REPORT
my $stringPrint2 = $all_cases{'all'}." features removed:\n";
$stringPrint2 .= $all_cases{'l1'}." features level1 (e.g. gene) removed\n";
$stringPrint2 .= $all_cases{'l2'}." features level2 (e.g. mRNA) removed\n";
$stringPrint2 .= $all_cases{'l3'}." features level3 (e.g. exon) removed\n";

print $ostreamReport $stringPrint2 if ($opt_output);
dual_print1 $stringPrint2;

# --- final messages ---
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

sub check_feature{
  my  ($feature, $level, $ptagList, $attListOk, $opt_test)=@_;

  my $removeit=undef;
  my $primary_tag=$feature->primary_tag;

  # check primary tag (feature type) to handle
  foreach my $ptag (@$ptagList){

    if($ptag eq "all"){
      $removeit = should_we_remove_feature($feature, $attListOk, $opt_test);
    }
    elsif(lc($ptag) eq $level){
      $removeit = should_we_remove_feature($feature, $attListOk, $opt_test);
    }
    elsif(lc($ptag) eq lc($primary_tag) ){
      $removeit = should_we_remove_feature($feature, $attListOk, $opt_test);
    }
  }
  return $removeit;
}

sub should_we_remove_feature{
  my ($feature, $attListOk, $opt_test)=@_;

  foreach my $att ( @{$attListOk} ){
    if ($feature->has_tag($att)){
      $opt_test ? return 0 : return 1;

    }
  }
  $opt_test ? return 1 : return 0;
}

__END__

=head1 NAME

agat_sp_filter_feature_by_attribute_presence.pl

=head1 DESCRIPTION

The script aims to filter features according to attribute presence (9th column).
If the attribute exists, the feature is discarded.
Attribute are stored in the 9th column and have this shape: tag=value
/!\ Removing a level1 or level2 feature will automatically remove all linked subfeatures, and
removing all children of a feature will automatically remove this feature too.

=head1 SYNOPSIS

    agat_sp_filter_feature_by_attribute_presence.pl --gff infile.gff -a <tag> [ --output outfile ]
    agat_sp_filter_feature_by_attribute_presence.pl --help

=head1 OPTIONS

=over 8

=item B<-f>, B<--reffile>, B<--gff>  or B<-ref>

Input GFF3 file that will be read

=item B<-p>,  B<--type> or  B<-l>

primary tag option, case insensitive, list. Allow to specied the feature types that will be handled.
You can specified a specific feature by given its primary tag name (column 3) as: cds, Gene, MrNa
You can specify directly all the feature of a particular level:
      level2=mRNA,ncRNA,tRNA,etc
      level3=CDS,exon,UTR,etc
By default all feature are taking into account. fill the option by the value "all" will have the same behaviour.


=item B<--attribute>, B<--att>, B<-a>

String - Attributes tag specified will be used to filter the feature type (feature type can also be specified by the option -p). 
List of attribute tags must be coma separated.

=item B<--flip>

BOLEAN - In order to flip the test and keep features that do have the attribute and filter those without

=item B<-o> or B<--output>

Output GFF file.  If no output file is specified, the output will be
written to STDOUT.

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
