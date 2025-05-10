#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use POSIX qw(strftime);
use Scalar::Util qw(looks_like_number);
use Pod::Usage;
use IO::File;
use AGAT::AGAT;

my $header = get_agat_header();
my $config;
my $cpu;
my $primaryTag=undef;
my $opt_output= undef;
my $opt_value = undef;
my $opt_keep_parental = undef;
my $opt_na_aside = undef;
my $opt_value_insensitive = undef;
my $opt_attribute = undef;
my $opt_test = "=";
my $opt_gff = undef;
my $opt_verbose = undef;
my $opt_help;

# OPTION MANAGMENT
my @copyARGV=@ARGV;
if ( !GetOptions( 'f|ref|reffile|gff=s' => \$opt_gff,
                  'value=s'             => \$opt_value,
                  'value_insensitive!'  => \$opt_value_insensitive,
                  'keep_parental!'      => \$opt_keep_parental,
                  'na_aside!'           => \$opt_na_aside, 
                  "p|type|l=s"          => \$primaryTag,
                  'a|attribute=s'       => \$opt_attribute,
                  't|test=s'            => \$opt_test,
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

if ( ! $opt_gff or ! defined($opt_value) or ! $opt_attribute ){
    pod2usage( {
           -message => "$header\nAt least 3 parameters are mandatory:\n1) Input reference gff file: --gff\n".
           "2) An attribute tag: -a\n3) A value (string or int) that will be used for filtering: --value\n\n",
           -verbose => 0,
           -exitval => 2 } );
}

# --- Manage config ---
initialize_agat({ config_file_in => $config, input => $opt_gff });
$CONFIG->{cpu} = $cpu if defined($cpu);

###############
# Test options
if($opt_test ne "<" and $opt_test ne ">" and $opt_test ne "<=" and $opt_test ne ">=" and $opt_test ne "=" and $opt_test ne "!"){
  print "The test to apply is Wrong: $opt_test.\nWe want something among this list: <,>,<=,>=,! or =.";exit;
}

###############
# Manage Output

## FOR GFF FILE
my $gffout_ok_file ;
my $fhout_discarded_file ;
my $ostreamReport_file;
my $fhout_semidDiscarded_file if $opt_na_aside;

if ($opt_output) {
  my ($outfile,$path,$ext) = fileparse($opt_output,qr/\.[^.]*/);

  # set file names
  $gffout_ok_file = $path.$outfile.$ext;
  $fhout_discarded_file = $path.$outfile."_discarded.gff";
  $ostreamReport_file = $path.$outfile."_report.txt";
  $fhout_semidDiscarded_file = $path.$outfile."_na.gff";
}

my $gffout_ok = prepare_gffout( $gffout_ok_file );
my $fhout_discarded = prepare_gffout( $fhout_discarded_file);
my $ostreamReport = prepare_fileout($ostreamReport_file);
my $fhout_semidDiscarded = prepare_gffout( $fhout_semidDiscarded_file) if $opt_na_aside;

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

# Transform value list into hash
my $value_hash = string_sep_to_hash({ string => $opt_value,
                                      separator => ","
                                    });
                                      
foreach my $value (keys %{$value_hash}){
  if( ! looks_like_number($value) ){
    if($opt_test ne "=" and $opt_test ne "!"){
      print "This test $opt_test is not possible with string value.\n";
      exit; 
    }
  }
}

# start with some interesting information
my $stringPrint = strftime "%m/%d/%Y at %Hh%Mm%Ss", localtime;
$stringPrint .= "\nusage: $0 @copyARGV\n";
$stringPrint .= "We will discard $print_feature_string that have the attribute $opt_attribute with the value $opt_test $opt_value";
if ($opt_value_insensitive){
  $stringPrint .= " case insensitive.\n";
}else{
   $stringPrint .= " case sensitive.\n";
}

if ($opt_output){
  print $ostreamReport $stringPrint;
  print $stringPrint;
}
else{ print $stringPrint; }
                          #######################
# >>>>>>>>>>>>>>>>>>>>>>>>#        MAIN         #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                          #######################
my %all_cases = ( 'left' => {'l1' => 0, 'l2' => 0, 'l3' => 0, 'all' => 0},
                  'discarded' => {'l1' => 0, 'l2' => 0, 'l3' => 0, 'all' => 0} );

######################
### Parse GFF input #
my ($hash_omniscient, $hash_mRNAGeneLink) =  slurp_gff3_file_JD({ input => $opt_gff });
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

	    $removeit = check_feature($feature_l1, 'level1');
			# we can remove feature L1 now because we are looping over $hash_sortBySeq not $hash_omniscient
	    if ($removeit){
        my $cases;
        if($removeit == 1){ 
          $cases = remove_l1_and_relatives($hash_omniscient, $feature_l1, $fhout_discarded);
          $all_cases{'discarded'}{'l1'} += $cases->{'l1'};
          $all_cases{'discarded'}{'l2'} += $cases->{'l2'};
          $all_cases{'discarded'}{'l3'} += $cases->{'l3'};
          $all_cases{'discarded'}{'all'} += $cases->{'all'};
        } elsif ($removeit == 2 and $opt_na_aside){ 
          $cases = remove_l1_and_relatives($hash_omniscient, $feature_l1, $fhout_semidDiscarded);
          $all_cases{'na'}{'l1'} += $cases->{'l1'};
          $all_cases{'na'}{'l2'} += $cases->{'l2'};
          $all_cases{'na'}{'l3'} += $cases->{'l3'};
          $all_cases{'na'}{'all'} += $cases->{'all'};
        }
				next;
	    }

	    #################
	    # == LEVEL 2 == #
	    #################
			my %list_l2_to_remove;
			my %list_l3_to_remove;
	    foreach my $tag_l2 (sort keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...

	      if ( exists_keys( $hash_omniscient, ('level2', $tag_l2, $id_l1) ) ){
	        my @list_fl2 = @{$hash_omniscient->{'level2'}{$tag_l2}{$id_l1}};
	        foreach my $feature_l2 ( @list_fl2 ) {

	          $removeit = check_feature($feature_l2,'level2');
	          if ($removeit){
              if ($removeit == 1){
                push ( @{$list_l2_to_remove{'discarded'}}, [$feature_l2, $tag_l1, $id_l1, $fhout_discarded]);
              } elsif ($removeit == 2 and $opt_na_aside){
                push ( @{$list_l2_to_remove{'na'}}, [$feature_l2, $tag_l1, $id_l1, $fhout_semidDiscarded]);
              }
	            next;
	          }
	          #################
	          # == LEVEL 3 == #
	          #################
	          my $id_l2 = lc($feature_l2->_tag_value('ID'));

	          foreach my $tag_l3 (sort keys %{$hash_omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...
	            if ( exists_keys( $hash_omniscient, ('level3', $tag_l3, $id_l2) ) ){
	              my @list_fl3 = @{$hash_omniscient->{'level3'}{$tag_l3}{$id_l2}};
	              foreach my $feature_l3 ( @list_fl3 ) {

	                $removeit = check_feature($feature_l3, 'level3');
                  if($removeit ==  1){
                    push ( @{$list_l3_to_remove{'discarded'}}, [$feature_l3, $tag_l1, $id_l1, $tag_l2, $id_l2, $fhout_discarded]); 
                  } elsif ( $removeit == 2 and $opt_na_aside ){
                    push ( @{$list_l3_to_remove{'na'}}, [$feature_l3, $tag_l1, $id_l1, $tag_l2, $id_l2, $fhout_semidDiscarded]);
	                }
	              }
	            }
	          }
	        }
	      }
	    }
			# Should be removed after looping over them to avoid problems
			foreach my $key ( keys %list_l2_to_remove ){
				foreach my $infos ( @{$list_l2_to_remove{$key}} ) {
					my $cases = remove_l2_and_relatives( $hash_omniscient, @$infos, $opt_keep_parental);
					$all_cases{$key}{'l1'} += $cases->{'l1'};
					$all_cases{$key}{'l2'} += $cases->{'l2'};
					$all_cases{$key}{'l3'} += $cases->{'l3'};
					$all_cases{$key}{'all'} += $cases->{'all'};
				}
			}
			foreach my $key ( sort keys %list_l3_to_remove ){
				foreach my $infos ( @{$list_l3_to_remove{$key}} ) {
					my $cases = remove_l3_and_relatives( $hash_omniscient, @$infos, $opt_keep_parental);
					$all_cases{$key}{'l1'} += $cases->{'l1'};
					$all_cases{$key}{'l2'} += $cases->{'l2'};
					$all_cases{$key}{'l3'} += $cases->{'l3'};
					$all_cases{$key}{'all'} += $cases->{'all'};
				}
			}
		}
  }
}

print_omniscient( {omniscient => $hash_omniscient, output => $gffout_ok} );

$stringPrint = "Feature discarded by applying the test (see $fhout_discarded_file file):\n";
$stringPrint .= $all_cases{'discarded'}{'all'}." features removed:\n";
$stringPrint .= $all_cases{'discarded'}{'l1'}." features level1 (e.g. gene) removed\n";
$stringPrint .= $all_cases{'discarded'}{'l2'}." features level2 (e.g. mRNA) removed\n";
$stringPrint .= $all_cases{'discarded'}{'l3'}." features level3 (e.g. exon) removed\n";

if($opt_na_aside){
  $stringPrint .= "Feature left out because the attribute is missing (see $fhout_semidDiscarded_file file):\n";
  $stringPrint .= $all_cases{'na'}{'all'}." features removed:\n";
  $stringPrint .= $all_cases{'na'}{'l1'}." features level1 (e.g. gene) removed\n";
  $stringPrint .= $all_cases{'na'}{'l2'}." features level2 (e.g. mRNA) removed\n";
  $stringPrint .= $all_cases{'na'}{'l3'}." features level3 (e.g. exon) removed\n";
}

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

sub check_feature{
  my  ($feature, $level)=@_;

  my $removeit=undef;
  my $primary_tag=$feature->primary_tag;

  # check primary tag (feature type) to handle
  foreach my $ptag (@ptagList){

    if($ptag eq "all"){
      $removeit = should_we_remove_feature($feature);
    }
    elsif(lc($ptag) eq $level){
      $removeit = should_we_remove_feature($feature);
    }
    elsif(lc($ptag) eq lc($primary_tag) ){
      $removeit = should_we_remove_feature($feature);
    }
  }
  return $removeit;
}

sub should_we_remove_feature{
  my ($feature)=@_;

  if ($feature->has_tag($opt_attribute)){

    # get list of values for the attribute
    my @file_values = $feature->get_tag_values($opt_attribute);

    # if we found among the values one pass the test we return 1
    foreach my $file_value (@file_values){

      foreach my $given_value (keys %{$value_hash}){
        # Deal with insensitive for template
        if ($opt_value_insensitive){
          $given_value = lc($given_value);
          $file_value = lc($file_value);
        }
        # for string values replace = by eq and ! by ne and avoid other type of test
        if ( ! looks_like_number ($given_value) or ! looks_like_number ($file_value)){
          print "String case\n" if $opt_verbose;
          if ($opt_test eq "="){
            if ($file_value eq $given_value) { print "equal\n" if $opt_verbose; return 1; }
            else { print "not equal\n" if $opt_verbose; }
          }
          elsif ($opt_test eq "!"){
            if ($file_value ne $given_value){ print "different\n" if $opt_verbose; return 1; }
            else { print "not different\n" if $opt_verbose; }
          }
        } 
        else{
          print "Number case\n" if $opt_verbose;
          if ($opt_test eq "="){
            if ($file_value == $given_value){return 1; }
          }
          elsif ($opt_test eq "!"){
            if ($file_value != $given_value){return 1; }
          }
          elsif ($opt_test eq ">"){
            if ($file_value > $given_value){return 1; }
          }
          elsif ($opt_test eq "<"){
            if ($file_value < $given_value){return 1; }
          }
          elsif ($opt_test eq "<="){
            if ($file_value <= $given_value){return 1; }
          }
          elsif ($opt_test eq ">="){
            if ($file_value >= $given_value){return 1; }
          }
        }
      }
    }
    return 0;
  } else {
    print "Attribute not found  case\n" if $opt_verbose;
    return 2;
  }
}

__END__

=head1 NAME

agat_sp_filter_feature_by_attribute_value.pl

=head1 DESCRIPTION

The script aims to filter features according to attribute value (9th column).
- If the attribute exists and the value do not pass the test, the feature is written into <output>.
- If the attribute exists and the value pass the test, the feature is discarded and written into <output>_discarded.gff.
- If the attribute tag is missing (test cannot be applyed), the feature will be written into <output> by default. If --na_aside parameter 
is activated then it will be written into <output>_na.gff.

Attribute are stored in the 9th column and have this shape: tag=value
/!\ Removing a level1 or level2 feature will automatically remove all linked subfeatures.
/!\ Removing all children of a feature will automatically remove this feature too (excepted if --keep_parental is activated).
/!\ If --keep_parental is not activated and --na_aside is activated, and all level3 features of a record are split between both <output>_na.gff and <output>_discarded.gff, 
then the parental level1 and level2 features are removed and will end up in the <output>_na.gff file only.

=head1 SYNOPSIS

    agat_sp_filter_feature_by_attribute_value.pl --gff infile.gff --value 1 -t "=" [ --output outfile ]
    agat_sp_filter_feature_by_attribute_value.pl --help

=head1 OPTIONS

=over 8

=item B<-f>, B<--reffile>, B<--gff>  or B<-ref>

Input GFF3 file that will be read

=item  B<-a> or B<--attribute>

Attribute tag to specify the attribute to analyse (attribute example: tag=value).

=item B<-p>,  B<--type> or  B<-l>

primary tag option, case insensitive, list. Allow to specied the feature types that will be handled.
You can specified a specific feature by given its primary tag name (column 3) as: cds, Gene, MrNa
You can specify directly all the feature of a particular level:
      level2=mRNA,ncRNA,tRNA,etc
      level3=CDS,exon,UTR,etc
By default all feature are taking into account. fill the option by the value "all" will have the same behaviour.

=item B<--value>

Value(s) to check in the attribute. Case sensitive. List of values must be coma separated.

=item B<--value_insensitive>

Bolean. Deactivated by default. When activated the values provided by the --value parameter are handled case insensitive.

=item B<--na_aside>

Bolean. Deactivated by default. By default if the attribute tag on which the filter is based is missing, the feature will be written into <output>.
When activated, such features will be written into a separate file called <output>_na.gff.

=item B<--keep_parental>

Bolean. Deactivated by default. When activated even if all child features have been removed, the parental one will be kept.

=item B<-t> or B<--test>

Test to apply (> < = ! >= <=). default value "=". 
If you use one of these two character >, <, please don't forget to quote the
parameter like that "<=" otherwise your terminal will complain.
Only = and ! tests can be used to compare string values.

=item B<-o> or B<--output>

Output GFF file. If no output file is specified, the output will be
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
