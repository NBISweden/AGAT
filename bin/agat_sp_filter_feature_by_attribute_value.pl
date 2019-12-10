#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use POSIX qw(strftime);
use List::MoreUtils qw(uniq);
use Scalar::Util qw(looks_like_number);
use Pod::Usage;
use Bio::Tools::GFF;
use IO::File;
use AGAT::Omniscient;

my $header = get_agat_header();
my $primaryTag=undef;
my $opt_output= undef;
my $opt_value = undef;
my $opt_attribute = undef;
my $opt_test = "=";
my $opt_gff = undef;
my $opt_verbose = undef;
my $opt_help;

# OPTION MANAGMENT
my @copyARGV=@ARGV;
if ( !GetOptions( 'f|ref|reffile|gff=s' => \$opt_gff,
                  'value=s'             => \$opt_value,
                  "p|type|l=s"          => \$primaryTag,
                  'a|attribute=s'       => \$opt_attribute,
                  't|test=s'            => \$opt_test,
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

if ( ! $opt_gff or ! defined($opt_value) or ! $opt_attribute ){
    pod2usage( {
           -message => "$header\nAt least 3 parameters are mandatory:\n1) Input reference gff file: -f\n".
           "2) An attribute tag: -a\n3) A value (string or int) that will be used for filtering: --value\n\n",
           -verbose => 0,
           -exitval => 2 } );
}

###############
# Test options
if($opt_test ne "<" and $opt_test ne ">" and $opt_test ne "<=" and $opt_test ne ">=" and $opt_test ne "="){
  print "The test to apply is Wrong: $opt_test.\nWe want something among this list: <,>,<=,>=,! or =.";exit;
}
if( ! looks_like_number($opt_value) ){
  if($opt_test eq "="){$opt_test="eq";}
  elsif($opt_test eq "!"){$opt_test="ne";}
  else{ print "This test $opt_test is not possible with string value.";exit; }
}

###############
# Manage Output

## FOR GFF FILE
my $gffout_ok; my $fhout_discarded ; my $ostreamReport;

if ($opt_output) {
  my ($outfile,$path,$ext) = fileparse($opt_output,qr/\.[^.]*/);

  # set file names
  my $outfile_ok = $path.$outfile.$ext;
  my $outfile_discarded = $path.$outfile."_discarded.txt";
  my $outfile_report = $path.$outfile."_report.txt";

  # check existence
  if(-f $outfile_ok){  print "File $outfile_ok already exist.\n";exit;}
  if(-f $outfile_discarded){  print "File $outfile_discarded already exist.\n";exit;}
  if(-f $outfile_report){  print "File $outfile_report already exist.\n";exit;}

  # create fh
  open( my $fh, '>', $outfile_ok) or die "Could not open file $outfile_ok $!";
  $gffout_ok = Bio::Tools::GFF->new(-fh => $fh, -gff_version => 3 );

  open( $fhout_discarded, '>', $outfile_discarded) or die "Could not open file $outfile_discarded $!";


  open($ostreamReport, '>', $outfile_report) or die "Could not open file $outfile_report $!";

}
else{
  $gffout_ok = Bio::Tools::GFF->new(-fh => \*STDOUT, -gff_version => 3);
  $fhout_discarded = Bio::Tools::GFF->new(-fh => \*STDOUT, -gff_version => 3);
  $ostreamReport = Bio::Tools::GFF->new(-fh => \*STDOUT, -gff_version => 3);
}

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

# start with some interesting information
my $stringPrint = strftime "%m/%d/%Y at %Hh%Mm%Ss", localtime;
$stringPrint .= "\nusage: $0 @copyARGV\n";
$stringPrint .= "We will discard $print_feature_string that have the attribute $opt_attribute with the value $opt_test $opt_value.\n";

if ($opt_output){
  print $ostreamReport $stringPrint;
  print $stringPrint;
}
else{ print $stringPrint; }
                                                      #######################
                                                      #        MAIN         #
#                     >>>>>>>>>>>>>>>>>>>>>>>>>       #######################       <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
my $cases=0;
######################
### Parse GFF input #
my ($hash_omniscient, $hash_mRNAGeneLink) =  slurp_gff3_file_JD({ input => $opt_gff,
                                                                  verbose => $opt_verbose
                                                                });
print("Parsing Finished\n");
### END Parse GFF input #
#########################

my $removeit=undef;
foreach my $tag_l1 ( sort keys %{$hash_omniscient->{'level1'}}){
  foreach my $id_l1 (sort keys %{$hash_omniscient->{'level1'}{$tag_l1} } ){

    my $feature_l1=$hash_omniscient->{'level1'}{$tag_l1}{$id_l1};

    $removeit = check_feature($feature_l1, 'level1', \@ptagList, $opt_attribute, $opt_test, $opt_value);
    if ($removeit){
      remove_l1_and_subfeature($hash_omniscient, $feature_l1, $tag_l1, $id_l1);
      next;
    }

    #################
    # == LEVEL 2 == #
    #################
    foreach my $tag_l2 (sort keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...

      if ( exists_keys( $hash_omniscient, ('level2', $tag_l2, $id_l1) ) ){
        my @list_fl2 = @{$hash_omniscient->{'level2'}{$tag_l2}{$id_l1}};
        foreach my $feature_l2 ( @list_fl2 ) {

          $removeit = check_feature($feature_l2,'level2', \@ptagList, $opt_attribute, $opt_test, $opt_value);
          if ($removeit){
            remove_l2_and_subfeature($hash_omniscient, $feature_l2, $tag_l1, $tag_l2, $id_l1);
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

                $removeit = check_feature($feature_l3, 'level3', \@ptagList, $opt_attribute, $opt_test, $opt_value);
                if ($removeit){
                  remove_l3($hash_omniscient, $feature_l3, $tag_l1, $tag_l2, $tag_l3, $id_l1, $id_l2);
                }
              }
            }
          }
        }
      }
    }
  }
}

print_omniscient($hash_omniscient, $gffout_ok); #print gene modified


$stringPrint = "$cases features removed \n";
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
  my  ($feature, $level, $ptagList, $opt_attribute, $opt_test, $opt_value)=@_;

  my $removeit=undef;
  my $primary_tag=$feature->primary_tag;

  # check primary tag (feature type) to handle
  foreach my $ptag (@$ptagList){

    if($ptag eq "all"){
      $removeit = should_we_remove_feature($feature, $opt_attribute, $opt_test, $opt_value);
    }
    elsif(lc($ptag) eq $level){
      $removeit = should_we_remove_feature($feature, $opt_attribute, $opt_test, $opt_value);
    }
    elsif(lc($ptag) eq lc($primary_tag) ){
      $removeit = should_we_remove_feature($feature, $opt_attribute, $opt_test, $opt_value);
    }
  }
  return $removeit;
}

sub should_we_remove_feature{
  my ($feature, $opt_attribute, $opt_test, $opt_value)=@_;

  if ($feature->has_tag($opt_attribute)){

    # get list of values for the attribute
    my @values = $feature->get_tag_values($opt_attribute);

    # if we found among the values one pass the test we return 1
    foreach my $value (@values){

      if ($opt_test eq "eq"){
        if ($value eq $opt_value){return 1; }
      }
      elsif ($opt_test eq "ne"){
        if ($value ne $opt_value){return 1; }
      }
      elsif ($opt_test eq "="){
       if ($value == $opt_value){return 1; }
      }
      elsif ($opt_test eq "!"){
        if ($value != $opt_value){return 1; }
      }
      elsif ($opt_test eq ">"){
        if ($value > $opt_value){return 1; }
      }
      elsif ($opt_test eq "<"){
        if ($value < $opt_value){return 1; }
      }
      elsif ($opt_test eq "<="){
        if ($value <= $opt_value){return 1; }
      }
      elsif ($opt_test eq ">="){
        if ($value >= $opt_value){return 1; }
      }
    }
  }
  return 0;
}

# remove from omniscient l1 feature and all subfeatures
sub remove_l1_and_subfeature{
  my ($omniscient, $feature, $tag_l1, $id_l1)=@_;

  foreach my $ptag_l2 (keys %{$omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...

    if ( exists_keys( $omniscient, ('level2', $ptag_l2, $id_l1) ) ){
      foreach my $feature_l2 ( @{$omniscient->{'level2'}{$ptag_l2}{$id_l1}}) {

        my $level2_ID = lc($feature_l2->_tag_value('ID'));

        foreach my $ptag_l3 (keys %{$omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...
          if ( exists_keys( $omniscient, ('level3', $ptag_l3, $level2_ID) ) ){
            foreach my $feature_l3 ( @{$omniscient->{'level3'}{$ptag_l3}{$level2_ID}}) {
              $cases++; print $fhout_discarded $feature_l3->_tag_value('ID')."\n";
            }
            delete $omniscient->{'level3'}{$ptag_l3}{$level2_ID} # delete level3
          }
        }
        $cases++;  print $fhout_discarded $feature_l2->_tag_value('ID')."\n";
      }
      delete $omniscient->{'level2'}{$ptag_l2}{$id_l1} # delete level2
    }
  }
  $cases++; print $fhout_discarded $feature->gff_string()."\n";
  delete $omniscient->{'level1'}{$tag_l1}{$id_l1}; # delete level1
}

# remove from omniscient l2 feature and all subfeatures
sub remove_l2_and_subfeature{
  my ($omniscient, $feature, $ptag_l1, $ptag_l2, $id_l1)=@_;

  my $level2_Parent_ID = lc($feature->_tag_value('Parent'));
  my $level2_ID = lc($feature->_tag_value('ID'));

  if ( exists_keys( $omniscient, ('level2', $ptag_l2, $id_l1) ) ){ # just extra security in case
    foreach my $feature_l2 ( @{$omniscient->{'level2'}{$ptag_l2}{$id_l1}}) {

      if($level2_ID eq lc($feature_l2->_tag_value('ID')) ){
        # let's delete all l3 subfeatures before to remove the l2
        foreach my $ptag_l3 (keys %{$omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...
          if ( exists_keys( $omniscient, ('level3', $ptag_l3, $level2_ID)  ) ){
            foreach my $feature_l3 ( @{$omniscient->{'level3'}{$ptag_l3}{$level2_ID}}) {
              $cases++; print $fhout_discarded $feature_l3->_tag_value('ID')."\n";
            }
            delete $omniscient->{'level3'}{$ptag_l3}{$level2_ID} # delete level3
          }
        }
      }
    }

    # delete level2 and the hash pointer if the list is empty (no uisoform left)
    my @id_concern_list=($id_l1);
    my @id_list_to_remove=($level2_ID);
    my @list_tag_key=('all');
    $cases++; print $fhout_discarded $feature->gff_string()."\n";
    remove_element_from_omniscient(\@id_concern_list, \@id_list_to_remove, $omniscient, 'level2','false', \@list_tag_key);


    if( ! exists_keys($omniscient, ('level2', $ptag_l2, $id_l1) ) ){
      #The list was empty so l2 has been removed, we can now remove l1
      if( exists_keys($hash_omniscient, ('level1', $ptag_l1, $id_l1) ) ){
        $cases++; print $fhout_discarded $id_l1."\n";
        delete $hash_omniscient->{'level1'}{$ptag_l1}{$id_l1};
      }
    }
  }
}

# remove from omniscient l2 feature and all subfeatures
sub remove_l3{
  my ($omniscient, $feature, $ptag_l1, $ptag_l2, $ptag_l3, $id_l1, $id_l2)=@_;

  my $level3_Parent_ID = lc($feature->_tag_value('Parent'));
  my $id_l3 = lc($feature->_tag_value('ID'));

  if ( exists_keys( $omniscient, ('level3', $ptag_l3, $id_l2) ) ){ # just extra security in case
    foreach my $feature_l3 ( @{$omniscient->{'level3'}{$ptag_l3}{$id_l2}}) {

      if($id_l3 eq lc($feature_l3->_tag_value('ID')) ){
        #remove one feature and pointer if no more feature left in the list
        my @id_concern_list=($level3_Parent_ID);
        my @id_list_to_remove=($id_l3);
        my @list_tag_key=('all');
        $cases++; print $fhout_discarded $feature->gff_string()."\n";
        remove_element_from_omniscient(\@id_concern_list, \@id_list_to_remove, $omniscient, 'level3','false', \@list_tag_key);
      }
    }
  }

  # List empty check if we remove l2 or other l3 linked to it
  if( ! exists_keys($omniscient, ('level3', $ptag_l3, $id_l2)) ){
    foreach my $tag_l3 (keys %{$hash_omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...
      if ( exists_keys( $hash_omniscient, ('level3', $tag_l3, $id_l2) ) ){
        return;
      }
    }
    print "Go there\n";
    # if we arrive here it means no more L3 feature is attached to L2
    # we remove the L2 parent properly (if isoforms they are kept)
    my @id_concern_list=($id_l1);
    my @id_list_to_remove=($id_l2);
    my @list_tag_key=('all');
    $cases++; print $fhout_discarded $id_l2."\n";
    remove_element_from_omniscient(\@id_concern_list, \@id_list_to_remove, $omniscient, 'level2','false', \@list_tag_key);

    # Should we remove l1 too?
    if( ! exists_keys($omniscient, ('level2', $ptag_l2, $id_l1) ) ){
      #The list was empty so l2 has been removed, we can now remove l1
      if( exists_keys($hash_omniscient, ('level1', $ptag_l1, $id_l1) ) ){
        $cases++; print $fhout_discarded $id_l1."\n";
        delete $hash_omniscient->{'level1'}{$ptag_l1}{$id_l1}
      }
    }
  }
}

__END__

=head1 NAME

agat_sp_select_feature_by_attribute_value.pl

=head1 DESCRIPTION

The script aims to filter features according to attribute value (9th column).
If the attribute tag is missing the feature will not be discarded.
Attribute are stored in the 9th column and have this shape: tag=value
/!\ Removing a level1 or level2 feature will automatically remove all linked subfeatures, and
removing all children of a feature will automatically remove this feature too.

Example case, removing all

=head1 SYNOPSIS

    ./agat_sp_select_feature_by_attribute_value.pl -f infile.gff --value 1 -t "=" [ --output outfile ]
    ./agat_sp_select_feature_by_attribute_value.pl --help

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

=item B<-v> or B<--value>

Value to check in the attribute

=item B<-t> or B<--test>
Test to apply (> < = >= <=). default value "=". If you us one of these two character >, <, please don't forget to quote you parameter liket that "<=". Else your terminal will complain.

=item B<-o> or B<--output>

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
