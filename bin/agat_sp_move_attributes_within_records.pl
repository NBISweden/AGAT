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
my $primaryTagCopy="level2";
my $primaryTagPaste="level3";
my $opt_output= undef;
my $attributes="all_attributes";
my $opt_gff = undef;
my $opt_verbose = undef;
my $opt_help;

# OPTION MANAGMENT
my @copyARGV=@ARGV;
if ( !GetOptions( 'f|ref|reffile|gff=s'  => \$opt_gff,
                  "feature_copy|fc=s"    => \$primaryTagCopy,
                  "feature_paste|fp=s"   => \$primaryTagPaste,
                  'o|output=s'           => \$opt_output,
                  "tag|att=s"            => \$attributes,
                  'v|verbose!'           => \$opt_verbose,
                  'c|config=s'           => \$config,
                    'thread|threads|cpu|cpus|core|cores|job|jobs=i' => \$cpu,
                  'h|help!'              => \$opt_help ) )
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

if ( ! $opt_gff  ){
    pod2usage( {
           -message => "$header\nAt least 1 parameter is mandatory: Input reference gff file: --gff\n\n",
           -verbose => 0,
           -exitval => 2 } );
}

# --- Manage config ---
initialize_agat({ config_file_in => $config, input => $opt_gff });
$CONFIG->{cpu} = $cpu if defined($cpu);

###############
# Manage Output

## FOR GFF FILE
my $gffout_ok_file ;

if ($opt_output) {
  my ($outfile,$path,$ext) = fileparse($opt_output,qr/\.[^.]*/);

  # set file names
  $gffout_ok_file = $path.$outfile.$ext;
}
my $gffout_ok = prepare_gffout( $gffout_ok_file);

# Manage $primaryTag for copy
my @ptagListCopy;
my $print_feature_string_copy;
if(! $primaryTagCopy or $primaryTagCopy eq "all"){
  $print_feature_string_copy = "all features";
  push(@ptagListCopy, "all");
}
elsif($primaryTagCopy =~/^level[123]$/){
  $print_feature_string_copy .= "$primaryTagCopy features ";
  push(@ptagListCopy, $primaryTagCopy);
}
else{
   @ptagListCopy= split(/,/, $primaryTagCopy);
   foreach my $tag (@ptagListCopy){
      if($tag =~/^level[123]$/){
        $print_feature_string_copy .= "$primaryTagCopy features ";
      }
      else{
        $print_feature_string_copy .= "$tag feature ";
      }
   }
}

# Manage $primaryTag to paste
my @ptagListPaste;
my $print_feature_string_paste;
if(! $primaryTagPaste or $primaryTagPaste eq "all"){
  $print_feature_string_paste = "all features";
  push(@ptagListPaste, "all");
}
elsif($primaryTagPaste =~/^level[123]$/){
  $print_feature_string_paste .= "$primaryTagPaste features ";
  push(@ptagListPaste, $primaryTagPaste);
}
else{
   @ptagListPaste= split(/,/, $primaryTagPaste);
   foreach my $tag (@ptagListPaste){
      if($tag =~/^level[123]$/){
        $print_feature_string_paste .= "$primaryTagPaste features ";
      }
      else{
        $print_feature_string_paste .= "$tag feature ";
      }
   }
}

# Manage attributes if given
### If attributes given, parse them:
my %attHashOk;
my @attListOk;
if ($attributes){

  if ($attributes eq "all_attributes"){
    print "All attributes will be used !\n";
    $attHashOk{"all_attributes"}++;
  }
  else{
    my @attList= split(/,/, $attributes);

    foreach my $attribute (@attList){

      if($attribute == 0){ # Attribute alone
        #check for ID attribute
        if(lc($attribute) eq "id" ){print "ID attribute cannot be modified !\n";exit;}
        #check for Parent attribute
        if(lc($attribute) eq "parent"){print "Parent attribute cannot be modified !\n";exit;}
        $attHashOk{$attribute}++;
        push(@attListOk, $attribute);
      }
    }
  }
  print "\n";
}

# start with some interesting information
my $stringPrint = strftime "%m/%d/%Y at %Hh%Mm%Ss", localtime;
$stringPrint .= "\nusage: $0 @copyARGV\n";
$stringPrint .= "The attributes @attListOk from the following feature types: $print_feature_string_copy will be copy pasted to the following feature types: $print_feature_string_paste.\n";
print $stringPrint;

                          #######################
# >>>>>>>>>>>>>>>>>>>>>>>>#        MAIN         #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                          #######################
my %all_cases = ('l1' => 0, 'l2' => 0, 'l3' => 0, 'all' => 0);
######################
### Parse GFF input #
my ($hash_omniscient, $hash_mRNAGeneLink) =  slurp_gff3_file_JD({ input => $opt_gff });

### END Parse GFF input #
#########################
# sort by seq id
my $hash_sortBySeq = gather_and_sort_l1_by_seq_id($hash_omniscient);

#################
# == LEVEL 1 == #
#################
foreach my $seqid (sort { (($a =~ /(\d+)$/)[0] || 0) <=> (($b =~ /(\d+)$/)[0] || 0) } keys %{$hash_sortBySeq}){ # loop over all the feature level1

	foreach my $tag_l1 (sort {$a cmp $b} keys %{$hash_omniscient->{'level1'}}){
		foreach my $feature_l1 ( @{$hash_sortBySeq->{$seqid}{$tag_l1}} ){
			my $id_l1 = lc($feature_l1->_tag_value('ID'));
		  my $copyl1 = check_feature($feature_l1, 'level1', \@ptagListCopy);

		  #################
		  # == LEVEL 2 == #
		  #################
		  foreach my $tag_l2 (sort keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...

		    if ( exists_keys( $hash_omniscient, ('level2', $tag_l2, $id_l1) ) ){
			    my @list_fl2 = @{$hash_omniscient->{'level2'}{$tag_l2}{$id_l1}};
	        foreach my $feature_l2 ( @list_fl2 ) {

            my $copyl2 = check_feature($feature_l2,'level2', \@ptagListCopy);
            my $pastel2 = check_feature($feature_l2,'level2', \@ptagListPaste);
            if($pastel2 and $copyl1){
              add_tags_from_to($feature_l1, $feature_l2);
            }

				    #################
				    # == LEVEL 3 == #
				    #################
			      my $id_l2 = lc($feature_l2->_tag_value('ID'));

				    foreach my $tag_l3 (sort keys %{$hash_omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...
				      if ( exists_keys( $hash_omniscient, ('level3', $tag_l3, $id_l2) ) ){
				       	my @list_fl3 = @{$hash_omniscient->{'level3'}{$tag_l3}{$id_l2}};
				       	foreach my $feature_l3 ( @list_fl3 ) {

					        my $copyl3 = check_feature($feature_l3, 'level3', \@ptagListCopy);
                  my $pastel3 = check_feature($feature_l3, 'level3', \@ptagListPaste);
                  #copy attributes from L1 to l3
                  if ($copyl1 and $pastel3){
                    add_tags_from_to($feature_l1, $feature_l3);
                  }
                   #copy attributes from L2 to l3
                  if ($copyl2 and $pastel3){
                    add_tags_from_to($feature_l2, $feature_l3);
                  }
                  # ------- Case L3 to L3 -------
                  foreach my $tag_l3_again (sort keys %{$hash_omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...
				            if ( exists_keys( $hash_omniscient, ('level3', $tag_l3_again, $id_l2) ) ){
				       	      my @list_fl3_again = @{$hash_omniscient->{'level3'}{$tag_l3_again}{$id_l2}};
				              foreach my $feature_l3_again ( @list_fl3_again ) {

                        # Check it is not the same feature
                        if( lc($feature_l3->_tag_value('ID')) ne lc($feature_l3_again->_tag_value('ID')) ){

                          my $pastel3_again = check_feature($feature_l3_again, 'level3', \@ptagListPaste);
                          #copy attributes from L1 to l3
                          if ($copyl3 and $pastel3_again){
                              add_tags_from_to($feature_l3, $feature_l3_again);
                          }
                        }
                      }
                    }
                  }
					      }
					    }
					  }
            # ------- Case L2 to L2 -------
            foreach my $tag_l2_again (sort keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...

              if ( exists_keys( $hash_omniscient, ('level2', $tag_l2_again, $id_l1) ) ){
                my @list_fl2_again = @{$hash_omniscient->{'level2'}{$tag_l2_again}{$id_l1}};

                foreach my $feature_l2_again ( @list_fl2_again ) {
                  # Check it is not the same feature
                  if( lc($feature_l2->_tag_value('ID')) ne lc($feature_l2_again->_tag_value('ID')) ){

                    my $pastel2_again = check_feature($feature_l2_again, 'level3', \@ptagListPaste);
                    #copy attributes from L1 to l3
                    if ($copyl2 and $pastel2_again){
                      add_tags_from_to($feature_l2, $feature_l2_again);
                    }
                  }
                }
              }
            }
				  }
		    }
			}
		}
	}
}

# create omniscient with only selected recoreds
print_omniscient( {omniscient => $hash_omniscient, output => $gffout_ok} );#print gene modified in file

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

# add tags and avoid Parent and ID
sub add_tags_from_to{
  my  ($feature_from, $feature_to)=@_;

  # for each tag of feature_from
  my @list_tags = $feature_from->get_all_tags();
  foreach my $tag (@list_tags){
    # Check tag is among those that have to be used
    if(exists_keys(\%attHashOk,("all_attributes")) or exists_keys(\%attHashOk,($tag))){
      if (lc($tag) ne "id" and lc($tag) ne "parent" ){
        my @tag_values = $feature_from->get_tag_values($tag);
        create_or_append_tag($feature_to, $tag, \@tag_values);
      }
    }
  }
}

sub check_feature{
  my  ($feature, $level, $ptagList)=@_;

  my $keepit=undef;
  my $primary_tag=$feature->primary_tag;
  # check primary tag (feature type) to handle
	foreach my $ptag (@$ptagList){

	  if($ptag eq "all"){
	    $keepit = 1 ;
	  }
	  elsif(lc($ptag) eq $level){
	    $keepit = 1 ;
	  }
	  elsif(lc($ptag) eq lc($primary_tag) ){
	    $keepit = 1 ;
	  }
	}
  return $keepit;
}

__END__

=head1 NAME

agat_sp_move_attributes_within_records.pl

=head1 DESCRIPTION

The script aims to keep move attributes within a record e.g. from Level1 to Level2 and/or Level3 features; and / or from Level2 to Level2 or Level3 features; and / or from Level3 to Level3 features.
Example of L1 feature: gene
Example of L2 featrue

=head1 SYNOPSIS

    agat_sp_move_attributes_within_records.pl --gff infile.gff --feature_copy mRNA  --feature_paste CDS --attribute Dbxref,Ontology [ --output outfile ]
    agat_sp_move_attributes_within_records.pl --help

=head1 OPTIONS

=over 8

=item B<-f>, B<--reffile>, B<--gff>  or B<-ref>

Input GFF3 file that will be read

=item B<--feature_copy> or B<--fc>

primary tag (feature type) option to list from which feature we will copy the attributes, case insensitive. 
You can specified a feature (or a coma separated list) by giving its primary tag / feature type (column 3) value as: cds, Gene, MrNa, etc
You can specify directly all the feature of a particular level:
      level2=mRNA,ncRNA,tRNA,etc
      level3=CDS,exon,UTR,etc
By default all level2 feature are used. 

=item B<--feature_paste> or B<--fp>

primary tag (feature type) option to list to which feature we will paste the attributes, case sensitive. 
You can specified a feature (or a coma separated list) by giving its primary tag / feature type (column 3) value as: cds, Gene, MrNa, etc
You can specify directly all the feature of a particular level:
      level2=mRNA,ncRNA,tRNA,etc
      level3=CDS,exon,UTR,etc
By default all feature level3 are used. 


=item  B<-a> or B<--attribute>

Attribute that will be copied and pasted. Case sensitive.
You can specified an attribute (or a coma separated list) by giving its attribute tag value (column9) as: Ontology, Dbxref, etc
Default: all_attributes
/!\ <all_attributes> is a specific parameter meaning all the attributes will be use.


=item B<-o> or B<--output>

Output GFF file.  If no output file is specified, the output will be
written to STDOUT.

=item B<-v>

Verbose option for debugging purpose.

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
