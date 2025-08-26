#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use POSIX qw(strftime);
use Pod::Usage;
use IO::File;
use AGAT::AGAT;

my $header = get_agat_header();
my ( $opt, $usage, $config ) = AGAT::AGAT::describe_script_options(
    $header,
    [ 'gff|f|ref|reffile=s',  'Input reference gff file', { required => 1 } ],
    [ 'feature_copy|fc=s',    'Feature types to copy from', { default => 'level2' } ],
    [ 'feature_paste|fp=s',   'Feature types to paste to', { default => 'level3' } ],
    [ 'attribute|tag|att|a=s', 'Comma-separated attributes to move', { default => 'all_attributes',
        callbacks => { valid => sub {
            my @att = split(/,/, shift);
            die 'ID attribute cannot be modified !'     if grep { lc($_) eq 'id' } @att;
            die 'Parent attribute cannot be modified !' if grep { lc($_) eq 'parent' } @att;
            return 1;
        } } } ],
);

my $opt_gff         = $opt->gff;
my $primaryTagCopy  = $opt->feature_copy;
my $primaryTagPaste = $opt->feature_paste;
my $attributes      = $opt->attribute;
my $opt_output      = $opt->out;
my $opt_verbose     = $config->{verbose};

my @copyARGV = @ARGV;

my $log;
if ( my $log_name = $config->{log_path} ) {
    open( $log, '>', $log_name ) or die "Can not open $log_name for printing: $!";
}
dual_print( $log, $header);

###############
# Manage Output

## FOR GFF FILE
my $gffout_ok_file ;

if ($opt_output) {
  my ($outfile,$path,$ext) = fileparse($opt_output,qr/\.[^.]*/);

  # set file names
  $gffout_ok_file = $path.$outfile.$ext;
}
my $gffout_ok = prepare_gffout($config, $gffout_ok_file);

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
    dual_print($log, "All attributes will be used !\n");
    $attHashOk{"all_attributes"}++;
  }
  else{
    my @attList= split(/,/, $attributes);
    foreach my $attribute (@attList){
      if($attribute == 0){
        $attHashOk{$attribute}++;
        push(@attListOk, $attribute);
      }
    }
  }
  dual_print($log, "\n");
}

# start with some interesting information
my $stringPrint = strftime "%m/%d/%Y at %Hh%Mm%Ss", localtime;
$stringPrint .= "\nusage: $0 @copyARGV\n";
$stringPrint .= "The attributes @attListOk from the following feature types: $print_feature_string_copy will be copy pasted to the following feature types: $print_feature_string_paste.\n";
dual_print($log, $stringPrint);

                          #######################
# >>>>>>>>>>>>>>>>>>>>>>>>#        MAIN         #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                          #######################
my %all_cases = ('l1' => 0, 'l2' => 0, 'l3' => 0, 'all' => 0);
######################
### Parse GFF input #
my ($hash_omniscient, $hash_mRNAGeneLink) =  slurp_gff3_file_JD({ input => $opt_gff,
                                                                  config => $config
                                                                });
dual_print($log, "Parsing Finished\n");
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
      my $pastel1 = check_feature($feature_l1, 'level1', \@ptagListPaste);

		  #################
		  # == LEVEL 2 == #
		  #################
		  foreach my $tag_l2 (sort keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...

		    if ( exists_keys( $hash_omniscient, ('level2', $tag_l2, $id_l1) ) ){
			    my @list_fl2 = @{$hash_omniscient->{'level2'}{$tag_l2}{$id_l1}};
	        foreach my $feature_l2 ( @list_fl2 ) {

            my $copyl2 = check_feature($feature_l2,'level2', \@ptagListCopy);
            my $pastel2 = check_feature($feature_l2,'level2', \@ptagListPaste);
            #copy attributes from L1 to l2
            if( $copyl1 and $pastel2 ){
              add_tags_from_to($feature_l1, $feature_l2);
            }
            #copy attributes from L2 to l1
            if($pastel1 and $copyl2){
              add_tags_from_to($feature_l2, $feature_l1);
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
                  #copy attributes from L3 to l1
                  if ($copyl3 and $pastel1){
                    add_tags_from_to($feature_l3, $feature_l1);
                  }
                  #copy attributes from L2 to l3
                  if ($copyl2 and $pastel3){
                    add_tags_from_to($feature_l2, $feature_l3);
                  }
                  #copy attributes from L3 to l2
                  if ($copyl3 and $pastel2){
                    add_tags_from_to($feature_l3, $feature_l2);
                  }
                  # ------- Case L3 to L3 -------
                  foreach my $tag_l3_again (sort keys %{$hash_omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...
				            if ( exists_keys( $hash_omniscient, ('level3', $tag_l3_again, $id_l2) ) ){
				       	      my @list_fl3_again = @{$hash_omniscient->{'level3'}{$tag_l3_again}{$id_l2}};
				              foreach my $feature_l3_again ( @list_fl3_again ) {

                        # Check it is not the same feature
                        if( lc($feature_l3->_tag_value('ID')) ne lc($feature_l3_again->_tag_value('ID')) ){

                          my $pastel3_again = check_feature($feature_l3_again, 'level3', \@ptagListPaste);
                          #copy attributes from L3 to l3
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
            foreach my $tag_l2_again (sort keys %{$hash_omniscient->{'level2'}}){ 

              if ( exists_keys( $hash_omniscient, ('level2', $tag_l2_again, $id_l1) ) ){
                my @list_fl2_again = @{$hash_omniscient->{'level2'}{$tag_l2_again}{$id_l1}};

                foreach my $feature_l2_again ( @list_fl2_again ) {
                  # Check it is not the same feature
                  if( lc($feature_l2->_tag_value('ID')) ne lc($feature_l2_again->_tag_value('ID')) ){

                    my $pastel2_again = check_feature($feature_l2_again, 'level3', \@ptagListPaste);
                    #copy attributes from L2 to l2
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
      foreach my $tag_l1_again (sort {$a cmp $b} keys %{$hash_omniscient->{'level1'}}){
        foreach my $feature_l1_again ( @{$hash_sortBySeq->{$seqid}{$tag_l1_again}} ){
          my $id_l1_again = lc($feature_l1_again->_tag_value('ID'));
          if( lc($id_l1) ne lc($id_l1_again) ){
            my $pastel1_again = check_feature($feature_l1_again, 'level1', \@ptagListPaste);
            #copy attributes from L1 to l1
            if ($copyl1 and $pastel1_again){
              add_tags_from_to($feature_l1, $feature_l1_again);
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


=item B<-a>, B<--tag>, B<--att> or B<--attribute>

Attribute that will be copied and pasted. Case sensitive.
You can specified an attribute (or a coma separated list) by giving its attribute tag value (column9) as: Ontology, Dbxref, etc
Default: all_attributes
/!\ <all_attributes> is a specific parameter meaning all the attributes will be use.


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
