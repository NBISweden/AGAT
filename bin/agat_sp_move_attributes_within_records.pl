#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use POSIX qw(strftime);
use Pod::Usage;
use IO::File;
use AGAT::AGAT;

start_script();
my $header = get_agat_header();
# -------------------------------- LOAD OPTIONS --------------------------------
my $primaryTagCopy="level2";
my $primaryTagPaste="level3";
my $opt_output= undef;
my $attributes="all_attributes";
my $opt_gff = undef;
my $opt_help;

# OPTION MANAGEMENT: split shared vs script options
my ($shared_argv, $script_argv) = split_argv_shared_vs_script(\@ARGV);

# Parse script-specific options
my $script_parser = Getopt::Long::Parser->new;
$script_parser->configure('bundling','no_auto_abbrev');
if ( ! $script_parser->getoptionsfromarray(
    $script_argv,
    'f|ref|reffile|gff=s'  => \$opt_gff,
    'feature_copy|fc=s'    => \$primaryTagCopy,
    'feature_paste|fp=s'   => \$primaryTagPaste,
    'o|output=s'           => \$opt_output,
    'a|tag|att|attribute=s'=> \$attributes,
    'h|help!'              => \$opt_help,
  ) )
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
my $shared_opts = parse_shared_options($shared_argv);
initialize_agat({ config_file_in => ( $shared_opts->{config} ), input => $opt_gff, shared_opts => $shared_opts });

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
    dual_print1 "All attributes will be used !\n";
    $attHashOk{"all_attributes"}++;
  }
  else{
    my @attList= split(/,/, $attributes);

    foreach my $attribute (@attList){

      if($attribute == 0){ # Attribute alone
        #check for ID attribute
        if(lc($attribute) eq "id" ){ die "ID attribute cannot be modified !\n"; }
        #check for Parent attribute
        if(lc($attribute) eq "parent"){ die "Parent attribute cannot be modified !\n"; }
        $attHashOk{$attribute}++;
        push(@attListOk, $attribute);
      }
    }
  }
  dual_print1 "\n";
}

# start with some interesting information
dual_print1 "The attributes @attListOk from the following feature types: $print_feature_string_copy will be copy pasted to the following feature types: $print_feature_string_paste.\n";

                          #######################
# >>>>>>>>>>>>>>>>>>>>>>>>#        MAIN         #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                          #######################
my %all_cases = ('l1' => 0, 'l2' => 0, 'l3' => 0, 'all' => 0);
######################
### Parse GFF input #
my ($hash_omniscient) =  slurp_gff3_file_JD({ input => $opt_gff });

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
		  foreach my $tag_l2 (sort keys %{$hash_omniscient->{'level2'}}){ 

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

				    foreach my $tag_l3 (sort keys %{$hash_omniscient->{'level3'}}){ 
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
                  foreach my $tag_l3_again (sort keys %{$hash_omniscient->{'level3'}}){ 
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

=item B<-h> or B<--help>

Display this helpful text.

=back

=head1 SHARED OPTIONS

Shared options are defined in the AGAT configuration file and can be overridden via the command line for this script only.
Common shared options are listed below; for the full list, please refer to the AGAT agat_config.yaml.

=over 8

=item B<--config>

String - Path to a custom AGAT configuration file.  
By default, AGAT uses `agat_config.yaml` from the working directory if present, otherwise the default file shipped with AGAT
(available locally via `agat config --expose`).

=item B<--cpu>, B<--core>, B<--job> or B<--thread>

Integer - Number of parallel processes to use for file input parsing (via forking).

=item B<-v> or B<--verbose>

Integer - Verbosity, choice are 0,1,2,3,4. 0 is quiet, 1 is normal, 2,3,4 is more verbose. Default 1.

=back

=head1 FEEDBACK

For questions, suggestions, or general discussions about AGAT, please use the AGAT community forum:
https://github.com/NBISweden/AGAT/discussions

=head1 BUG REPORTING

Bug reports should be submitted through the AGAT GitHub issue tracker:
https://github.com/NBISweden/AGAT/issues

=cut

AUTHOR - Jacques Dainat
