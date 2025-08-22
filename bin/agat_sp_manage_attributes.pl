#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Clone 'clone';
use List::MoreUtils qw(uniq);
use AGAT::AGAT;

my $header = get_agat_header();
my ( $opt, $usage, $config ) = AGAT::AGAT::describe_script_options( $header,
    [ 'gff|f=s',      'Input reference gff file', { required => 1 } ],
    [ 'add!',         'Add attribute(s)' ],
    [ 'overwrite!',   'Overwrite attribute(s)' ],
    [ 'cp!',          'Copy attribute value(s)' ],
    [ 'p|type|l=s',   'Primary tag/level to operate on' ],
    [ 'tag|att=s',    'Attribute tag list' ],
);

my $gff        = $opt->gff;
my $add        = $opt->add;
my $overwrite  = $opt->overwrite;
my $cp         = $opt->cp;
my $primaryTag = $opt->p;
my $attributes = $opt->tag;

my $outfile = $config->{output};
my $log;
if ( my $log_name = $config->{log_path} ) {
    open( $log, '>', $log_name )
      or die "Can not open $log_name for printing: $!";
    dual_print( $log, $header, 0 );
}
my $verbose = $config->{verbose};

my $gffout = prepare_gffout( $config, $outfile );

# Manage $primaryTag
my @ptagList;
if ( !$primaryTag or $primaryTag eq "all" ) {
    dual_print( $log, "We will work on attributes from all features\n", $verbose );
    push( @ptagList, "all" );
}
elsif ( $primaryTag =~ /^level[123]$/ ) {
    dual_print( $log, "We will work on attributes from all the $primaryTag features\n", $verbose );
    push( @ptagList, $primaryTag );
}
else {
    @ptagList = split( /,/, $primaryTag );
    foreach my $tag (@ptagList) {
        if ( $tag =~ /^level[123]$/ ) {
            dual_print( $log, "We will work on attributes from all the $tag features\n", $verbose );
        }
        else {
            dual_print( $log, "We will work on attributes from $tag feature.\n", $verbose );
        }
    }
}

# Manage attributes if given
### If attributes given, parse them:
my %attListOk;
my @attListPair;
if ($attributes){

  if ($attributes eq "all_attributes"){
    if($add){
      dual_print( $log, "You cannot use the all_attributes value with the add option. Please change the parameters !\n", $verbose );
      exit;
    }
    dual_print( $log, "All attributes will be removed except ID and Parent attributes !\n", $verbose );
    $attListOk{"all_attributes"}++;
  }
  else{
    @attListPair= split(/,/, $attributes);
    my $nbAtt=$#attListPair+1;

    foreach my $attributeTuple (@attListPair){
      my @attList= split(/\//, $attributeTuple);
      if($#attList == 0){ # Attribute alone
        #check for ID attribute
        if(lc($attList[0]) eq "id" and ! $add){
            dual_print( $log, "It's forbidden to remove the ID attribute in a gff3 file !\n", $verbose );
            exit;
        }
        #check for Parent attribute
        if(lc($attList[0]) eq "parent" and ! $add){
          foreach my $tag (@ptagList){
            if($tag ne "gene" and $tag ne "level1"){
              dual_print( $log, "It's forbidden to remove the $attList[0] attribute to a $tag feature in a gff3 file !\n", $verbose );
              exit;
            }
          }
        }
        $attListOk{$attList[0]}="null";
        if($add){
          dual_print( $log, "$attList[0] attribute will be added. The value will be empty.\n", $verbose );
        }
        else{
          dual_print( $log, "$attList[0] attribute will be removed.\n", $verbose );
        }
      }
      else{ # Attribute will be replaced/copied with a new tag name
        $attListOk{$attList[0]}=$attList[1];
        dual_print( $log, "$attList[0] attribute will be replaced by $attList[1].\n", $verbose );
      }
    }
  }
  dual_print( $log, "\n", $verbose );
}


                #####################
                #     MAIN          #
                #####################


######################
### Parse GFF input #
my ($hash_omniscient, $hash_mRNAGeneLink) = slurp_gff3_file_JD({ input => $gff,
                                                                 config => $config
                                                              });
dual_print( $log, "GFF3 file parsed\n", $verbose );


foreach my $tag_l1 (keys %{$hash_omniscient->{'level1'}}){
  foreach my $id_l1 (keys %{$hash_omniscient->{'level1'}{$tag_l1}}){

    my $feature_l1=$hash_omniscient->{'level1'}{$tag_l1}{$id_l1};

    manage_attributes($feature_l1, 'level1', \@ptagList,\%attListOk);

    #################
    # == LEVEL 2 == #
    #################
    foreach my $tag_l2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...

      if ( exists ($hash_omniscient->{'level2'}{$tag_l2}{$id_l1} ) ){
        foreach my $feature_l2 ( @{$hash_omniscient->{'level2'}{$tag_l2}{$id_l1}}) {

          manage_attributes($feature_l2,'level2', \@ptagList,\%attListOk);
          #################
          # == LEVEL 3 == #
          #################
          my $level2_ID = lc($feature_l2->_tag_value('ID'));

          foreach my $tag_l3 (keys %{$hash_omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...
            if ( exists ($hash_omniscient->{'level3'}{$tag_l3}{$level2_ID} ) ){
              foreach my $feature_l3 ( @{$hash_omniscient->{'level3'}{$tag_l3}{$level2_ID}}) {
                manage_attributes($feature_l3, 'level3', \@ptagList,\%attListOk);
              }
            }
          }
        }
      }
    }
  }
}
#print "We added $nbNameAdded Name attributes\n";

print_omniscient( {omniscient => $hash_omniscient, output => $gffout} );

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

sub  manage_attributes{
  my  ($feature, $level, $ptagList, $attListOk)=@_;

  my $primary_tag=$feature->primary_tag;

  # check primary tag (feature type) to handle
  foreach my $ptag (@$ptagList){

    if($ptag eq "all"){
      remove_tag_from_list($feature,$attListOk);
    }
    elsif(lc($ptag) eq $level){
      remove_tag_from_list($feature,$attListOk);
    }
    elsif(lc($ptag) eq lc($primary_tag) ){
      remove_tag_from_list($feature,$attListOk);
    }
  }
}

sub remove_tag_from_list{
  my  ($feature, $attListOk)=@_;

  if (exists ($attListOk{"all_attributes"} ) ){ # all attributes removed except ID and Parent
    my @list_att = $feature->get_all_tags;
    foreach my $tag (@list_att){
      if(lc($tag) ne "id" and lc($tag) ne "parent"){
        $feature->remove_tag($tag);
      }
    }
  }
  else{
    foreach my $att (keys %{$attListOk}){

      if ($feature->has_tag($att)){

        if ($attListOk{$att} eq "null" ){ # the attribute name is kept inctact
					if(!$add){
						$feature->remove_tag($att);
					}
					elsif($add and $overwrite){
						create_or_replace_tag($feature,$att,'empty');
					}
				}
        else{ # We replace the attribute name

          my @values=$feature->get_tag_values($att);
          my $newAttributeName=$attListOk{$att};
					if($overwrite){
            create_or_replace_tag($feature,$newAttributeName, @values);
					}
					else{#if new attribute exist we do no overwrite it
						if(! $feature->has_tag($newAttributeName)){
							create_or_replace_tag($feature,$newAttributeName, @values);
						}
						#else we do not change the existing value
					}
					if(! $cp){
            $feature->remove_tag($att); #remove old attribute if it is not the cp option
          }
        }
      }

      elsif($add){
        if ($attListOk{$att} eq "null" ){ # the attribute name is kept inctact
          create_or_replace_tag($feature,$att,'empty');
        }
      }
    }
  }
}


__END__

=head1 NAME

agat_sp_manage_attributes.pl

=head1 DESCRIPTION

The script removes choosen attributes of selected features. It can also create new
attribute with 'empty' value, or copy paste an existing attribute using a new specified tag.
Attribute in a gff file have this shape (2 attributes here): tag=value;tag=value and
are stored within the 9th column.

=head1 SYNOPSIS

    agat_sp_manage_attributes.pl -gff file.gff  -att locus_tag,product,name/NewName -p level2,cds,exon [ -o outfile ]
    agat_sp_manage_attributes.pl --help

=head1 OPTIONS

=over 8

=item B<--gff> or B<-f>

Input GTF/GFF file.

=item B<-p>,  B<--type> or  B<-l>

primary tag option, case insensitive, list. Allow to specied the feature types that will be handled.
You can specified a specific feature by given its primary tag name (column 3) as: cds, Gene, MrNa
You can specify directly all the feature of a particular level:
      level2=mRNA,ncRNA,tRNA,etc
      level3=CDS,exon,UTR,etc
By default all feature are taking in account. fill the option by the value "all" will have the same behaviour.

=item B<--tag>, B<--att>

Attributes with the tag specified will be removed from the feature type specified by the option p (primary tag). List of tag must be coma separated.
/!\\ You must use "" if name contains spaces.
Instead to remove an attribute, you can replace its Tag by a new Tag using this formulation tagName/newTagName.
To remove all attributes non mandatory (only ID and Parent are mandatory) you can use the option with <all_attributes> parameter.

=item B<--add>

Attribute with the tag specified will be added if doesn't exist. The value will be 'empty'.

=item B<--cp>

Bolean. When tags specied are with this form: tagName/newTagName.
By using this <cp> parameter, the attribute with the tag tagName will be duplicated
with the new tag newTagName if no attribute with the tag newTagName already exits.

=item B<--overwrite>

When using -add parameter, if an attribute with the specificed tag already exists, it will not be modified.
When using --cp parameter, if an attribute with the specificed newTagName already exists, it will not be modified.
So using the --overwrite parameter allows to overwrite the value of the existing attribute.

=item B<-o> , B<--output> , B<--out> or B<--outfile>

Output GFF file.  If no output file is specified, the output will be
written to STDOUT.

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
