#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Clone 'clone';
use Getopt::Long::Descriptive;
use Pod::Usage;
use List::MoreUtils qw(uniq);
use AGAT::AGAT;

my $header    = get_agat_header();
my $start_run = time();

my ( $opt, $usage, $config ) = AGAT::AGAT::describe_script_options(
    $header,
    [ 'gff|f=s', 'Input GTF/GFF file', { required => 1 } ],
    [ 'tag|att=s', 'Attribute tag(s) to manage', { required => 1 } ],
    [ 'type|p|l=s',
        'Comma-separated list of feature types or level1/2/3',
        { default => 'all' }
    ],
    [ 'add!',       'Add attribute if it does not exist' ],
    [ 'cp!',        'Duplicate attribute using tag/newTag pairs' ],
    [ 'overwrite!', 'Overwrite existing attribute when adding or copying' ],
    [ 'value=s',    'Process only attributes matching this value' ],
    [
        'strategy=s',
        'Strategy for --value (equal or match)',
        {
            default   => 'equal',
            callbacks => {
                valid => sub {
                    $_[0] =~ /^(?:equal|match)$/
                      or die 'Strategy must be equal or match';
                    return 1;
                }
            }
        }
    ],
);

my $gff        = $opt->gff;
my $attributes = $opt->tag;
my $primaryTag = $opt->type;
my $add        = $opt->add;
my $cp         = $opt->cp;
my $overwrite  = $opt->overwrite;
my $value      = $opt->value;
my $strategy   = $opt->strategy;
my $outfile    = $config->{output};
my $cpt_case   = 0;

my $log;
if ( my $log_name = $config->{log_path} ) {
  open( $log, '>', $log_name )
    or die "Can not open $log_name for printing: $!";
  dual_print( $log, $header, 0 );
}

# --- Manage output
my $gffout = prepare_gffout( $config, $outfile );

# deal with strategy input
$strategy=lc($strategy);
if( ($strategy ne "equal") and ($strategy ne "match") ){
  dual_print($log, "Strategy must be <equal> or <match>. Wrong value provided: <$strategy>\n", $config->{verbose});
  exit;
}

# Manage $primaryTag
my @ptagList;
if(! $primaryTag or $primaryTag eq "all"){
  dual_print($log, "We will work on attributes from all features\n", $config->{verbose});
  push(@ptagList, "all");
}elsif($primaryTag =~/^level[123]$/){
  dual_print($log, "We will work on attributes from all the $primaryTag features\n", $config->{verbose});
  push(@ptagList, $primaryTag);
}else{
   @ptagList= split(/,/, $primaryTag);
   foreach my $tag (@ptagList){
      if($tag =~/^level[123]$/){
        dual_print($log, "We will work on attributes from all the $tag features\n", $config->{verbose});
      }
      else{
       dual_print($log, "We will work on attributes from $tag feature.\n", $config->{verbose});
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
      dual_print($log, "You cannot use the all_attributes value with the add option. Please change the parameters !\n", $config->{verbose});
      exit;
    }
    dual_print($log, "All attributes will be removed except ID and Parent attributes !\n", $config->{verbose});
    $attListOk{"all_attributes"}++;
  }
  else{
    @attListPair= split(/,/, $attributes);
    my $nbAtt=$#attListPair+1;

    foreach my $attributeTuple (@attListPair){
      my @attList= split(/\//, $attributeTuple);
      if($#attList == 0){ # Attribute alone
        #check for ID attribute
        if(lc($attList[0]) eq "id" and ! $add){dual_print($log, "It's forbidden to remove the ID attribute in a gff3 file !\n", $config->{verbose}); exit;}
        #check for Parent attribute
        if(lc($attList[0]) eq "parent" and ! $add){
          foreach my $tag (@ptagList){
            if($tag ne "gene" and $tag ne "level1"){
              dual_print($log, "It's forbidden to remove the $attList[0] attribute to a $tag feature in a gff3 file !\n", $config->{verbose});
              exit;
            }
          }
        }
        $attListOk{$attList[0]}="null";
        if($add){
          dual_print($log, "$attList[0] attribute will be added. The value will be empty.\n", $config->{verbose});
        }
        else{
          if($value){
              dual_print($log, "$attList[0] attribute will be removed if it has the value:$value.\n", $config->{verbose});
          }
          else{
            dual_print($log, "$attList[0] attribute will be removed.\n", $config->{verbose});
          }
        }
      }
      else{ # Attribute will be replaced/copied with a new tag name
        $attListOk{$attList[0]}=$attList[1];
        dual_print($log, "$attList[0] attribute will be replaced by $attList[1].\n", $config->{verbose});
      }
    }
  }
  dual_print($log, "\n", $config->{verbose});
}

my $hash_info= get_levels_info();
my $hash_level = $hash_info->{'other'}{'level'};


                #####################
                #     MAIN          #
                #####################


# Manage gff file
my $format = $config->{force_gff_input_version};
if(! $format ){ $format = select_gff_format($gff); }
my $ref_in = AGAT::BioperlGFF->new(-file => $gff, -gff_version => $format);

#time to calcul progression
my $startP=time;
my $nbLine=`wc -l < $gff`;
$nbLine =~ s/ //g;
chomp $nbLine;
dual_print($log, "$nbLine line to process...\n", $config->{verbose});

my $line_cpt=0;
my %hash_IDs;
my %featCount;
my %mapID;
while (my $feature = $ref_in->next_feature() ) {
  $line_cpt++;

  manage_attributes($feature, \@ptagList,\%attListOk);
  $gffout->write_feature($feature);

  #####################
  #Display progression
  if ((30 - (time - $startP)) < 0) {
    my $done = ($line_cpt*100)/$nbLine;
    $done = sprintf ('%.0f', $done);
        dual_print($log, "\rProgression : $done % processed.\n", $config->{verbose});
    $startP= time;
  }
}

##Last round
my $end_run = time();
my $run_time = $end_run - $start_run;

if($add){
  dual_print($log, "$cpt_case attribute added\n", $config->{verbose});
}
elsif($cp){
  dual_print($log, "$cpt_case attribute copied\n", $config->{verbose});

}else{
  dual_print($log, "$cpt_case attribute removed\n", $config->{verbose});
}
dual_print($log, "Job done in $run_time seconds\n", $config->{verbose});


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
  my  ($feature, $ptagList, $attListOk)=@_;

  my $primary_tag=$feature->primary_tag;

  # check primary tag (feature type) to handle
  foreach my $ptag (@$ptagList){

    if($ptag eq "all"){
      remove_tag_from_list($feature,$attListOk);
    }
    elsif( exists_keys($hash_level,(lc($ptag),lc($primary_tag) ) ) ){
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
        if(!$value or ( $strategy eq "equal" and lc($value) eq lc($feature->_tag_value($tag) ) )
        or ( $strategy eq "match" and lc($feature->_tag_value($tag) ) =~ lc($value) ) ){
          $feature->remove_tag($tag);
          $cpt_case++;
        }
      }
    }
  }

  else{
    foreach my $att (keys %{$attListOk}){

      if ($feature->has_tag($att)){

        if(!$value or ( $strategy eq "equal" and lc($value) eq lc($feature->_tag_value($att) ) )
        or ( $strategy eq "match" and lc($feature->_tag_value($att) ) =~ lc($value) ) ){

          if ($attListOk{$att} eq "null" ){ # the attribute name is kept inctact
  					if(!$add){
  						  $feature->remove_tag($att);
                $cpt_case++;
  					}
  					elsif($add and $overwrite){
  						create_or_replace_tag($feature,$att,'empty');
              $cpt_case++;
  					}
  				}
          else{ # We replace the attribute name

            my @values=$feature->get_tag_values($att);
            my $newAttributeName=$attListOk{$att};
  					if($overwrite){
              create_or_replace_tag($feature,$newAttributeName, @values);
              $cpt_case++;
  					}
  					else{#if new attribute exist we do no overwrite it
  						if(! $feature->has_tag($newAttributeName)){
  							create_or_replace_tag($feature,$newAttributeName, @values);
                $cpt_case++;
  						}
  						#else we do not change the existing value
  					}
  					if(! $cp){
              $feature->remove_tag($att); #remove old attribute if it is not the cp option
              $cpt_case++;
            }
          }
        }

        elsif($add){
          if ($attListOk{$att} eq "null" ){ # the attribute name is kept inctact
            create_or_replace_tag($feature,$att,'empty');
            $cpt_case++;
          }
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

    agat_sq_manage_attributes.pl --gff file.gff  --att locus_tag,product,name/NewName -p level2,cds,exon [ -o outfile ]
    agat_sq_manage_attributes.pl --help

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
By default all feature are taking in account.

=item B<--tag>, B<--att>

Attributes with the tag specified will be removed from the feature type specified by the option p (primary tag). List of tag must be coma separated.
/!\\ You must use "" if name contains spaces.
Instead to remove an attribute, you can replace its Tag by a new Tag using this formulation tagName/newTagName.
To remove all attributes non mandatory (only ID and Parent are mandatory) you can use the option with <all_attributes> parameter.

=item B<--add>

Attribute with the tag specified will be added if doesn't exist. The value will be 'empty'.

=item B<--cp>

When tags specied are with this form: tagName/newTagName.
By using this <cp> parameter, the attribute with the tag tagName will be duplicated
with the new tag newTagName if no attribute with the tag newTagName already exits.

=item B<--overwrite>

When using --add parameter, if an attribute with the specificed tag already exists, it will not be modified.
When using --cp parameter, if an attribute with the specificed newTagName already exists, it will not be modified.
So using the --overwrite parameter allows to overwrite the value of the existing attribute.

=item B<--value>

String. When a value is provided the attribute is taken into account only if
the attribute contains (or match) a specific value

=item B<--strategy>

String. [Default equal]. Strategy to use when --value parameter is in use. Can be equal or match.
Equal => the attribute value must be identical. Match => the attribute value must match

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
