#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Clone 'clone';
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use IO::File;
use List::MoreUtils qw(uniq);
use AGAT::AGAT;

my $header = get_agat_header();
my $config;
my %handlers;
my $gff = undef;
my $one_tsv = undef;
my $opt_help= undef;
my $primaryTag=undef;
my $attributes=undef;
my $outfile=undef;
my $outInOne=undef;
my $doNotReportEmptyCase=undef;

if ( !GetOptions(
    'c|config=s'               => \$config,
    "h|help" => \$opt_help,
    "gff|f=s" => \$gff,
    "d!" => \$doNotReportEmptyCase,
    "m|merge!" => \$one_tsv,
    "p|t|l=s" => \$primaryTag,
    "attribute|a|att=s" => \$attributes,
    "output|outfile|out|o=s" => \$outfile))

{
    pod2usage( { -message => 'Failed to parse command line',
                 -verbose => 1,
                 -exitval => 1 } );
}

# Print Help and exit
if ($opt_help) {
    pod2usage( { -verbose => 99,
                 -exitval => 0,
                 -message => "$header\n" } );
}

if ( ! $gff or ! $attributes ){
    pod2usage( {
           -message => "$header\nAt least 2 parameter is mandatory:\nInput reference gff file (--gff)\n".
           "Attribute tag to investigate --att \n\n",
           -verbose => 0,
           -exitval => 2 } );
}

# --- Manage config ---
$config = get_agat_config({config_file_in => $config});

# If one output file we can create it here
my $outfile_pref; my $path ; my $ext;
if ($outfile) {
    ($outfile_pref,$path,$ext) = fileparse($outfile,qr/\.[^.]*/);
}

if($one_tsv){
	$outInOne = prepare_fileout($outfile);
}

# Manage $primaryTag
my @ptagList;
if(! $primaryTag or $primaryTag eq "all"){
  print "We will work on attributes from all features\n";
  push(@ptagList, "all");
}elsif($primaryTag =~/^level[123]$/){
  print "We will work on attributes from all the $primaryTag features\n";
  push(@ptagList, $primaryTag);
}else{
   @ptagList= split(/,/, $primaryTag);
   foreach my $tag (@ptagList){
      if($tag =~/^level[123]$/){
        print "We will work on attributes from all the $tag features\n";
      }
      else{
       print "We will work on attributes from $tag feature.\n";
      }
   }
}

# Manage attributes if given
### If attributes given, parse them:
my @attListOk;
if ($attributes){
  my @attList = split(/,/, $attributes); # split at comma as separated value

  foreach my $attribute (@attList){
      push @attListOk, $attribute;
      print "$attribute attribute will be processed.\n";

  }
  print "\n";
}


                #####################
                #     MAIN          #
                #####################


######################
### Parse GFF input #
my ($hash_omniscient, $hash_mRNAGeneLink) = slurp_gff3_file_JD({ input => $gff,
                                                                 config => $config
                                                              });
print ("GFF3 file parsed\n");


foreach my $tag_l1 (sort keys %{$hash_omniscient->{'level1'}}){
  foreach my $id_l1 (sort keys %{$hash_omniscient->{'level1'}{$tag_l1}}){

    my $feature_l1=$hash_omniscient->{'level1'}{$tag_l1}{$id_l1};

    manage_attributes($feature_l1, 'level1', \@ptagList,\@attListOk);

    #################
    # == LEVEL 2 == #
    #################
    foreach my $tag_l2 (sort keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...

      if ( exists ($hash_omniscient->{'level2'}{$tag_l2}{$id_l1} ) ){
        foreach my $feature_l2 ( @{$hash_omniscient->{'level2'}{$tag_l2}{$id_l1}}) {

          manage_attributes($feature_l2,'level2',, \@ptagList,\@attListOk);
          #################
          # == LEVEL 3 == #
          #################
          my $level2_ID = lc($feature_l2->_tag_value('ID'));

          foreach my $tag_l3 (sort keys %{$hash_omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...
            if ( exists ($hash_omniscient->{'level3'}{$tag_l3}{$level2_ID} ) ){
              foreach my $feature_l3 ( @{$hash_omniscient->{'level3'}{$tag_l3}{$level2_ID}}) {
                manage_attributes($feature_l3, 'level3', \@ptagList,\@attListOk);
              }
            }
          }
        }
      }
    }
  }
}
#print "We added $nbNameAdded Name attributes\n";


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
      tag_from_list($feature,$attListOk);
    }
    elsif(lc($ptag) eq $level){
      tag_from_list($feature,$attListOk);
    }
    elsif(lc($ptag) eq lc($primary_tag) ){
      tag_from_list($feature,$attListOk);
    }
  }
}

sub tag_from_list{
  my  ($feature, $attListOk)=@_;

  my $tags_string = undef;
  foreach my $att ( @{$attListOk} ){

    # create handler if needed (on the fly)
    if (! $one_tsv){
      if(! exists ( $handlers{$att} ) ) {
        my $out = IO::File->new();
        if ($outfile) {
          my $file_name =  $path.$outfile_pref."_".$att.$ext;
          open($out, '>', $file_name) or die "Could not open file $file_name $!";
        }
        else{
          $out->fdopen( fileno(STDOUT), 'w' );
        }
        $handlers{$att}=$out;
      }
    }


    if ($feature->has_tag($att)){

      # get values of the attribute
      my @values = $feature->get_tag_values($att);

      # print values of one attribute per file
      if (! $one_tsv){
        my $out = $handlers{$att};
        print $out join(",", @values), "\n";
      }
      else{ # put everything in one tsv
        $tags_string .= join(",", @values)."\t";
      }
    }
    else{
      if (! $one_tsv){
        my $out = $handlers{$att};
        print $out ".\n" if (! $doNotReportEmptyCase);
      }
      else{ # put everything in one tsv
        if (! $doNotReportEmptyCase){
          $tags_string .= ".\t";
        }
        else{
          $tags_string .= "\t";
        }
      }
    }
  }
  if($tags_string){
    chop $tags_string;
    print $outInOne  $tags_string."\n";
  }
}


__END__

=head1 NAME

agat_sp_extract_attributes.pl

=head1 DESCRIPTION

The script takes a gtf/gff file as input.
The script allows to extract choosen attributes of all or specific feature types.
The 9th column of a gff/gtf file contains a list of attributes.
An attribute (gff3) looks like that tag=value

=head1 SYNOPSIS

    agat_sp_extract_attributes.pl --gff file.gff  --att locus_tag,product,name -p level2,cds,exon [ -o outfile ]
    agat_sp_extract_attributes.pl --help

=head1 OPTIONS

=over 8

=item B<--gff> or B<-f>

Input GTF/GFF file.

=item B<-p>,  B<-t> or  B<-l>

primary tag option, case insensitive, list. Allow to specied the feature types that will be handled.
You can specified a specific feature by given its primary tag name (column 3) as: cds, Gene, MrNa
You can specify directly all the feature of a particular level:
      level2=mRNA,ncRNA,tRNA,etc
      level3=CDS,exon,UTR,etc
By default all feature are taking in account. fill the option by the value "all" will have the same behaviour.

=item B<--attribute>, B<--att>, B<-a>

attribute tag. The value of the attribute tag specified will be extracted from the feature type specified by the option -p. List of attributes must be coma separated.

=item B<--merge> or B<-m>

By default the values of each attribute tag is writen in its dedicated file. To write the values of all tags in only one file use this option.

=item B<-d>

By default when an attribute is not found for a feature, a dot (.) is reported. If you don't want anything to be printed in such case use this option.

=item B<-o> , B<--output> , B<--out> or B<--outfile>

Output GFF file.  If no output file is specified, the output will be
written to STDOUT.

=item B<-c> or B<--config>

String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any, 
otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
The --config option gives yo the possibility to use your own AGAT config file (located elsewhere or named differently).

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
