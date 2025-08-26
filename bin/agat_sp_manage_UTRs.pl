#!/usr/bin/env perl

use strict;
use warnings;
use POSIX qw(strftime);
use Carp;
use Try::Tiny;
use File::Basename;
use IO::File;
use Pod::Usage;
use Getopt::Long::Descriptive;
use AGAT::AGAT;

my $header = get_agat_header();
my $config;
my $DefaultUTRnb=5;

my @copyARGV = @ARGV;
my ( $opt, $usage, $cfg ) = AGAT::AGAT::describe_script_options( $header,
    [ 'gff|f|ref|reffile=s', 'Input GTF/GFF file', { required => 1 } ],
    [ 'number|n|t|nb=i',     'Threshold of exon\'s number of the UTR' ],
    [ 'mode' => hidden => { one_of => [
            [ 'three|3|three_prime_utr' => 'Apply threshold on the 3\'UTR' ],
            [ 'five|5|five_prime_utr'   => 'Apply threshold on the 5\'UTR' ],
            [ 'both|b|bs'               => 'Apply threshold on both UTRs' ],
        ] } ],
    [ 'p|plot!', 'Allows to create an histogram in pdf of UTR sizes distribution' ],
);

my $opt_reffile = $opt->gff;
my $opt_nbUTR   = $opt->number;
my $mode        = $opt->mode;
my $opt_plot    = $opt->{plot};
my $opt_output  = $opt->out;
$config         = $cfg;
my $opt_verbose = $config->{verbose};

my $log;
if ( my $log_name = $config->{log_path} ) {
  open( $log, '>', $log_name ) or die "Can not open $log_name for printing: $!";
}
dual_print( $log, $header);

my $opt_utr3 = ($mode && $mode eq 'three') ? 1 : 0;
my $opt_utr5 = ($mode && $mode eq 'five') ? 1 : 0;
my $opt_bst  = ($mode && $mode eq 'both') ? 1 : 0;

if ( ! defined($opt_reffile ) or ! ($opt_utr3 or $opt_utr5 or $opt_bst or $opt_plot) ) {
    pod2usage( {
           -message => "$header\nMust specify at least 2 parameters:\nReference data gff3 file (--gff)\nOne UTR option (3, 5 , both, plot)",
           -verbose => 0,
           -exitval => 1 } );
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    PARAMS    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
my $ostreamReport_file;
if (defined($opt_output) ) {
  my ($name,$path,$ext) = fileparse($opt_output,qr/\.[^.]*/);
  $opt_output = $path.$name;
  if (-d $opt_output){
    my $msg = "The output directory choosen already exists. Please geve me another Name.\n";
    dual_print($log, $msg, 3);
    warn $msg if $opt_verbose;
    exit();
  }
  else{
    mkdir $opt_output;
  }
  $ostreamReport_file = $opt_output."/report.txt";
}

my $ostreamReport = prepare_fileout($ostreamReport_file);

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    EXTRA     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

my $string1 = strftime "%m/%d/%Y at %Hh%Mm%Ss", localtime;
$string1 .= "\n\nusage: $0 @copyARGV\n\n";

if (! $opt_nbUTR){
  $opt_nbUTR=$DefaultUTRnb;
}elsif(!($opt_utr3 or $opt_utr5 or $opt_bst)){$string1 .= "The value $opt_nbUTR of the parameter <n> will no be taken into account. Indeed no UTRs option called. (three, five, both).\n";}
if($opt_utr3 or $opt_utr5 or $opt_bst){
  $string1 .= "Genes with more than $opt_nbUTR UTRs will be reported.\n";
}

print $ostreamReport $string1;
dual_print($log, $string1, $opt_output ? $opt_verbose : 3);

# Check if dependencies for plot are available
if($opt_plot){
	if ( ! may_i_plot() ) {
		$opt_plot = undef;
	}
}
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    PARAMS    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# #################################
# # Manage Ouput Directory / File #
# #################################

my $ostreamUTR_file;
my $ostreamUTRdiscarded_file;

if (defined($opt_output) ) {

  my ($file_in,$path,$ext) = fileparse($opt_reffile,qr/\.[^.]*/);

  #manage name output
  my $utr_type_under=undef;
  my $utr_type_over=undef;
  # case no filter so we don't create discarded file output.
  if (! ($opt_utr3 or $opt_utr5 or $opt_bst)){

    $utr_type_under = $file_in;
    $utr_type_over = $file_in;
    $ostreamUTR_file = $opt_output."/".$file_in.".gff";
  }
  else{  # case with filter so we create discarded file output and a of file output.
    if ($opt_utr3){
      $utr_type_under = $file_in."_UTR3_under".$opt_nbUTR;
      $utr_type_over = $file_in."_UTR3_overORequal".$opt_nbUTR;
    }
    if ($opt_utr5){
      if($utr_type_under){
         $utr_type_under.="_and_UTR5_under".$opt_nbUTR;
         $utr_type_over.="_and_UTR5_overORequal".$opt_nbUTR;
      }
      else{
        $utr_type_under=$file_in."_UTR5_under".$opt_nbUTR;
        $utr_type_over=$file_in."_UTR5_overORequal".$opt_nbUTR;
      }
    }
    if ($opt_bst){
      if($utr_type_under){
         $utr_type_under.="_and_bothSides_under".$opt_nbUTR;
         $utr_type_over.="_and_bothSides_overORequal".$opt_nbUTR;
      }
      else{
        $utr_type_under=$file_in."_bothSides_under".$opt_nbUTR;
        $utr_type_over=$file_in."_bothSides_overORequal".$opt_nbUTR;
      }
    }
    $ostreamUTR_file = $opt_output."/".$utr_type_under.".gff";
    $ostreamUTRdiscarded_file = $opt_output."/".$utr_type_over.".gff";
  }
}

my $ostreamUTR = prepare_gffout($config, $ostreamUTR_file);
my $ostreamUTRdiscarded = prepare_gffout($config, $ostreamUTRdiscarded_file);

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>     MAIN     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

######################
### Parse GFF input #
my ($hash_omniscient, $hash_mRNAGeneLink) = slurp_gff3_file_JD({ input => $opt_reffile,
                                                                 config => $config
                                                              });
dual_print($log, "Parsing Finished\n\n");
### END Parse GFF input #
#########################

my %UTRdistribution;
my %UTRbymRNA;
my %UTRoverview;
# #########################################################
# # A.1) Assign utr side if they are not
# #########################################################

foreach my $tag_l3 ( keys %{$hash_omniscient->{'level3'}} ) {
   if($tag_l3 =~ /utr/){
     if ($tag_l3 ne 'three_prime_utr' and $tag_l3 ne 'five_prime_utr') {

       foreach my $id_l2 ( keys %{$hash_omniscient->{'level3'}{$tag_l3}} ){

         my $geneID = $hash_mRNAGeneLink->{$id_l2};
         my $feature_l2 = get_feature_l2_from_id_l2_l1($hash_omniscient, $id_l2, $geneID);
         my $strand = $feature_l2->strand();
         my $cds_feature_example =  $hash_omniscient->{'level3'}{'cds'}{$id_l2}[0]; #if utr exists, cds should exists

         foreach my $feature_l3 ( @{$hash_omniscient->{'level3'}{$tag_l3}{$id_l2}}){
           if($feature_l3->start <= $cds_feature_example->start){
             if ($strand eq "+" or $strand eq "1"){
               $feature_l3->primary_tag= "five_prime_utr";
             }
             else{$feature_l3->primary_tag= "three_prime_utr";}
           }
           else{
             if ($strand eq "+" or $strand eq "1"){
               $feature_l3->primary_tag= "three_prime_utr";
             }
             else{
               $feature_l3->primary_tag= "five_prime_utr";
             }
           }
         }
       }
     }
   }
}

# #########################################################
# # A.1) Count utr exon by side and total
# #########################################################

foreach my $tag_l2 (keys %{$hash_omniscient->{'level2'}}) {
  foreach my $id_l1 (keys %{$hash_omniscient->{'level2'}{$tag_l2}}) {
    foreach my $feature_l2 (@{$hash_omniscient->{'level2'}{$tag_l2}{$id_l1}}) {

      my $has_an_utr=undef;
      my $id_l2= lc($feature_l2->_tag_value('ID'));

      foreach my $tag_l3 (keys %{$hash_omniscient->{'level3'}}) {
        if($tag_l3 =~ /utr/){
          if(exists ($hash_omniscient->{'level3'}{$tag_l3}{$id_l2})){
            $has_an_utr="yes";

            my $nbUTR = $#{$hash_omniscient->{'level3'}{$tag_l3}{$id_l2}}+1; #nb utrs
            $UTRdistribution{$tag_l3}{$nbUTR}++;
            $UTRbymRNA{$tag_l3}{$id_l2}=$nbUTR;
            $UTRoverview{$tag_l3}++;
            if(! exists $UTRbymRNA{'both'}{$id_l2}){
              $UTRbymRNA{'both'}{$id_l2}=$nbUTR;
            }
            else{
              $UTRbymRNA{'both'}{$id_l2}+=$nbUTR;
            }
          }
        }
      }
      if(!$has_an_utr){
         $UTRbymRNA{'both'}{$id_l2}=0;
      }
    }
  }
}

###########################
# compute if on of UTR option called
###########################

###########################
# Overview of UTRs per sides
if($opt_utr3 or $opt_utr5 or $opt_bst){
  # print preliminary results
  my $stringPrint="";
  foreach my $key (keys %UTRoverview) {
    $stringPrint.="There are ".scalar $UTRoverview{$key}." $key\n";
    my $total=0;
    foreach my $value  ( sort {$b <=> $a} keys %{$UTRdistribution{$key}}){
      if($value >= $opt_nbUTR){
        $total+=$UTRdistribution{$key}{$value};
      }
      else{last;}
    }
    $stringPrint.= "Among them $total have over or equal $opt_nbUTR exons.\n";
  }

  ###########################
  # Overview of UTRs both
  # print preliminary results
  $stringPrint.="There are ".scalar %{$UTRbymRNA{'both'}}." features that have UTRs (some at both sides some only at one extremity)\n";
  my $total=0;
  foreach my $id  ( keys %{$UTRbymRNA{'both'}}){
    if($UTRbymRNA{'both'}{$id} >= $opt_nbUTR){
      $total++;
    }
  }
  $stringPrint.= "Among them $total have over or equal UTR (5' and/or 3') exons.\n\n";

  ###########################
  # Main compute
  my @listIDl2discarded;
  my @listIDlok;
  my %geneName;
  my %geneName_ok;
  foreach my $tag (keys %UTRbymRNA) {
    foreach my $id_level2 (keys %{$UTRbymRNA{$tag}}){

      # case only opt_utr3
      if ($opt_utr3 and $tag eq "three_prime_utr"){
        if ($UTRbymRNA{$tag}{$id_level2} >= $opt_nbUTR){
          push @listIDl2discarded, $id_level2 ;
          $geneName{$hash_mRNAGeneLink->{$id_level2}}++;
        }
        else{
          push @listIDlok, $id_level2 ;
          $geneName_ok{$hash_mRNAGeneLink->{$id_level2}}++;
        }
      }
      # case only opt_utr5
      if ($opt_utr5 and $tag eq "five_prime_utr"){
        if ($UTRbymRNA{$tag}{$id_level2} >=  $opt_nbUTR){
          push @listIDl2discarded, $id_level2 ;
          $geneName{$hash_mRNAGeneLink->{$id_level2}}++;
        }
        else{
          push  @listIDlok, $id_level2;
          $geneName_ok{$hash_mRNAGeneLink->{$id_level2}}++;
        }
      }                                           ### REMOVE OPTION BOTH ?
      # case both side together (when added) should be over $opt_nbUTR)
      if ($opt_bst and $tag eq "both"){
        if ($UTRbymRNA{$tag}{$id_level2} >=  $opt_nbUTR){
          push @listIDl2discarded, $id_level2;
          $geneName{$hash_mRNAGeneLink->{$id_level2}}++;
        }
        else{
          push @listIDlok, $id_level2 ;
          $geneName_ok{$hash_mRNAGeneLink->{$id_level2}}++;
        }
      }
      # case both side independant (side 3 and and5 should be over $opt_nbUTR)

      # case no option print all so put all in @listIDlok
      if(! $opt_utr3 and ! $opt_utr5 and ! $opt_bst) { # in case where no option, We do by default side3 side5 idependant. On sufficiant to discard the mRNA
          push @listIDlok, $id_level2 ;
          $geneName_ok{$hash_mRNAGeneLink->{$id_level2}}++;
      }
    }
  }

  # remove duplicate in case several option tends to give the same case
  if(@listIDl2discarded){
    my $sizeList= @listIDl2discarded;
    my $nbGene = keys %geneName;
    $stringPrint.= "According to the parameters $sizeList RNA discarded from $nbGene genes\n";
    my @listIDl2discardedUniq = uniq(@listIDl2discarded);
    my $omniscient_discarded = create_omniscient_from_idlevel2list($hash_omniscient, $hash_mRNAGeneLink, \@listIDl2discarded);
    print_omniscient( {omniscient => $omniscient_discarded, output => $ostreamUTRdiscarded} );

  }
  if(@listIDlok){
    my $sizeList= @listIDlok;
    my $nbGeneOk = keys %geneName_ok;
    $stringPrint.= "$sizeList RNA from $nbGeneOk genes pass the filter (under the UTR Threshold).\n";
    my @listIDlokUniq = uniq(@listIDlok);
    my $omniscient_ok = create_omniscient_from_idlevel2list($hash_omniscient, $hash_mRNAGeneLink, \@listIDlokUniq);
    print_omniscient( {omniscient => $omniscient_ok, output => $ostreamUTR} );
  }
  if(@listIDl2discarded and @listIDlok){
    my $union=0;
    foreach my $name (keys %geneName){
      if (exists ($geneName_ok{$name}) ){
        $union++;
      }
    }
     $stringPrint.= "$union genes have RNA isoform discarded (over the UTR Threshold).\n";
  }

  #Print Info OUtput
  print $ostreamReport $stringPrint;
  dual_print($log, $stringPrint, $opt_output ? $opt_verbose : 3);
}

############################
# Plot distribution if asked
if ($opt_plot){

  foreach my $utr_type (keys %UTRdistribution) {

    my $txtFile;
    my $outPlot;
    my $txtFileOver;
    my $outPlotOver;
    if($opt_output){
      if($opt_utr3 or $opt_utr5 or $opt_bst){
        $txtFileOver = $opt_output."/".$utr_type."_overORequal".$opt_nbUTR.".txt";
        $outPlotOver = $opt_output."/".$utr_type."_overORequal".$opt_nbUTR.".pdf";
        $txtFile = $opt_output."/".$utr_type."_under".$opt_nbUTR.".txt";
        $outPlot = $opt_output."/".$utr_type."_under".$opt_nbUTR.".pdf";
      }
      else{
        $txtFile = $opt_output."/".$utr_type.".txt";
        $outPlot = $opt_output."/".$utr_type.".pdf";
      }
    }else{
      $txtFile = $utr_type.".txt";
      $outPlot = $utr_type.".pdf";
      if($opt_utr3 or $opt_utr5 or $opt_bst){
        $txtFileOver = $utr_type."_over".$opt_nbUTR."txt";
        $outPlotOver = $utr_type."_over".$opt_nbUTR."pdf";
      }
    }
    #print file thtat will be read by R
    open(FH, ">".$txtFile) || die "Erreur E/S:$!\n";
    if($opt_utr3 or $opt_utr5 or $opt_bst){
      open(FH_filter, ">".$txtFileOver) || die "Erreur E/S:$!\n";
    }
    my $firstLine="yes";
    my $firstLineOver="yes";
    foreach my $value (keys %{$UTRdistribution{$utr_type}}) {

      if($opt_utr3 or $opt_utr5 or $opt_bst){ #we have a filter
        if($value >= $opt_nbUTR){ #print utr over threshold
          if($firstLineOver){
            print FH_filter $value."\t".$UTRdistribution{$utr_type}{$value};
            $firstLineOver=undef;
          }
          else{
            print FH_filter "\n".$value."\t".$UTRdistribution{$utr_type}{$value};
          }
        }
        else{ #print utr under threshold
          if($firstLine){
            print FH $value."\t".$UTRdistribution{$utr_type}{$value};
            $firstLine=undef;
          }
          else{
            print FH "\n".$value."\t".$UTRdistribution{$utr_type}{$value};
          }
        }

      }
      else{ # no filter we print everything
        if($firstLine){
          print FH $value."\t".$UTRdistribution{$utr_type}{$value};
          $firstLine=undef;
        }
        else{
          print FH "\n".$value."\t".$UTRdistribution{$utr_type}{$value};
        }
      }
    }
    close FH;

  my $R = Statistics::R->new() or die "Problem with R : $!\n";

  #R command
  $R->send(
        qq`
        listValues=as.matrix(read.table("$txtFile", sep="\t", he=F)) ##///!!!\\\\\
        legendToDisplay=paste("Number of value used : ",length(listValues))
        listValueMoreThan <- listValues[listValues[,1]>5,]

        pdf("$outPlot")
        plot(listValues[,2]~listValues[,1], xlab="Contig size", ylab="Frequency", main="Size distribution of $utr_type")
        dev.off()

        pdf("$outPlotOver")
        plot(listValueMoreThan[,2]~listValueMoreThan[,1], xlab="Contig size", ylab="Frequency", main="Size distribution of $utr_type over 5")
        dev.off()`
            );

  # Close the bridge
  $R->stopR();

  # Delete temporary file
  unlink "$txtFile";
  if($opt_utr3 or $opt_utr5 or $opt_bst){
    unlink "$txtFileOver";
  }
  }

}

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

sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}

__END__


=head1 NAME

agat_sp_manage_UTRs.pl

=head1 DESCRIPTION

Detect the genes containing too much UTR's exon according to a choosen threshold.
If no UTR option (3, 5, 3 and 5, both) is given the threshold will be not used.
option 3 and 5 together is different of "both". In the first case the gene is discarded if either the 3' or the 5' UTR contains more exon than the threshold given.
In the second case, will be discarded only the genes where the addition of UTR's exon of both side is over the threshold given.

=head1 SYNOPSIS

    agat_sp_manage_UTRs.pl --ref infile --three --five -p --out outFile
    agat_sp_manage_UTRs.pl --help

=head1 OPTIONS

=over 8

=item B<--gff>, B<--ref>, B<--reffile> or B<-f>

Input GTF/GFF file.

=item B<-n>, B<-t>, B<--nb> or B<--number>

Threshold of exon's number of the UTR. Over or equal to this threshold, the UTR will be discarded. Default value is 5.

=item B<-3>, B<--three> or B<--tree_prime_utr>

The threshold of the option <n> will be applied on the 3'UTR.

=item B<-5>, B<--five> or B<--five_prime_utr>

The threshold of the option <n> will be applied on the 5'UTR.

=item B<-b>, B<--both> or B<--bs>

The threshold of the option <n> will be applied on genes where the number of UTR exon (3' and 5' additioned) is over it.

=item  B<--p> or B<--plot>

Allows to create an histogram in pdf of UTR sizes distribution.

=item  B<--out>, B<--output> or B<-o>

Output gff3 file where the gene incriminated will be write.

=item B<-c> or B<--config>

String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any, 
otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

=item B<--help> or B<-h>

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
