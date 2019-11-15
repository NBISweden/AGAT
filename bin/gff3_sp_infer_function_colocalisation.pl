#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Statistics::R;
use Pod::Usage;
use Bio::Tools::GFF;
use List::MoreUtils qw(uniq);
use Agat::Omniscient;

my $header = get_agat_header();
my $output = undef;
my $ref = undef;
my $tar = undef;
my $inv = undef;
my $_dblr = undef;
my $_ablt= undef;
my $liftOver=undef;
my $overlapT = undef;
my $featureType=undef;
my $checkMultiOverlap=undef;
my $help= 0;

if ( !GetOptions(
    "help|h" => \$help,
    "ref|r=s" => \$ref,
    "mapped|m|tar=s" => \$tar,
    "inv|inverse|oppposite" => \$inv,
    "value|threshold|overlap=i" => \$overlapT,
    "t|transfert|lift" => \$liftOver,
    "feature=s" => \$featureType,
    "dblr" => \$_dblr,
    "ablt" => \$_ablt,
    "cmo" => \$checkMultiOverlap,
    "outdir|out|o=s" => \$output))

{
    pod2usage( { -message => "$header"."Failed to parse command line.",
                 -verbose => 1,
                 -exitval => 1} );
}

# Print Help and exit
if ($help) {
    pod2usage( { -message => "$header",
                 -verbose => 2,
                 -exitval => 2 } );
}

if ( ! (defined($ref)) or  ! (defined($tar)) ){
    pod2usage( {
           -message => "$header\nAt least 2 parameters are mandatory.\n",
           -verbose => 0,
           -exitval => 1 } );
}

######################
# Manage output file #
my($filename_tar, $dirs, $suffix) = fileparse($tar,qr/\.[^.]*/);
my($filename_ref, $dirs, $suffix) = fileparse($ref,qr/\.[^.]*/);

my $gffout;
my $outReport;

if ($output) {
  $output=~ s/.gff//g;
  if (-d $output ){
    print "Directory $output already exists.\n";exit;
  }
  else{
    mkdir $output;

    my $out=">".$output."/".$filename_tar."_FuncLiftOn_".$filename_ref.".gff";
    open(my $fh, $out) or die "Could not open file '$out' $!";
    $gffout= Bio::Tools::GFF->new(-fh => $fh, -gff_version => 3 );

    my $out=$output."/".$filename_tar."_FuncLiftOn_".$filename_ref."-report.txt";
    open($outReport, '>', $out) or die "Could not open file '$out' $!";
  }
}
else{
  $gffout = Bio::Tools::GFF->new(-fh => \*STDOUT, -gff_version => 3);
}

##################
# Manage overlap threshold
if (! $overlapT ){
  $overlapT = 20;
}
print "We will consider gene with overlap over $overlapT\n";

##################
# Manage feature Type (Level2)
if(! $featureType){
  $featureType='mRNA';
}
$featureType=lc($featureType);
#####################################
# END Manage Options                #
#####################################

###########
# DEFINE NBIS PATTERN OF NAMES
###########
my $nbis_suffix_p=qr/_([0-9]*(_iso[0-9]+)?$|iso[0-9]+$)/o;
my $nbis_suffix_d=qr/_partial_part.*/;

                #####################
                #     MAIN          #
                #####################
############
# Parse GFF reference #
print ("Parse file $ref\n");
my ($refhash_omniscient, $refhash_mRNAGeneLink) = slurp_gff3_file_JD({ input => $ref
                                                              });
print ("$ref file parsed\n");

##############
# Manage gene name reference
my %hash_geneName;
my %hash_geneNameLab;
if(! $_dblr){
  # save gene names
  foreach my $tag_level1 (keys %{$refhash_omniscient->{'level1'}}){
    foreach my $geneID (keys %{$refhash_omniscient->{'level1'}{$tag_level1}} ) {
      my $gene_feature = $refhash_omniscient->{'level1'}{$tag_level1}{$geneID};
      my $tag=undef;
      if($gene_feature->has_tag('Name')){
        $tag='Name';
      }
      elsif($gene_feature->has_tag('gene_name')){
        $tag='gene_name';
      }
      if($tag){
        my @tmp=$gene_feature->get_tag_values($tag);
        my $name=lc(shift @tmp);
        $hash_geneNameLab{$name}++;
        $name =~ s/$nbis_suffix_p//; # We remove what has been added by NBIS during gene name annotation
        $hash_geneName{$name}++;
      }
    }
  }
}

############
# Parse GFF target #
print ("Parse file $tar\n");
my ($tarhash_omniscient, $tarhash_mRNAGeneLink) = slurp_gff3_file_JD({ input => $tar
                                                              });
print ("$tar file parsed\n");

#count level1 in reference
my $count_ref_level1=0;
foreach my $tag_level1 (keys %{$refhash_omniscient->{'level1'}}) {
    my $nbKeys= keys  %{$refhash_omniscient->{'level1'}{$tag_level1}};
    $count_ref_level1=$count_ref_level1+ $nbKeys;
}

#count level1 in tar
my $count_tar_level1=0;
foreach my $tag_level1 (keys %{$tarhash_omniscient->{'level1'}}) {
    my $nbKeys= keys  %{$tarhash_omniscient->{'level1'}{$tag_level1}};
    $count_tar_level1=$count_tar_level1+ $nbKeys;
}

#############
# sort by seq id and featuretype
print ("Sort files by seqid\n");
my $refhash_sortBySeq = sort_by_seq($refhash_omniscient, $featureType) ;
my $tarhash_sortBySeq = sort_by_seq($tarhash_omniscient, $featureType) ;
print ("Sorting files terminated\n");

#count level1 in ref sort by seq
my $count_ref_level1_sortBySeq=0;
foreach my $tag_level1 (keys %{$refhash_sortBySeq}) {
  foreach my $Contig (keys %{$refhash_sortBySeq->{$tag_level1}}) {
    my $nbKeys= @{$refhash_sortBySeq->{$tag_level1}{$Contig}};
    $count_ref_level1_sortBySeq=$count_ref_level1_sortBySeq + $nbKeys;
  }
}

#count level1 in tar sort by seq
my $count_tar_level1_sortBySeq=0;
foreach my $tag_level1 (keys %{$tarhash_sortBySeq}) {
  foreach my $Contig (keys %{$tarhash_sortBySeq->{$tag_level1}}) {
    my $nbKeys= @{$tarhash_sortBySeq->{$tag_level1}{$Contig}};
    $count_tar_level1_sortBySeq=$count_tar_level1_sortBySeq + $nbKeys;
  }
}

print "Now we are analyzing the contigs containing annotations from both annotation builds:\n";
my $OverlapingB=0;
my $OverlapingA=0;
my $nbNoOverlapingA=0;
my $nbTotalClusterCaseA=0;
my $nbTotalClusterCaseB=0;
my %clusterCase;
my %o_one2one;
my %o_one2many;
my %o_many2one;
my %o_many2many;
my %split_omniscient;
my %fusion_omniscient;

######
# LOOP over all the contig
foreach my $tag_level1 (keys %{$refhash_sortBySeq}) {
  foreach my $ContigName (keys %{$refhash_sortBySeq->{$tag_level1}}) {

    ######
    # LOOP over all gene feature
    # USE A COPY - use of of temporary variable to be sure to loop over all element.
    my @copyContig = @{$refhash_sortBySeq->{$tag_level1}{$ContigName}}; #
    foreach my $copyGene (@copyContig) {

      my @copyGeneNameList=$copyGene->get_tag_values('ID');
      my $GeneName=lc(shift @copyGeneNameList);
      #print "\nStudy of ".$GeneName.":\n";

      ###### Test if gene already studied (started by another gene but due to overlap it is already studied)
      my $letStudyThatGene=undef;
      foreach my $gene_feature (@{$refhash_sortBySeq->{$tag_level1}{$ContigName}}){
        my @GeneNameList=$gene_feature->get_tag_values('ID');
        my $geneNameOriginal=lc(shift @GeneNameList);
        if ($GeneName eq $geneNameOriginal){$letStudyThatGene="yes";} # already studied if we cannot find it among refhash_sortBySeq
      }
      ##### End test if already studied. If not, we continue


      if ($letStudyThatGene){

        # Declare table which I will work with
        my @ListRefOverlapAtotest_list;my @ListOverlapAtested_list;my @ListNoOverlapA_list; my @ListOverlapAtestnoneed_list;
        my @ListRefOverlapBtotest_list;my @ListOverlapBtested_list;my @ListNoOverlapB_list; my @ListOverlapBtestnoneed_list;

        my $ListBtotest=\@ListRefOverlapBtotest_list; my $ListOverlapAtested=\@ListOverlapAtested_list; my $ListNoOverlapA=\@ListNoOverlapA_list;
                                                                my $ListOverlapBtested=\@ListOverlapBtested_list; my $ListNoOverlapB=\@ListNoOverlapB_list;

        #### Initialize list of gene to test
        my @LinkTocurrentGeneFeature;
        push (@LinkTocurrentGeneFeature, $copyGene);  #### >>>>>> BASE OF THE FEATURE TESTED FOR OVERLAP !!! CURRENTLY WE CHECK THE GENE FEATURE !
        my $ListAtotest=\@LinkTocurrentGeneFeature;


        my $lap=0;
        while (@{$ListAtotest} != 0 or @{$ListBtotest} != 0 ){
            $lap++;

            #########
            # TEST A side
            if (@{$ListAtotest} != 0){
                #print "START START START TEST AAAAAAAA\n";
                my ($ListToTestX, $ListOverlapTested, $ListNoOverlapTested) = retrieveAllOverlap( $lap, $ContigName,
                                                                              $refhash_sortBySeq, $tarhash_sortBySeq, $ListAtotest, $ListOverlapAtested, $ListNoOverlapA,1,
                                                                               $refhash_omniscient, $tarhash_omniscient, $featureType);

                $ListBtotest = $ListToTestX;
                $ListOverlapAtested = $ListOverlapTested;
                $ListNoOverlapA = $ListNoOverlapTested;

                # Reinitialisation empty
                my @list_empty;
                $ListAtotest = \@list_empty;

                #print "END END END TEST AAAAAAAAA\n";
                next(); #stop here and avoid test B
            }

            ###########
            # TEST B side
            if (@{$ListBtotest} != 0){
                # test every B
                #print "START TEST BBBBBBBBBBBBBBBB\n";
                my ($ListToTestX, $ListOverlapTested, $ListNoOverlapTested) = retrieveAllOverlap( $lap, $ContigName,
                                                                                $tarhash_sortBySeq, $refhash_sortBySeq, $ListBtotest, $ListOverlapBtested,$ListNoOverlapB,2,
                                                                                 $tarhash_omniscient, $refhash_omniscient, $featureType);
                $ListAtotest = $ListToTestX;
                $ListOverlapBtested = $ListOverlapTested;
                $ListNoOverlapB = $ListNoOverlapTested;

                # Reinitialisation empty
                my @list_empty;
                $ListBtotest =\@list_empty;

                #print "END TEST BBBBBBBBBBBBBBBBBB\n";
                next(); #stop here
            }
        }

        $nbNoOverlapingA=$nbNoOverlapingA+@{$ListNoOverlapA};
        my $nbOverlpInA=@{$ListOverlapAtested};
        my @OverlapInA=(@{$ListOverlapAtested});
        $OverlapingA+=$nbOverlpInA;

        my $nbOverlpInB=@{$ListOverlapBtested};
        my @OverlapInB=(@{$ListOverlapBtested});
        $OverlapingB+=$nbOverlpInB;

        my $nbFragment = $nbOverlpInA+$nbOverlpInB;

    ##################################
    # Manage different cases         #
    ##################################
    #print ("OVERLAP step 1 END frgt nb:$nbFragment\n");
        if ( $nbFragment < 2 ){
           #HEre can be printed => Single gene from file A and PerfectMatch from A or B depending to an iption like my $contigFusionA = "ok";  < /!\ >
        }
        elsif ( $nbFragment == 2){ #case can be stretched
          @{$o_one2one{$GeneName}} = (@OverlapInA,@OverlapInB);
        }
        else{ # $nbFragment > 2)
          if (@OverlapInA == 1){
            @{$o_one2many{$GeneName}} = (@OverlapInA,[@OverlapInB]);
          }
          elsif(@OverlapInB == 1){
            @{ $o_many2one{$GeneName}} = ([@OverlapInA],@OverlapInB);
          }  # Split in the target Build
          else{ #

           @{ $o_many2many{$GeneName}} = ([@OverlapInA],[@OverlapInB]);
            $nbTotalClusterCaseA=$nbTotalClusterCaseA+$nbOverlpInA;
            $nbTotalClusterCaseB=$nbTotalClusterCaseB+$nbOverlpInB;
          }
        }
      }
    }
  }
}
print "Now managing overlaping cases found\n";
##############################
#Manage overlap one2one
my $nb_one2one = keys (%o_one2one);
#print "We have $nb_one2one cases one 2 one\n\n";
my $nbLifted_one2one=0;
my $nbNameChanged_one2one=0;
my $nbNewName_one2one=0;
my $overlapOK=0;

if ($liftOver){
  ($overlapOK, $nbNameChanged_one2one, $nbNewName_one2one) = manage_one2one(\%o_one2one,$refhash_omniscient,$tarhash_omniscient,$overlapT, \%hash_geneName);
}
$nbLifted_one2one=$nbNameChanged_one2one + $nbNewName_one2one;

##############################
#Manage overlap one2many
my $nb_one2many=0;
my $nb_fusion_ok=0;
my $nb_fusion_notV=0;
my $nbLifted_one2many=0;
my $nbNameChanged_one2many=0;
my $nbNewName_one2many=0;
my $overlapOK_one2many=0;
my %tmp_o_one2one;

$nb_one2many = keys (%o_one2many);
#print "We have $nb_one2many cases one2many\n\n";
#Check if same gene that is split in 2
foreach my $key (keys %o_one2many){

  my $geneA_feature = $o_one2many{$key}[0];
  my @geneB = @{$o_one2many{$key}[1]};

#  print $geneA_feature->gff_string."\n";
#  print $geneB[0]->gff_string."\n";
#  print $geneB[1]->gff_string."\n\n";

  # prepare in case we have to save the feature
  my $geneA_name = lc($geneA_feature->_tag_value('ID'));
  my @IDlist_A=($geneA_name);

  my $sameGene=check_feature_same_names(\%hash_geneName, \@geneB, 'Bs');

  if($sameGene eq "none"){# GENE DIFFERENT in B (cannot really compare because one is missing) - only one long in A
    $nb_fusion_notV++;
    fill_omniscient_from_other_omniscient_level1_id(\@IDlist_A, $refhash_omniscient, \%fusion_omniscient);
  }
  elsif(! $sameGene){ # GENE DIFFERENT in B - only one long in A
      $nb_fusion_ok++;
      fill_omniscient_from_other_omniscient_level1_id(\@IDlist_A, $refhash_omniscient, \%fusion_omniscient);
  }
  elsif($liftOver){ # Gene are the same / Option lift given / we try to lift name from one gene choosen randomly
    @{$tmp_o_one2one{$geneA_name}} = (($geneA_feature),($geneB[0]));
  }


  ### force try to change name for better tuple
  if(((! $sameGene) or ($sameGene eq "none") ) and ($checkMultiOverlap) and ($liftOver)){

    my $overlap_percent_1 = testOverlaplevel3($geneA_feature, $geneB[0], $refhash_omniscient, $tarhash_omniscient, $featureType);
    my $overlap_percent_2 = testOverlaplevel3($geneA_feature, $geneB[1], $refhash_omniscient, $tarhash_omniscient, $featureType);

    ## We don't consider if both have 100% overlap
    my $bestFeature = undef;
    if ($overlap_percent_1 > $overlap_percent_2){
      $bestFeature=$geneB[0];
    }
    elsif($overlap_percent_1 < $overlap_percent_2){
      $bestFeature=$geneB[1];
    }

    if ($bestFeature){
      @{$tmp_o_one2one{$geneA_name}} = (($geneA_feature),($bestFeature));
    }
  }
}
($overlapOK_one2many, $nbNameChanged_one2many, $nbNewName_one2many) = manage_one2one(\%tmp_o_one2one,$refhash_omniscient,$tarhash_omniscient,$overlapT, \%hash_geneName);
$nbLifted_one2many = $nbNameChanged_one2many + $nbNewName_one2many;

##############################
#Manage overlap many2one
my $nb_many2one=0;
my $nb_split_ok=0;
my $nb_split_notV=0;
my $nbLifted_many2one=0;
my $nbNameChanged_many2one=0;
my $nbNewName_many2one=0;
my $overlapOK_many2one=0;
my %tmp_o_one2one;

$nb_many2one = keys (%o_many2one);
#print "We have $nb_many2one cases many2one\n";

foreach my $key (keys %o_many2one){

  my @geneA = @{$o_many2one{$key}[0]};
  my $geneB_feature = $o_many2one{$key}[1];

  my $sameGene=check_feature_same_names(\%hash_geneName, \@geneA, 'As');

  # prepare in case we have to save the feature
  my $geneA0_name=lc($geneA[0]->_tag_value('ID'));
  my $geneA1_name=lc($geneA[1]->_tag_value('ID'));
  my @IDlist=($geneA0_name, $geneA1_name);

  if(! $sameGene){ # NAME DIFFERENT in A - only one long in B
      #print "Looks like $geneA0_name and $geneA1_name should be only one gene which has been split.\n";
      fill_omniscient_from_other_omniscient_level1_id(\@IDlist, $refhash_omniscient, \%split_omniscient);
  }
  elsif($sameGene eq "none"){ # NAME DIFFERENT in A (cannot really compare because one is missing) - only one long in B
    $nb_split_notV++;
    fill_omniscient_from_other_omniscient_level1_id(\@IDlist, $refhash_omniscient, \%split_omniscient);
  }
  else{$nb_split_ok++;
    if ($liftOver) { #same gene we take one randomly to use the name
      @{$tmp_o_one2one{$geneA[0]}} = ($geneA[0],$geneB_feature);
      @{$tmp_o_one2one{$geneA[1]}} = ($geneA[1],$geneB_feature);
    }
  }

   ### force try to change name for better tuple
  if(((! $sameGene) or ($sameGene eq "none") ) and ($checkMultiOverlap) and ($liftOver)){
    @{$tmp_o_one2one{$geneA[0]}} = (($geneA[0]),($geneB_feature));
    @{$tmp_o_one2one{$geneA[1]}} = (($geneA[1]),($geneB_feature));

      #print $geneA[0]->gff_string."\n";
      #print $geneA[1]->gff_string."\n";
      #print $geneB_feature->gff_string."\n\n";
      #print $bestFeature->gff_string."\n";
      #print "$overlap_percent_1 $overlap_percent_2 \n\n";
  }
}

($overlapOK_many2one, $nbNameChanged_many2one, $nbNewName_many2one) = manage_one2one(\%tmp_o_one2one,$refhash_omniscient,$tarhash_omniscient,$overlapT, \%hash_geneName);
$nbLifted_many2one=$nbNameChanged_many2one, $nbNewName_many2one;

##############################
#Manage overlap many2many
my $nb_many2many = 0;
my $nbLifted_many2many=0;
my $nbNameChanged_many2many=0;
my $nbNewName_many2many=0;
my $overlapOK_many2many=0;
my %tmp_o_one2one;

$nb_many2many = keys (%o_many2many);

if( ($checkMultiOverlap) and ($liftOver) ){
  foreach my $key (keys %o_many2many){
    my @List_geneA_feature = @{$o_many2many{$key}[0]};
    my @List_geneB_feature = @{$o_many2many{$key}[1]};

    foreach my $featureA (@List_geneA_feature){
      my $best_overlap_featureB=undef;
      my $best_overlap_percent=-1;
      foreach my $featureB (@List_geneB_feature){
        if(defined(testOverlaplevel3($featureA, $featureB, $refhash_omniscient, $tarhash_omniscient, $featureType))){
          my $overlap_percent = testOverlaplevel3($featureA, $featureB, $refhash_omniscient, $tarhash_omniscient, $featureType);

          if($overlap_percent > $best_overlap_percent){
            $best_overlap_percent=$overlap_percent;
            $best_overlap_featureB=$featureB;
          }
        }
      }

      #Save the best overlaping gene B to Manage downstream
      if($best_overlap_featureB){
        @{$tmp_o_one2one{lc($featureA->_tag_value('ID'))}} = (($featureA),($best_overlap_featureB));
      }
    }
  }
}

( $overlapOK_many2many, $nbNameChanged_many2many, $nbNewName_many2many ) = manage_one2one(\%tmp_o_one2one,$refhash_omniscient,$tarhash_omniscient,$overlapT, \%hash_geneName);
$nbLifted_many2many=$nbNameChanged_many2many + $nbNewName_many2many;

#################
# Display results
my $totalA=$OverlapingA+$nbNoOverlapingA;

my $resultToPrint.= "\n\n######### RESULTS #########:\n\n";
$resultToPrint.=  "File1 ($ref):\nTotal gene=$count_ref_level1\n".
                  "Total gene Sort by seq $count_ref_level1_sortBySeq\n".
                  "Total gene checked : $totalA \n".
                  "   => $OverlapingA genes overlap genes from file 2\n".
                  "   => $nbNoOverlapingA genes seem to be \"orphan\"\n\n";

$resultToPrint.=  "File2 ($tar):\nTotal gene = $count_tar_level1\nTotal gene Sort by seq $count_tar_level1_sortBySeq\n";
$resultToPrint.=  "Total gene checked = $OverlapingB\n";
$resultToPrint.=  "=> $OverlapingB genes overlap genes from file 1\n\n";

$resultToPrint.=  "Further details:\n================\n\n";
$resultToPrint.=  "$nb_one2one one2one cases. (A one2one case is 1 gene from $ref overlaping 1 genes of $tar. \n";
if($liftOver){
  $resultToPrint.=  "Among the one2one cases $overlapOK genes overlap correctly. We lifted annotation for $nbLifted_one2one of these cases. Among them:\nWe modified previous annotation for $nbNameChanged_one2one case(s).\n".
  "We newly annotated $nbNewName_one2one previously unannotated case(s).\n";

}

$resultToPrint.=  "\n$nb_one2many one2many (putative fusion) cases. (A fusion case is 1   gene  from $ref overlaping >=2 genes of $tar. \n".
                  "  => $nb_fusion_ok of them are clearly right and reported in output (If we suppose that Target is true)($tar genes have identical names).\n".
                  "  => $nb_fusion_notV of them cannot be checked (do not have any name).\n";
if($checkMultiOverlap or $liftOver){
  $resultToPrint.= "=> Among the one2many cases $overlapOK_one2many genes overlap correctly.  We lifted annotation for $nbLifted_one2many of these cases. Among them:\nWe modified previous annotation(s) for $nbNameChanged_one2many case(s).\n".
  "We newly annotated $nbNewName_one2many previously unannotated case(s).\n";
}

$resultToPrint.=  "\n$nb_many2one many2one (putative split) cases. (A split  case is >=2 genes from $ref overlaping 1   gene  of $tar. \n".
                  "  => $nb_split_ok of them are clearly right and reported in output ($ref genes have identical names).\n".
                  "  => $nb_split_notV of them cannot be checked (do not have any name)\n"; #The others have to be checked because can be close duplicate/paralog.\n\n";
if($checkMultiOverlap or $liftOver){
  $resultToPrint.= "=> Among the many2one cases $overlapOK_many2one genes overlap correctly.  We lifted annotation for $nbLifted_many2one of these cases. Among them:\nWe modified previous annotation(s) for $nbNameChanged_many2one case(s).\n".
  "We newly annotated $nbNewName_many2one previously unannotated case(s).\n";
}

$resultToPrint.=  "\nNumber many2many : $nb_many2many\n";
$resultToPrint.=  "Among the many2many cases there is a total of $nbTotalClusterCaseA genes from $ref and $nbTotalClusterCaseB from $tar\n";
if($checkMultiOverlap){
  $resultToPrint.= "=> Among the many2many cases $overlapOK_many2many genes overlap correctly.  We lifted annotation for $nbLifted_many2many of these cases. Among them:\nWe modified previous annotation(s) for $nbNameChanged_many2many case(s).\n".
  "We newly annotated $nbNewName_many2many previously unannotated case(s).\n";
}

$resultToPrint.= "\n!! Split and fusion results don't take in account possible rearrangement that could occur independently during the evolution of each genome investigated.\n".
                 "!! Please consider that Split, Fusion and Cluster cases an also be mostly symply independant gene overlaping (Common in prokaryotes). <= Don't forget that using the option \"all\" you compare any kind of features\n";


$resultToPrint.="\nCommand line executed:";
foreach (@ARGV) { $resultToPrint.= "$_ " };$resultToPrint.="\n";
## OUTPUT
if($output){

  $resultToPrint.=  "=> Result are written in gff3 format in $output directory\n\n";
  print $outReport $resultToPrint;

  if($nb_one2one != 0){
    print "results with one2one cases\n";
    print_omniscient($refhash_omniscient, $gffout);
  }
  if($nb_one2many != 0) {
    print "one2many cases\n";
    open(my $fh, '>', $output."/".$filename_tar."_FuncLiftOn_".$filename_ref."_fusion.gff") or die "Could not open file '$output' $!";
    my $gffout_fusion= Bio::Tools::GFF->new(-fh => $fh, -gff_version  => 3 );
    print_omniscient(\%fusion_omniscient, $gffout_fusion);
  }
  if($nb_many2one != 0) {
    print "many2one cases\n";
    open(my $fh, '>', $output."/".$filename_tar."_FuncLiftOn_".$filename_ref."_split.gff") or die "Could not open file '$output' $!";
    my $gffout_split= Bio::Tools::GFF->new(-fh => $fh, -gff_version  => 3 );
    print_omniscient(\%split_omniscient, $gffout_split);
  }
}

print $resultToPrint;
#######################################################################################################################
        ####################
         #     METHODS    #
          ################
           ##############
            ############
             ##########
              ########
               ######
                ####
                 ##


sub retrieveAllOverlap {
  my ($lap, $ContigName, $pieceStudiedHashA, $pieceStudiedHashB,  $List_geneRef_to_test, $ListA_tested, $ListA_NoOverlap, $ok, $ref_omniscient, $oppo_omniscient, $featureType) = @_;
  my @newListToTest;
  my @listGeneToRemove;

  foreach my $geneRef_f (@{$List_geneRef_to_test}) {

      my $startA = $geneRef_f->start;
      my $endA = $geneRef_f->end;
      my $Overlaped_FinalAnswer=undef; # default value

      # Now test of all the opposite pieces
      foreach my $tagB (keys %{$pieceStudiedHashB}){ # allow to work only on gene on the contig
        foreach my $geneOppo_f( @{$pieceStudiedHashB->{$tagB}{$ContigName}}){ # allow to work only on gene on the contig

          my $startB=$geneOppo_f->start;
          my $endB=$geneOppo_f->end;

          my $resuOverlap = testOverlap($startA, $endA, $startB, $endB);  ## <====== Call Overlap method

          if ($resuOverlap){

            if(defined(testOverlaplevel3($geneRef_f, $geneOppo_f, $ref_omniscient, $oppo_omniscient, $featureType))){

              my @tmp=$geneOppo_f->get_tag_values('ID');
              my $name_geneOppo=lc(shift @tmp);

              # check if doublon
              my $tmp_list = pushIfNotExit(\@newListToTest, $geneOppo_f);
              @newListToTest = @{$tmp_list};

              $Overlaped_FinalAnswer="yes";
              #print "OVERLAP FOUND => No NeedToverify\n";
            }
          }
        }
      }
      ## END ALL OPPO CHECKED

     if(! $Overlaped_FinalAnswer){
        if ($lap == "1"){
          #print "NO Overlap\n";
          push (@{$ListA_NoOverlap}, $geneRef_f);
        }
        else{
          #print "No More Overlap\n";
          $ListA_tested = pushIfNotExit($ListA_tested, $geneRef_f);
        }
      }
     else{
        $ListA_tested = pushIfNotExit($ListA_tested, $geneRef_f);
      }

      # REMOVE
      # Need to be deleted to not re-use it for retrieve overlap => Because we will found again the same // The removed one will be display through @ListA_Overlap_tested
      removeElementInList($pieceStudiedHashA, $geneRef_f, $ContigName);
  }

  return (\@newListToTest, $ListA_tested, $ListA_NoOverlap);
}

sub  testOverlaplevel3{

  my ($geneA_feature, $geneB_feature, $refhash_omniscient, $tarhash_omniscient, $featureType)=@_;
  my $overlap=undef;

  my $name_geneA=lc($geneA_feature->_tag_value('ID'));
  my $name_geneB=lc($geneB_feature->_tag_value('ID'));

  if($featureType eq "all"){

    foreach my $tag_level2_A (keys %{$refhash_omniscient->{'level2'}}){
      if (exists $refhash_omniscient->{'level2'}{$tag_level2_A}{$name_geneA}){

        foreach my $feature_level2_A (@{$refhash_omniscient->{'level2'}{$tag_level2_A}{$name_geneA}}){
          my $name_feature_level2_A=lc($feature_level2_A->_tag_value('ID'));

          foreach my $tag_level2_B (keys %{$tarhash_omniscient->{'level2'}}){
            if (exists $tarhash_omniscient->{'level2'}{$tag_level2_B}{$name_geneB}){

              foreach my $feature_level2_B (@{$tarhash_omniscient->{'level2'}{$tag_level2_B}{$name_geneB}}){
                my $name_feature_level2_B=lc($feature_level2_B->_tag_value('ID'));

                my $ref_list_to_checkA=undef;
                if(exists($refhash_omniscient->{'level3'}{'exon'}{$name_feature_level2_A})){
                  $ref_list_to_checkA=$refhash_omniscient->{'level3'}{'exon'}{$name_feature_level2_A};
                }
                elsif(exists($refhash_omniscient->{'level3'}{'cds'}{$name_feature_level2_A})){
                  $ref_list_to_checkA=$refhash_omniscient->{'level3'}{'cds'}{$name_feature_level2_A}
                }

                my $ref_list_to_checkB=undef;
                if(exists($tarhash_omniscient->{'level3'}{'exon'}{$name_feature_level2_B})){
                  $ref_list_to_checkB=$tarhash_omniscient->{'level3'}{'exon'}{$name_feature_level2_B};
                }
                elsif(exists($tarhash_omniscient->{'level3'}{'cds'}{$name_feature_level2_B})){
                  $ref_list_to_checkB=$tarhash_omniscient->{'level3'}{'cds'}{$name_feature_level2_B};
                }

                if($ref_list_to_checkA and $ref_list_to_checkB){ # we found at least exon or cds for both
                  if(defined(featuresList_overlap($ref_list_to_checkA, $ref_list_to_checkB))){
                    $overlap=featuresList_overlap($ref_list_to_checkA, $ref_list_to_checkB);
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  else{
    #if (exists_keys($refhash_omniscient, ('level2', $featureType, $name_geneA)) ){
    if (exists $refhash_omniscient->{'level2'}{$featureType}{$name_geneA}){

      foreach my $feature_level2_A (@{$refhash_omniscient->{'level2'}{$featureType}{$name_geneA}}){
        my $name_feature_level2_A=lc($feature_level2_A->_tag_value('ID'));

      #  if (exists_keys($refhash_omniscient, ('level2', $featureType, $name_geneB)) ){
        if (exists $tarhash_omniscient->{'level2'}{$featureType}{$name_geneB}){

          foreach my $feature_level2_B (@{$tarhash_omniscient->{'level2'}{$featureType}{$name_geneB}}){
            my $name_feature_level2_B=lc($feature_level2_B->_tag_value('ID'));

            my $ref_list_to_checkA=undef;
            if(exists($refhash_omniscient->{'level3'}{'exon'}{$name_feature_level2_A})){
              $ref_list_to_checkA=$refhash_omniscient->{'level3'}{'exon'}{$name_feature_level2_A};
            }
            elsif(exists($refhash_omniscient->{'level3'}{'cds'}{$name_feature_level2_A})){
              $ref_list_to_checkA=$refhash_omniscient->{'level3'}{'cds'}{$name_feature_level2_A}
            }

            my $ref_list_to_checkB=undef;
            if(exists($tarhash_omniscient->{'level3'}{'exon'}{$name_feature_level2_B})){
              $ref_list_to_checkB=$tarhash_omniscient->{'level3'}{'exon'}{$name_feature_level2_B};
            }
            elsif(exists($tarhash_omniscient->{'level3'}{'cds'}{$name_feature_level2_B})){
              $ref_list_to_checkB=$tarhash_omniscient->{'level3'}{'cds'}{$name_feature_level2_B};
            }

            if($ref_list_to_checkA and $ref_list_to_checkB){ # we found at least exon or cds for both
              if(defined(featuresList_overlap($ref_list_to_checkA, $ref_list_to_checkB))){
                $overlap=featuresList_overlap($ref_list_to_checkA, $ref_list_to_checkB);
              }
            }
          }
        }
      }
    }
  }
  return $overlap;
}

sub pushIfNotExit{

  my ($list, $feature)=@_;

  my @listEx=@$list;

  my @tmp = $feature->get_tag_values('ID');
  my $feature_id = shift @tmp;

  foreach my $f (@listEx){

    my @tmp = $f->get_tag_values('ID');
    my $f_id = shift @tmp;

    if ($feature_id eq $f_id){
      return $list;
    }
  }

  push (@listEx, $feature);
  return \@listEx;
}

sub removeElementInList {
  my ($hash_sortbyseq, $feature, $ContigName) = @_;

  foreach my $tagB (keys %{$hash_sortbyseq}){ # allow to work only on gene on the contig
    if(exists ($hash_sortbyseq->{$tagB}{$ContigName})){ # allow to work only on gene on the contig
    #  print "size before ".$#{$hash_sortbyseq->{$tagB}{$ContigName}}."\n";
      @{$hash_sortbyseq->{$tagB}{$ContigName}}= grep { $_ ne $feature } @{$hash_sortbyseq->{$tagB}{$ContigName}}; #remove element of the list
    #  print "size after ".$#{$hash_sortbyseq->{$tagB}{$ContigName}}."\n";
    my @List=($feature);
    }
    else{print "problem when deleting...!\n";}
  }
}

sub sort_by_seq{
  my ($hash_omniscient, $featureType)=@_;

  my %hash_sortBySeq;

  foreach my $tag_level1 (keys %{$hash_omniscient->{'level1'}}){
    foreach my $level1_id (keys %{$hash_omniscient->{'level1'}{$tag_level1}}){
      my $level1_feature = $hash_omniscient->{'level1'}{$tag_level1}{$level1_id};

  #    foreach my $tag_level2 (keys %{$hash_omniscient->{'level2'}}){
        if (exists_keys($hash_omniscient, ('level2', $featureType, $level1_id)) ){ # check if they have mRNA avoiding autovivifcation
          my $firstFeature=$hash_omniscient->{'level2'}{$featureType}{$level1_id}[0];

          if($firstFeature->has_tag('ID')){
            my @mrna_values = $firstFeature->get_tag_values('ID');
           my $mrna_id = shift @mrna_values;
          }
          else{print "tag missing\n";exit;}

          my $position=$level1_feature->seq_id."".$level1_feature->strand;
          push (@{$hash_sortBySeq{$tag_level1}{$position}}, $level1_feature);
        }
    #  }
    }
  }
  return \%hash_sortBySeq;
}

sub testOverlap {
  my ($startA, $endA, $startB, $endB) = @_;

  my $overlap = undef;

  if($startA <= $endB and $endA >= $startB){
      $overlap="yes";
  }

  return $overlap;
}

sub manage_one2one{
  my ($o_one2one,$refhash_omniscient,$reftar_omniscient, $overlapT, $hash_geneName)=@_;

  my $nameChanged=0;
  my $newName=0;
  my $overlapOK=0;

  foreach my $key (keys %{$o_one2one}){

    my $geneA_feature = $o_one2one->{$key}[0];
    my $geneB_feature = $o_one2one->{$key}[1];

    #print "\n\nCheck the gene ".$geneA_feature->gff_string."\nagainst        ".$geneB_feature->gff_string." \n\n";
    my $best_overlap=0;

    my $name_geneA=lc($geneA_feature->_tag_value('ID'));
    my $name_geneB=lc($geneB_feature->_tag_value('ID'));

    foreach my $tag_level2_A (keys %{$refhash_omniscient->{'level2'}}){

      if (exists $refhash_omniscient->{'level2'}{$tag_level2_A}{$name_geneA}){

        foreach my $feature_level2_A (@{$refhash_omniscient->{'level2'}{$tag_level2_A}{$name_geneA}}){

          my $name_feature_level2_A=lc($feature_level2_A->_tag_value('ID'));

          foreach my $tag_level2_B (keys %{$tarhash_omniscient->{'level2'}}){
            if (exists $tarhash_omniscient->{'level2'}{$tag_level2_B}{$name_geneB}){
              foreach my $feature_level2_B (@{$tarhash_omniscient->{'level2'}{$tag_level2_B}{$name_geneB}}){
                my $name_feature_level2_B=lc($feature_level2_B->_tag_value('ID'));

                my $ref_list_to_checkA=undef;
                if(exists($refhash_omniscient->{'level3'}{'exon'}{$name_feature_level2_A})){
                  $ref_list_to_checkA=$refhash_omniscient->{'level3'}{'exon'}{$name_feature_level2_A};
                }
                elsif(exists($refhash_omniscient->{'level3'}{'cds'}{$name_feature_level2_A})){
                  $ref_list_to_checkA=$refhash_omniscient->{'level3'}{'cds'}{$name_feature_level2_A}
                }

                my $ref_list_to_checkB=undef;
                if(exists($tarhash_omniscient->{'level3'}{'exon'}{$name_feature_level2_B})){
                  $ref_list_to_checkB=$tarhash_omniscient->{'level3'}{'exon'}{$name_feature_level2_B};
                }
                elsif(exists($tarhash_omniscient->{'level3'}{'cds'}{$name_feature_level2_B})){
                  $ref_list_to_checkB=$tarhash_omniscient->{'level3'}{'cds'}{$name_feature_level2_B};
                }

                my $overlap_percent=0;
                if(! $inv){
                  if($ref_list_to_checkA and $ref_list_to_checkB){
                    $overlap_percent = featuresList_overlap($ref_list_to_checkA, $ref_list_to_checkB);
                  }
                }
                else{
                  if($ref_list_to_checkA and $ref_list_to_checkB){
                    $overlap_percent = featuresList_overlap($ref_list_to_checkB, $ref_list_to_checkA);
                  }
                }

                if($overlap_percent > $best_overlap){
                  $best_overlap=$overlap_percent;
                }
              }
            }
          }
        }
      }
    }
    #######
    #NOW CHECK IF WE ARE UNDER THE THRESHOLD TO CONSIDER IT
    #print "My best overlap is: $best_overlap\n";
    if($best_overlap >= $overlapT){
      $overlapOK++;
      my ($nameChanged_u, $newName_u) = liftGeneName($geneA_feature, $geneB_feature, $refhash_omniscient, $hash_geneName);
      $nameChanged=$nameChanged+$nameChanged_u;
      $newName=$newName+$newName_u;
    }
  }
return $overlapOK, $nameChanged, $newName;
}

sub liftGeneName{
  my ($geneA_feature, $geneB_feature, $refhash_omniscient, $hash_geneName)=@_;
    my $nameChanged=0;
    my $newName=0;

    my $id_geneA = lc($geneA_feature->_tag_value('ID'));
    my $id_geneB = lc($geneB_feature->_tag_value('ID'));

    my @list_features=($geneA_feature, $geneB_feature);

    foreach my $tag_level1 (keys %{$refhash_omniscient->{'level1'}}){
      if (exists $refhash_omniscient->{'level1'}{$tag_level1}{$id_geneA}){
        my $featureA = $refhash_omniscient->{'level1'}{$tag_level1}{$id_geneA};

        #######################
        # MANAGE NAME (ONLY IF DIFERENT)
        my $same_gene = check_feature_same_names($hash_geneName, \@list_features, 'AB');
        if((! $same_gene) or ($same_gene eq "none")){
          my $tag=undef;
          if($geneB_feature->has_tag('Name')){
            $tag='Name';
          }
          elsif($geneB_feature->has_tag('gene_name')){
            $tag='gene_name';
          }
          elsif($geneB_feature->has_tag('gene_Name')){
            $tag='gene_Name';
          }
          else{print "No tag Name or gene_name for geneB feature !!".$geneB_feature->gff_string."\n";}

          if($tag){
            my $name_featureA="Name Absent";
            if ($geneA_feature->has_tag('Name')){
              $name_featureA = uc($geneA_feature->_tag_value('Name'));
            }
            my $name_featureB=uc($geneB_feature->_tag_value($tag));

            #####
            # We take add label if necessary if not deactivated
            if (! $_dblr){
              if(exists ($hash_geneNameLab{$name_featureB}) ){ # Name already
                my $cpt_label=1;
                my $nameB=$name_featureB."_".$cpt_label;
                while($hash_geneNameLab{$nameB}){
                  $cpt_label++;
                  $nameB=$name_featureB."_".$cpt_label;
                }
              }
            }

            ####################
            # Display info old vs new name
                       #print "\nP1311_101 gene: ".$geneA_feature->gff_string."\n";            #print "K12 gene: ".$geneB_feature->gff_string."\n";
            if($output){
              print $outReport "old name = $name_featureA <=> New name = $name_featureB\n";
            }else{print "old name = $name_featureA <=> New name = $name_featureB\n";}

            create_or_replace_tag($geneA_feature, 'Name', $name_featureB);

            # Track modified gene names
            if($name_featureA eq "Name Absent"){
              $newName=1
            }else{$nameChanged=1;}

            ################
            ## take care of "description" attribute
            my $new_description=get_attribute_value($refhash_omniscient, $geneB_feature, $id_geneB, 'description');
            # Now we try to transfert description on genefeature A
            if($geneA_feature->has_tag('description') and !$new_description){
              $geneA_feature->remove_tag('description');
            }
            elsif($new_description){
              create_or_replace_tag($geneA_feature, 'description', $new_description);
            }

            ################
            ## take care of "product" attribute
            my $new_product=get_attribute_value($refhash_omniscient, $geneB_feature, $id_geneB, 'product');
            # Now we try to transfert product on genefeature A
            if($geneA_feature->has_tag('product') and !$new_product){
              $geneA_feature->remove_tag('product');
            }
            elsif($new_product){
              create_or_replace_tag($geneA_feature, 'product', $new_product);
            }

            #########
            #change now info of all mRNA of geneA
            foreach my $tag_level2 (keys %{$refhash_omniscient->{'level2'}}){
              if (exists $refhash_omniscient->{'level2'}{$tag_level2}{$id_geneA}){
                foreach my $feature (@{$refhash_omniscient->{'level2'}{$tag_level2}{$id_geneA}}){
                  create_or_replace_tag($feature, 'Name', $name_featureB);

                  # Now we try to transfert description on feature level2
                  if($feature->has_tag('description') and !$new_description){
                    $feature->remove_tag('description');
                  }
                  elsif($new_description){
                    create_or_replace_tag($feature, 'description', $new_description);
                  }

                  # Now we try to transfert product on feature level2
                  if($feature->has_tag('product') and !$new_product){
                    $feature->remove_tag('product');
                  }
                  elsif($new_product){
                    create_or_replace_tag($feature, 'product', $new_product);
                  }
                }
              }
            }
          }
        }
        ###########################
        # ADD ORTHOLOGY INFORMATION at gene level
        my $name_geneB=$geneB_feature->_tag_value('ID');
        $name_geneB =~ s/$nbis_suffix_d//;
        create_or_replace_tag($geneA_feature, 'orthology', $name_geneB);
      }
    }

return $nameChanged, $newName;
}

sub get_attribute_value{

  my ($refhash_omniscient, $gene_feature, $id_gene, $attribute)=@_;

  my $value=undef;

  #first we try to get if from genefeature
  if($gene_feature->has_tag($attribute)){ #get it at gene level
    $value=$gene_feature->_tag_value($attribute);
  }
  else{ #Not found at gene level, we try to find it at level2
    foreach my $tag_level2 (keys %{$refhash_omniscient->{'level2'}}){
      if (exists $refhash_omniscient->{'level2'}{$tag_level2}{$id_gene}){
        foreach my $feature (@{$refhash_omniscient->{'level2'}{$tag_level2}{$id_gene}}){
          if($feature->has_tag($attribute)){
            $value=$gene_feature->_tag_value($attribute);
            return $value;
          }
        }
      }
    }
  }
  return $value;
}

#check name of a list of features
sub check_feature_same_names{

  my ($hash_geneName, $ListFeatures, $type)= @_;
  my $sameName="yes";
  my $nameBefore;
  my $gene_cpt=0;

  foreach my $gene_feature (@{$ListFeatures}){
    my $tag=undef;
    if($gene_feature->has_tag('Name')){
      $tag='Name';
    }
    elsif($gene_feature->has_tag('gene_name')){
      $tag='gene_name';
    }
    else{
    #  print "Warning: No name found for ".$gene_feature->gff_string."\n";
      $sameName="none";last;
    }

    if($tag){
      my @tmp=$gene_feature->get_tag_values($tag);
      my $name_gene_feature=lc(shift @tmp);
      $gene_cpt++;

      ####
      #
      my $typeCheck;
      if ($type ne 'AB'){
        $typeCheck = $type;
      }
      else{
        if($gene_cpt == 1){$typeCheck = 'As';} else {$typeCheck = 'Bs';}
      }
      ####
      #  MANAGE IF Name on target comes from nbis functional annotation. We have to remove _1 _2 etc (Not DEFAULT behavior)
      if($_ablt and  $typeCheck eq "Bs") # We remove what has been added by NBIS during gene name annotation
      {
        $name_gene_feature =~ s/$nbis_suffix_p//;
      }
      ####
      #  MANAGE IF Name on reference comes from nbis functional annotation. We have to remove _1 _2 etc (DEFAULT behavior)
      if(! $_dblr and  $typeCheck eq "As") # We remove what has been added by NBIS during gene name annotation
      {
        $name_gene_feature =~ s/$nbis_suffix_p//;
      }

      if ($gene_cpt >= 2){
        if ($nameBefore ne $name_gene_feature){ # the two names are different !
          $sameName=undef;
        }
      }
      $nameBefore=$name_gene_feature;
    }
  }
  return $sameName;
}
=head1 NAME

infer_function_from_synteny.pl

From 2 gffs file from the same assembly, the tool will lift functional information from a target file (--tar file) on the reference file (--ref file).
So the tool the detect the genes from target that overlap, according to the threshold choosen, the genes of the reference file.
Thus, one2one,one2many (fusion), many2one (split) and many2many cases are detected.
Function for one2one cases are liftover, and other cases are reported in corresponding output file.

=head1 SYNOPSIS

    ./infer_function_from_synteny.pl --ref=infile --tar=infile [Options]
    ./infer_function_from_synteny.pl --help

=head1 OPTIONS

=over 8

=item B<--ref>, B<--reffile> or B<-f>

Input GFF3 file correponding to the reference file where function will be written on.

=item B<--tar>, B<--mapped> or B<--m>

Input GFF3 file corresponding to the target file containing function that will be lift over on the reference.

=item B<--inv>, B<--inverse> or B<--oppposite>

By default the ref file genes are tested against the tar file genes for the overlaping percentage. To do the opposite just set this parameter (Do not correspond to a shift between reference and target files).

=item B<--transfert>, B<--lift> or B<-t>

Lift the function of genes from target file to gene of reference file is they overlap on2one with length percentage over "value" parameter.
Particular cases: - one2many => if the many genes from tar have the same name we lift it to the ref according to the overlap percentage value threshold.
                  - many2one => if the many genes from ref have the same name we try to lift the ref name according to the overlap percentage value threshold.

=item B<--cmo>

cno = checkMultiOverlap => When this option is activated, multi-overlaping genes (one2many, many2one (cases not takken in account by the --lift option), and many2many) will be checked.
The most overlaping couple of gene will be taken to try to lift the name according to the overlap percentage value threshold.
If a gene has the same overlapping value within several couple we do not lift-over the function because we cannot define wich one is the real 'ortholog' one.

=item B<--value>, B<--threshold> or B<--overlap>

Define the percentage of overlaping to use to consider genes as ortholog. Usefull only if "lift" option activated.
/!\\ Dont forget your kraken data set has already been selected by an overlapping filer. So if the previous kraken scipt you choose to keep gene mapped over 20 percent;
Here choose a value of 20 will mean you accept to consider the gene as ortholog even if 20% of 20% is mapped (~4%).
Exon features are considered. If there is no exon, cds will be used instead.

=item B<--dblr>

Deactivate NBIS Label Reference. By default to compare names from two files we remove a potential label _(0-9)* at the end of geme names in the reference file because they can have been added by NBIS during the functional annotation process.
If we lift a new name to the reference annotation, we first check that the name is already existing elsewhere in the annotation. If it exists, we also add the labbel according to the number of gene with that name.
If you don't want to take in account the labbel, you can deactivate the behaviour by calling that otpion (--dbl).

=item B<--ablt>

Activate NBIS Label Target. By default we don't look for label in target file. Most of time has not been annotated by NBIS. But in case where the target file is also annotated by NBIS,
it possible to take in account the possible label at the end of gene names by activating the option.

=item B<--feature>

The script checks the overlap using the cds or exon of level2 features.  By default we are considering only mRNA as level2 feature. If you want to change the level2 feature type you can set up that option and defining the new feature to consider. --feature tRNA.
An option "all" has to be fully implemented.

=item  B<--out>, B<--output> or B<-o>

Output directory where diffrent output files will be written.

=item B<--help> or B<-h>

Getting help.
Display the full information.

=back

=cut
