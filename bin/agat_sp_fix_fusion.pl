#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Clone 'clone';
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use List::MoreUtils qw(uniq);
use Bio::DB::Fasta;
use Bio::SeqIO;
use Bio::Tools::CodonTable;
use AGAT::AGAT;


my $start_run = time();
my $startP=time;
my $SIZE_OPT=21;
my $PREFIX_CPT_EXON=1;
my $PREFIX_CPT_MRNA=1;

my $header = get_agat_header();
my $config;
my $outfile = undef;
my $gff = undef;
my $file_fasta=undef;
my $opt_codonTableID=1;
my $stranded=undef;
my $threshold=undef;
my $verbose=undef;
my $opt_help= 0;

my @copyARGV=@ARGV;
if ( !GetOptions(
    'c|config=s'               => \$config,
    "h|help"           => \$opt_help,
    "gff=s"            => \$gff,
    "fasta|fa=s"       => \$file_fasta,
    "stranded|s"       => \$stranded,
    "table|codon|ct=i" => \$opt_codonTableID,
    "verbose|v"        => \$verbose,
    "threshold|t=i"    => \$threshold,
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

if ( ! (defined($gff)) or !(defined($file_fasta)) ){
    pod2usage( {
           -message => "$header\nAt least 2 parameter is mandatory:\nInput reference gff file (--gff) and Input fasta file (--fasta)\n\n",
           -verbose => 0,
           -exitval => 1 } );
}

# --- Manage config ---
$config = get_agat_config({config_file_in => $config});

######################
# Manage output file #
my $gffout_file;
my $gffout2_file;
my $gffout3_file;
my $logout_file;
if ($outfile) {
  my ($filename,$path,$ext) = fileparse($outfile,qr/\.[^.]*/);

  $gffout_file  = $path.$filename."-intact.gff";
  $gffout2_file = $path.$filename."-only_modified.gff";
  $gffout3_file = $path.$filename."-all.gff";
  $logout_file = $path.$filename."-report.txt";
}

my $gffout  = prepare_gffout($config, $gffout_file);
my $gffout2 = prepare_gffout($config, $gffout2_file);
my $gffout3 = prepare_gffout($config, $gffout3_file);
my $logout = prepare_fileout($logout_file);

$opt_codonTableID = get_proper_codon_table($opt_codonTableID);

if(!$threshold){
  $threshold=100;
}
print "Minimum protein length taken in account = $threshold AA\n";

if($stranded){
  $stranded=1;
  print "You say that annotation has been done using stranded RNA. So, most probable fusion will be between close gene in same direction. We will focuse on that !\n";
}
else{ print "You didn't use the option stranded. We will look for fusion in all strand (+ and -)!\n";}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>     MAIN     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

######################
### Parse GFF input #
my ($hash_omniscient, $hash_mRNAGeneLink) = slurp_gff3_file_JD({ input => $gff,
                                                                 config => $config
                                                              });
print ("GFF3 file parsed\n");

####################
# index the genome #
my $db = Bio::DB::Fasta->new($file_fasta);
print ("Fasta file parsed\n");

####################

#counters
my $geneCounter=0;
my $mRNACounter_fixed=0;

my %omniscient_modified_gene; 
initialize_omni_from(\%omniscient_modified_gene, $hash_omniscient);
my @intact_gene_list;

# create the hash temp
my $tmpOmniscient={};
my @mRNAlistToTakeCareR;
my $mRNAlistToTakeCare=\@mRNAlistToTakeCareR;

# manage progression bar variables
my $TotalFeatureL1 = nb_feature_level1($hash_omniscient);
my $featureChecked = 0;
local $| = 1; # Or use IO::Handle; STDOUT->autoflush; Use to print progression bar

foreach my $primary_tag_key_level1 ( keys %{$hash_omniscient->{'level1'}}){ # primary_tag_key_level1 = gene or repeat etc...
  foreach my $gene_id ( keys %{$hash_omniscient->{'level1'}{$primary_tag_key_level1}}){

    #Display progression
    $featureChecked++;
    if ((10 - (time - $startP)) < 0) {
        my $done = ($featureChecked*100)/$TotalFeatureL1;
        $done = sprintf ('%.0f', $done);
        if($verbose) { print "Progress : $done %"; }
        else{ print "\rProgress : $done %"; }
        $startP= time;
    }

    my $gene_feature=$hash_omniscient->{'level1'}{$primary_tag_key_level1}{$gene_id};
    my $oneMRNAmodified=undef;
    my $mrna_pseudo=0;
    my @list_mrna_pseudo;
    my $one_level2_modified; # check if one of the level2 feature will be modified

    # COPY gene and subfeatures in tmpOmniscient.
    $tmpOmniscient = {}; # empty the hash
    @$mRNAlistToTakeCare = (); # empty the list
    my @tmpListID=($gene_id);
    fill_omniscient_from_other_omniscient_level1_id(\@tmpListID,$hash_omniscient,$tmpOmniscient);

    foreach my $primary_tag_key_level2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
      if (exists_keys($hash_omniscient, ('level2', $primary_tag_key_level2, $gene_id)) ){ # check if they have mRNA avoiding autovivifcation
        foreach my $level2_feature ( @{$hash_omniscient->{'level2'}{$primary_tag_key_level2}{$gene_id}}) {

          $PREFIX_CPT_MRNA=1;

          # get multiple info
          my $id_level2 = lc($level2_feature->_tag_value('ID'));
          push (@$mRNAlistToTakeCare,$id_level2);

          ##############################
          #If UTR3 #
          my $oneRoundAgain="yes";
          my $nbNewUTR3gene=0;
          if ( exists_keys($hash_omniscient, ('level3', 'three_prime_utr', $id_level2)) ){

            while($oneRoundAgain){
              if($verbose) {print "\nNew round three_prime_utr\n";}
              my ($breakRound, $nbNewUTRgene, $mRNAlistToTakeCare) = take_care_utr('three_prime_utr', $tmpOmniscient, $mRNAlistToTakeCare, $stranded, $gffout);
              $oneRoundAgain = $breakRound;
              $nbNewUTR3gene += $nbNewUTRgene;
            }
          }
          ##############################
          #If UTR5 #
          $oneRoundAgain="yes";
          my $nbNewUTR5gene=0;
          if ( exists_keys($hash_omniscient, ('level3', 'five_prime_utr', $id_level2)) ){

            while($oneRoundAgain){
                if($verbose) { print "\nNew round five_prime_utr\n";}
                my ($breakRound, $nbNewUTRgene, $mRNAlistToTakeCare) = take_care_utr('five_prime_utr', $tmpOmniscient, $mRNAlistToTakeCare, $stranded, $gffout);
                $oneRoundAgain = $breakRound;
                $nbNewUTR5gene += $nbNewUTRgene;
              }
          }
          ##########################
          #If UTR not well defined #
          if ( exists_keys ($hash_omniscient, ('level3', 'utr', $id_level2) ) ){
            print "Sorry but we need to know which utr it is ... 5 or 3 ?\n";exit;
          }

          #############
          # CHECK AFTER ALL UTR ANALIZED
          my $totalNewUTRgene=$nbNewUTR3gene+$nbNewUTR5gene;
          if($totalNewUTRgene > 0){
            $oneMRNAmodified="yes";
            $mRNACounter_fixed++; # Count only mRNA modified
          }
          @$mRNAlistToTakeCare = (); # empty the list
        } # End foreach mRNA
      }
      if($oneMRNAmodified){
        $geneCounter++;
        $oneMRNAmodified=undef;
        #save remodelate gene name
        merge_omniscients(\%omniscient_modified_gene, $tmpOmniscient);
      }
      else{push(@intact_gene_list, $gene_id);}
    }
  }
}
# end progreesion bar
if($verbose) { print "Progress : 100 %\n"; }
else{print "\rProgress : 100 %\n"; }

###
# Fix frame
fil_cds_frame(\%omniscient_modified_gene, $db, $opt_codonTableID);
fil_cds_frame($hash_omniscient, $db, $opt_codonTableID);

#####################################
# Manage modified gene to be sure they not overlap already existing gene. If yes => we give the same gene ID and remove one.
print "Managing spurious labelling at gene level\n";
# 1) create a hash omniscient intact
my $hash_omniscient_intact={}; initialize_omni_from($hash_omniscient_intact, $hash_omniscient);
fill_omniscient_from_other_omniscient_level1_id(\@intact_gene_list, $hash_omniscient, $hash_omniscient_intact);
delete $hash_omniscient->{$_} for (keys %{$hash_omniscient});

# 2) print the intact one
print "print intact...\n";
print_omniscient( {omniscient => $hash_omniscient_intact, output => $gffout} );

# 3) Sort by seq_id - review all newly created gene
my $hash_sortBySeq = gather_and_sort_l1_by_seq_id_and_strand($hash_omniscient_intact);
my $overlap=0;

foreach my $tag_l1 ( keys %{$omniscient_modified_gene{'level1'}} ){ # primary_tag_key_level1 = gene or repeat etc...
  foreach my $id_l1 ( keys %{$omniscient_modified_gene{'level1'}{$tag_l1}} ) {
    #find_overlap_between_geneFeature_and_sortBySeqId can remove the level1 ID so we have to check it is still present
    if( exists_keys ( \%omniscient_modified_gene, ('level1', $tag_l1, $id_l1 ) ) ){
      my $geneFeature = $omniscient_modified_gene{'level1'}{$tag_l1}{$id_l1};
      if (find_overlap_between_geneFeature_and_sortBySeqId($geneFeature, \%omniscient_modified_gene, $hash_omniscient_intact, $hash_sortBySeq) ){
        $overlap++
      }
    }
  }
}

# 4) special case where two newly created gene from to different gene are overlapping
# Be careful If you by testing 2 identical omniscient, the method could remove element haven't yet been loop over. So check the gene exists before to analyse it !
$hash_sortBySeq = gather_and_sort_l1_by_seq_id_and_strand(\%omniscient_modified_gene);

foreach my $tag_l1 ( keys %{$omniscient_modified_gene{'level1'}} ){ # primary_tag_key_level1 = gene or repeat etc...
  foreach my $id_l1 ( keys %{$omniscient_modified_gene{'level1'}{$tag_l1}} ) {
    if( exists_keys( \%omniscient_modified_gene, ('level1', $tag_l1, $id_l1 ) ) ) {
      my $geneFeature = $omniscient_modified_gene{'level1'}{$tag_l1}{$id_l1};
      if (find_overlap_between_geneFeature_and_sortBySeqId($geneFeature, \%omniscient_modified_gene, \%omniscient_modified_gene, $hash_sortBySeq) ){
        $overlap++
      }
    }
  }
}

# 5) Print modified genes
print "print modified...\n";
if (exists_undef_value(\%omniscient_modified_gene)){print"there is an undef value";exit;}

print_omniscient( {omniscient => \%omniscient_modified_gene, output => $gffout2} );

# 6) Print all together
merge_omniscients_fuse_l1duplicates($hash_omniscient_intact, \%omniscient_modified_gene);
print "print all together...\n";
print_omniscient( {omniscient => $hash_omniscient_intact, output => $gffout3} );

if ($overlap and $verbose){print "We found $overlap case gene overlapping at CDS level wihout the same ID, we fixed them.\n";}
# End manage overlaping name
#####################################

#END
my $string_to_print="usage: $0 @copyARGV\n";
$string_to_print .="Results:\n";
$string_to_print .="$geneCounter genes affected and $mRNACounter_fixed mRNA.\n";
my $end_run = time();
my $run_time = $end_run - $start_run;
$string_to_print .= "Job done in $run_time seconds\n";

if($outfile){
  print $logout $string_to_print
}
print $string_to_print;
print "Bye Bye.\n";
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

sub merge_omniscients_fuse_l1duplicates {
  my ($hash_omniscient1, $hash_omniscient2)=@_;

  # == LEVEL 1 == #

  foreach my $tag_l1 (keys %{$hash_omniscient2->{'level1'}}){
    foreach my $id_l1_2 (keys %{$hash_omniscient2->{'level1'}{$tag_l1}}){

      if ( ! exists_keys ( $hash_omniscient1, ('level1', $tag_l1, $id_l1_2 ) ) ){
        my $feature = $hash_omniscient2->{'level1'}{$tag_l1}{lc($id_l1_2)};
        $hash_omniscient1->{'level1'}{$tag_l1}{$id_l1_2} = $hash_omniscient2->{'level1'}{$tag_l1}{$id_l1_2}; # save feature level1
      }

      # == LEVEL 2 == #

      foreach my $tag_l2 (keys %{$hash_omniscient2->{'level2'}}){
        if (exists_keys ($hash_omniscient2, ('level2', $tag_l2, $id_l1_2) ) ){
          foreach my $feature_l2 ( @{$hash_omniscient2->{'level2'}{$tag_l2}{$id_l1_2}}) {

            my $id_l2 = lc($feature_l2->_tag_value('ID'));

            # == LEVEL 3 == #

            foreach my $tag_l3 (keys %{$hash_omniscient2->{'level3'}}){

              if (exists_keys ($hash_omniscient2, ('level3', $tag_l3, $id_l2) ) ){
                $hash_omniscient1->{'level3'}{$tag_l3}{$id_l2} = delete $hash_omniscient2->{'level3'}{$tag_l3}{$id_l2}; # save @l3
              }
            }
          }

          if ( ! exists_keys ( $hash_omniscient1, ('level2', $tag_l2, $id_l1_2 ) ) ){
            $hash_omniscient1->{'level2'}{$tag_l2}{$id_l1_2} = delete $hash_omniscient2->{'level2'}{$tag_l2}{$id_l1_2}; # save @l2
          }
          else{
            push @{$hash_omniscient1->{'level2'}{$tag_l2}{$id_l1_2}}, @{$hash_omniscient2->{'level2'}{$tag_l2}{$id_l1_2}};
          }
        }
      }
    }
  }
  return $hash_omniscient1;
}

# @Purpose: The hash of reference will be the Hash target (HashT). The kept name will come from the hash of reference.
# When an overlap is found, the ID/parent are fixed and we return 1 as a success !
# @input: 4 => object(gene feature), hash(omniscient), hash(omniscient), hash(sortBySeq)
# @output: 1 => undef || integer(1)
sub find_overlap_between_geneFeature_and_sortBySeqId {
  my ($geneFeature, $hash_source, $hashT, $hashT_sortBySeq )=@_;

  my $tag = lc($geneFeature->primary_tag);
  my $seqid = $geneFeature->seq_id;
  my $strand = $geneFeature->strand;
  my $gene_idS = $geneFeature->_tag_value('ID');

  #find overlap
  my $total_overlap=0;
  my $nb_feat_overlap=0;
  my @ListOverlapingGene=();

  foreach my $gene_featureT ( @{$hashT_sortBySeq->{"$seqid$strand"}{$tag}}){

    my $gene_idT = $gene_featureT->_tag_value('ID');

    if($gene_idT eq $gene_idS){ next;} # avoid to compare same feature if we are checking same omniscient

    my ($start1,$end1) = get_most_right_left_cds_positions($hashT,$gene_idT); # look at CDS because we want only overlapinng CDS
    if( !$start1 or !$end1){next;} # No cds so no cds locations
    my ($start2,$end2) = get_most_right_left_cds_positions($hash_source,$gene_idS); # look at CDS becaus we want only ioverlapinng CDS
    if( !$start2 or !$end2){next;} # No cds so no cds locations

    if( ($start2 <= $end1) and ($end2 >= $start1) ){ #feature overlap considering extrem start and extrem stop. It's just to optimise the next step. Avoid to do the next step every time. So at the end, that test (current one) could be removed
                                                     # Even if true, they do not necessarly overlap on the spreded features
          #now check at each CDS feature independently
          if (_two_features_overlap_two_hashes($hash_source,$gene_idS, $hashT, $gene_idT)){
            #print "These two features overlap without same id ! :\n".$geneFeature->gff_string."\n".$gene_featureT->gff_string."\n";
            $nb_feat_overlap++;

            push(@ListOverlapingGene, $gene_featureT);
          }
      }
  }

   # Now manage name if some feature overlap
  if( $nb_feat_overlap > 0){
    my $reference_feature = shift(@ListOverlapingGene);
      push(@ListOverlapingGene, $geneFeature);
      #print "$nb_feat_overlap overlapping feature found ! We will treat them now:\n";
      #print "We decided to keep that one: ".$reference_feature->gff_string."\n";

      my $gene_id_ref  = $reference_feature->_tag_value('ID');

      #change level2 parent for feature of level2 that have a feature of level1 in $ListToRemove list
      foreach my $featureToRemove (@ListOverlapingGene){

        my $gene_id_to_remove  = lc($featureToRemove->_tag_value('ID'));

        #######
        #which hash the feature come from ?
        my $currentHash=undef;
        foreach my $tag_l1 (keys %{$hash_source->{'level1'}} ){ # primary_tag_key_level1 = gene or repeat etc...
        if($hash_source->{'level1'}{$tag_l1}{$gene_id_to_remove} ){
          $currentHash = $hash_source;
        }
      }
      if(! $currentHash){$currentHash = $hashT;}
      # ok now hash is choosen
      ################

        foreach my $tag_level2 (keys %{$currentHash->{'level2'}}){

            if (exists_keys($currentHash, ('level2',$tag_level2,$gene_id_to_remove)) ){ # check if they have cds avoiding autovivification.

            my @list_tmp_features = @{$currentHash->{'level2'}{$tag_level2}{$gene_id_to_remove}}; # As we will remove element of the list we cannot loop over it directly, we have to save the list in a temporary list;
            foreach my $level2_feature (@list_tmp_features){ #replace Parent of each feature level2 by the new level1 reference
              # Change parent feature
              create_or_replace_tag($level2_feature,'Parent',$gene_id_ref);

              #add it in other list
              push (@{$currentHash->{'level2'}{$tag_level2}{lc($gene_id_ref)}},$level2_feature);

              #remove mRNA from list <= not mandatory
              my $mrna_id_to_remove = $level2_feature->_tag_value('ID');
              my @tag_list=('all');
              my @id_list=($gene_id_to_remove);my @id_list2=($mrna_id_to_remove);

              remove_element_from_omniscient(\@id_list, \@id_list2, $currentHash, 'level2', 'false', \@tag_list);

            }
          }
        }

        foreach my $tag_level1 (keys %{$currentHash->{'level1'}}){ # remove the old feature level1 now
          my $new_l1_feature = clone($reference_feature);
          delete $currentHash->{'level1'}{$tag_level1}{$gene_id_to_remove}; # delete level1
          $currentHash->{'level1'}{$tag_level1}{lc($gene_id_ref)} = $new_l1_feature;
        }

      } #END FEATURE TO HANDLE
      ###
      # check end and start of the new feature
      my $gene_id=lc($reference_feature->_tag_value('ID'));
      check_level1_positions( { omniscient => $hashT, feature => $reference_feature } );
      return 1;
  }
  else{return undef;}
}

# @Purpose: Check if two genes have at least one mRNA isoform which overlap at cds level.
# @input: 4 => hash(omniscient), string(gene identifier), hash(omniscient), string(gene identifier)
# @output: 1 => undef || string(yes)
sub _two_features_overlap_two_hashes{
  my  ($hash1, $gene_id1, $hash2, $gene_id2)=@_;
  my $resu=undef;

  #check full CDS for each mRNA
  foreach my $mrna_feature (@{$hash1->{'level2'}{'mrna'}{lc($gene_id1)}}){
    foreach my $mrna_feature2 (@{$hash2->{'level2'}{'mrna'}{lc($gene_id2)}}){

      my $mrna_id1 = $mrna_feature->_tag_value('ID');
      my $mrna_id2 = $mrna_feature2->_tag_value('ID');

      #check all cds pieces
      foreach my $cds_feature1 (@{$hash1->{'level3'}{'cds'}{lc($mrna_id1)}}){
        foreach my $cds_feature2 (@{$hash2->{'level3'}{'cds'}{lc($mrna_id2)}}){

          if(($cds_feature2->start <= $cds_feature1->end) and ($cds_feature2->end >= $cds_feature1->start )){ # they overlap
            $resu="yes";last;
          }
        }
        if($resu){last;}
      }
      if($resu){last;}
    }
    if($resu){last;}
  }
  return $resu;
}

sub take_care_utr{

  my ($utr_tag, $tmpOmniscient, $mRNAlistToTakeCare, $stranded, $gffout)=@_;

  my $oneRoundAgain=undef;
  my $nbNewUTRgene=0;

  foreach my $primary_tag_key_level1 (keys %{$tmpOmniscient->{'level1'}}){ # primary_tag_key_level1 = gene or repeat etc...
    foreach my $gene_id (sort keys %{$tmpOmniscient->{'level1'}{$primary_tag_key_level1}}){

    my $gene_feature=$tmpOmniscient->{'level1'}{$primary_tag_key_level1}{$gene_id};
    my $gene_id = lc($gene_feature->_tag_value('ID'));
    #print "\ntake care utr GeneID = $gene_id\n";

      foreach my $primary_tag_key_level2 (sort keys %{$tmpOmniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
        if (exists_keys($tmpOmniscient, ('level2', $primary_tag_key_level2, $gene_id)) ){ # check if they have mRNA avoiding autovivifcation
          foreach my $level2_feature ( @{$tmpOmniscient->{'level2'}{$primary_tag_key_level2}{$gene_id}}) {

            my $id_level2=lc($level2_feature->_tag_value('ID'));
            foreach my $mRNAtoTakeCare (@{$mRNAlistToTakeCare}){

              if($mRNAtoTakeCare eq $id_level2){ # ok is among the list of those to analyze
                if($verbose) { print "id_level2 -- $id_level2 ***** to take_care -- $mRNAtoTakeCare  \n";}
                if(exists_keys ($tmpOmniscient, ('level3', $utr_tag) ) and exists_keys ($tmpOmniscient, ('level3', $utr_tag,$id_level2) ) ){

                  ##################################################
                  # extract the concatenated exon and cds sequence #
                  my $oppDir=undef;
                  my $original_strand=$level2_feature->strand;

                  ###############
                  # Manage UTRS #
                  my @utr_feature_list = sort {$a->start <=> $b->start} @{$tmpOmniscient->{'level3'}{$utr_tag}{$id_level2}}; # be sure that list is sorted
                  my ($utrExtremStart, $utr_seq, $utrExtremEnd) = concatenate_feature_list(\@utr_feature_list);

                  #If UTR shorter than the minimum DNA size expected, we skip it => WE SAVE TIME
                  if(length($utr_seq) < ($threshold*3) ){
                    next;
                  }

                  #create the utr object
                  my $utr_obj = Bio::Seq->new(-seq => $utr_seq, -alphabet => 'dna' );

                  #Reverse complement according to strand
                  if ($original_strand == -1 or $original_strand eq "-"){
                    $utr_obj = $utr_obj->revcom();
                  }

                  # get the revcomp
                  my $opposite_utr_obj = $utr_obj->revcom();


                  my $longest_ORF_prot_obj;
                  my $orf_utr_region;
                  #################################
                  # Get the longest ORF positive ## record ORF = start, end (half-open), length, and frame
                  my ($longest_ORF_prot_obj_p, $orf_utr_region_p) = translate_JD($utr_obj,
                                                                              -orf => 'longest',
                                                                              -codontable_id => $opt_codonTableID);
                  ########################################
                  # Get the longest ORF opposite strand ## record ORF = start, end (half-open), length, and frame
                  my $length_longest_ORF_prot_obj_n=0;
                  my $longest_ORF_prot_obj_n;
                  my $orf_utr_region_n;

                  if(! $stranded){
                    ($longest_ORF_prot_obj_n, $orf_utr_region_n) = translate_JD($opposite_utr_obj,
                                                                                -orf => 'longest',
                                                                                -codontable_id => $opt_codonTableID);
                    $length_longest_ORF_prot_obj_n = $longest_ORF_prot_obj_n->length();
                  }

                  #################
                  # Choose the best
                  if($longest_ORF_prot_obj_p->length() >= $length_longest_ORF_prot_obj_n){
                    $longest_ORF_prot_obj= $longest_ORF_prot_obj_p;
                    $orf_utr_region= $orf_utr_region_p;
                  }
                  else{
                    $longest_ORF_prot_obj= $longest_ORF_prot_obj_n;
                    $orf_utr_region= $orf_utr_region_n;
                    $oppDir=1;
                    #my @cds_feature_list = sort {$a->start <=> $b->start} @{$tmpOmniscient->{'level3'}{'cds'}{$id_level2}}; # be sure that list is sorted
                    #($cdsExtremStart, $cds_dna_seq, $cdsExtremEnd) = concatenate_feature_list($cds_feature_list); # we have to change these value because it was not predicted as same direction as mRNA
                  }


                  ########################
                  # prediction is longer than threshold#
                  if($longest_ORF_prot_obj->length() > $threshold){

                    if($verbose) {print "Longer AA in utr = ".$longest_ORF_prot_obj->length()."\n".$longest_ORF_prot_obj->seq."\n";}

                    my @exons_features = sort {$a->start <=> $b->start} @{$tmpOmniscient->{'level3'}{'exon'}{$id_level2}};# be sure that list is sorted
                    my ($exonExtremStart, $mrna_seq, $exonExtremEnd) = concatenate_feature_list(\@exons_features);

                    my @cds_feature_list = sort {$a->start <=> $b->start} @{$tmpOmniscient->{'level3'}{'cds'}{$id_level2}}; # be sure that list is sorted
                    my ($cdsExtremStart, $cds_dna_seq, $cdsExtremEnd) = concatenate_feature_list(\@cds_feature_list);

                    # set real start and stop to orf
                    my $realORFstart;
                    my $realORFend;
                    #print "mRNA length: ".length($mrna_seq)."  UTR length: ".length($utr_seq)."\n";
                    #print "start in UTR piece ".$orf_utr_region->[0]." end ".$orf_utr_region->[1]."\n";

                    ####################################
                    # Recreate position of start in mRNA positive strand
                    my $startUTRinMRNA = length($mrna_seq) - length($utr_seq);
                    if ($utr_tag eq 'three_prime_utr' ){
                      if($original_strand == 1 or $original_strand eq "+" ){
                        if(! $oppDir){
                          $orf_utr_region->[0] += $startUTRinMRNA;
                        }
                        else{ #opposite direction
                          $orf_utr_region->[0] = length($mrna_seq) - $orf_utr_region->[1];
                        }
                      }
                      else{ #minus strand
                          if(! $oppDir){
                            $orf_utr_region->[0] = length($utr_seq) - $orf_utr_region->[1]; #flip position
                          }
                      }
                    }
                    elsif ($utr_tag eq 'five_prime_utr'){
                      if($original_strand == 1 or $original_strand eq "+"){
                        if($oppDir){
                          $orf_utr_region->[0]=length($utr_seq) - $orf_utr_region->[1];
                        }
                      }
                      else{ #minus strand
                        if(! $oppDir){
                          $orf_utr_region->[0] = (length($utr_seq) - $orf_utr_region->[1])+$startUTRinMRNA;
                        }
                        else{ #opposite direction
                           $orf_utr_region->[0] += $startUTRinMRNA;
                        }
                      }
                    }

                    #calcul the real start end stop of utr in genome
                    ($realORFstart, $realORFend) = calcul_real_orf_end_and_start($orf_utr_region, \@exons_features);

                    #save the real start and stop
                    $orf_utr_region->[0]=$realORFstart;
                    $orf_utr_region->[1]=$realORFend;

                    # Now manage splitting the old gene to obtain two genes
                    $mRNAlistToTakeCare = split_gene_model($tmpOmniscient, $gene_feature, $level2_feature, \@exons_features, \@cds_feature_list, $cdsExtremStart, $cdsExtremEnd, $realORFstart, $realORFend, $oppDir, $mRNAlistToTakeCare, $gffout);

                    $oneRoundAgain="yes";
                    $nbNewUTRgene++;
                  } # We predict something in UTR
                  else{ if($verbose) { print "Nothing predicted over threshold :". $longest_ORF_prot_obj->length()." ! Next\n";} }
                } # End there is UTR
                else{ if($verbose) {print "There is no UTR ! Next\n";} }
              }
              #else{print "Not among the list mRNAtoTakeCare. Next \n";}
            }
          }
        }
      }
    }
  }
  return $oneRoundAgain, $nbNewUTRgene, $mRNAlistToTakeCare;
}

############
# P.S: when a gene is newly created, it has a new name even if it overlap at CDS level an another gene that is not part of the current temporary omniscient studied. So, an extra step at the end will catch and fix those kind of cases.
sub split_gene_model{

   my ($tmpOmniscient, $gene_feature, $level2_feature, $exons_features, $cds_feature_list, $cdsExtremStart, $cdsExtremEnd, $realORFstart, $realORFend, $oppDir, $mRNAlistToTakeCare, $gffout)=@_;

      my $gene_id = $gene_feature->_tag_value('ID');
      my $id_level2 = lc($level2_feature->_tag_value('ID'));
      my $newcontainerUsed=0;

                  ######################
                  # Recreate exon list #
                  my $bolean_original_is_first;
                  my $first_end;
                  my $second_start;
                  #if new prediction after on the sequence
                  if($realORFstart >= $cdsExtremEnd){
                    $bolean_original_is_first="true";
                    $first_end=$cdsExtremEnd;
                    $second_start=$realORFstart;
                  }else{ # ($realORFend < $cdsExtremStart)
                    $bolean_original_is_first="false";
                    $first_end=$realORFend;
                    $second_start=$cdsExtremStart;
                  }
                  my ($newOrignal_exon_list, $newPred_exon_list) = create_two_exon_lists($tmpOmniscient, $exons_features, $first_end, $second_start, $bolean_original_is_first, $oppDir);

        ####################################
        # Remodelate ancient gene
        ####################################
                  if($verbose) { print "Remodelate ancient gene\n"; }
                  #############################################################
                  #  Remove all level3 feature execept cds
                  my @tag_list=('cds');
                  my @l2_id_list=($id_level2);
                  remove_tuple_from_omniscient(\@l2_id_list, $tmpOmniscient, 'level3', 'false', \@tag_list);


                  #############
                  # Recreate original exon
                  @{$tmpOmniscient->{'level3'}{'exon'}{$id_level2}}=@$newOrignal_exon_list;

                  #########
                  #RE-SHAPE last/first exon if less than 3 nucleotides (1  or 2 must be romved) when the CDS finish 1 or 2 nuclotide before... because cannot be defined as UTR
                  shape_exon_extremity($newOrignal_exon_list,$cds_feature_list);

                  ########
                  # calcul utr
                  if($verbose) { print "Remodelate ancient gene ($gene_id)".$gene_feature->start." ".$gene_feature->end."\n";}

                  my ($original_utr5_list, $variable_not_needed, $original_utr3_list) = modelate_utr_and_cds_features_from_exon_features_and_cds_start_stop($newOrignal_exon_list, $cdsExtremStart, $cdsExtremEnd);
                  @{$tmpOmniscient->{'level3'}{'five_prime_utr'}{$id_level2}}=@$original_utr5_list;
                  @{$tmpOmniscient->{'level3'}{'three_prime_utr'}{$id_level2}}=@$original_utr3_list;


                  ####
                  # Check existance
                  my ($new_gene, $new_mrna, $overlaping_gene_ft, $overlaping_mrna_ft) = must_be_a_new_gene_new_mrna($tmpOmniscient, $cds_feature_list, $newOrignal_exon_list);


                  if ($new_mrna){
                    #########
                    #RE-SHAPE mrna extremities
                    check_mrna_positions({ l2_feature => $level2_feature,
                                           exon_list => $newOrignal_exon_list});

                  }
                  else{
                    if($verbose) { print "*** remove IT *** because exon and CDS IDENTIK ! $id_level2 \n"; }
                    my @l2_feature_list=($level2_feature);
                    remove_omniscient_elements_from_level2_feature_list($tmpOmniscient, \@l2_feature_list);
                  }

                  #########
                  #RE-SHAPE gene extremities
                  check_level1_positions( { omniscient => $tmpOmniscient, feature => $gene_feature } );

        ###################################
        # Remodelate New Prediction
        ###################################
                  if($verbose) { print "Remodelate New Prediction\n"; }
                  # If newPred_exon_list list is empty we skipt the new gene modeling part
                  #if(!@$newPred_exon_list){
                  #  next;
                  #}
                  ###############################################
                  # modelate level3 features for new prediction #
                  my ($new_pred_utr5_list, $new_pred_cds_list, $new_pred_utr3_list) = modelate_utr_and_cds_features_from_exon_features_and_cds_start_stop($newPred_exon_list, $realORFstart, $realORFend);

                  ####################################
                  #RE-SHAPE last/first exon if less than 3 nucleotides (1  or 2 must be removed) when the CDS finish 1 or 2 nuclotide before... because cannot be defined as UTR
                  if(shape_exon_extremity($newPred_exon_list, $new_pred_cds_list)){
                    #we reshaped the exon, it means that the UTR are not correct anymore, we have to recalculate them
                    ($new_pred_utr5_list, $new_pred_cds_list, $new_pred_utr3_list) = modelate_utr_and_cds_features_from_exon_features_and_cds_start_stop($newPred_exon_list, $realORFstart, $realORFend);
                  }

                  my @level1_list;
                  my @level2_list;
                  my @level3_list;
                  my $transcript_id = $newPred_exon_list->[0]->_tag_value('Parent');
                  #############################################
                  # Modelate gene features for new prediction #

                  # $containerUsed exist when we already use the gene container. So in the case where we have only one mRNA, the split will give 2 mRNA. One is linked to the original gene container (done before)
                  # The second must be linked to a new gene container. So, even if must_be_a_new_gene method say no, we must create it because the original one has been already used.
                  ($new_gene, $new_mrna, $overlaping_gene_ft, $overlaping_mrna_ft) = must_be_a_new_gene_new_mrna($tmpOmniscient, $new_pred_cds_list, $newPred_exon_list);
                  if ( $new_gene ){
                    $newcontainerUsed++;
                    $gene_id = take_care_gene_id($gene_id, $tmpOmniscient);

                    my $new_gene_feature = Bio::SeqFeature::Generic->new(-seq_id => $newPred_exon_list->[0]->seq_id, -source_tag => $newPred_exon_list->[0]->source_tag, -primary_tag => 'gene' , -start => $newPred_exon_list->[0]->start,  -end => $newPred_exon_list->[$#{$newPred_exon_list}]->end, -frame => $newPred_exon_list->[0]->frame, -strand => $newPred_exon_list->[0]->strand , -tag => { 'ID' => $gene_id }) ;
                    @level1_list=($new_gene_feature);
                    #print "create_a_new_gene for ".$transcript_id." !!!! - ".$new_gene_feature->gff_string."\n";

                  }
                  else{ #the new mRNA still overlap an isoform. So we keep the link with the original gene

                    # change gene ID
                    $gene_id = $overlaping_gene_ft->_tag_value('ID');
                    #print "We use $gene_id\n";
                    check_level1_positions( { omniscient => $tmpOmniscient, feature => $overlaping_gene_ft } );
                    @level1_list=($overlaping_gene_ft);
                  }

                  #############################################
                  # Modelate mRNA features for new prediction #
                  if ( $new_mrna ){
                    my $new_mRNA_feature = Bio::SeqFeature::Generic->new(-seq_id => $newPred_exon_list->[0]->seq_id, -source_tag => $newPred_exon_list->[0]->source_tag, -primary_tag => $level2_feature->primary_tag() , -start => $newPred_exon_list->[0]->start,  -end => $newPred_exon_list->[$#{$newPred_exon_list}]->end, -frame => $newPred_exon_list->[0]->frame, -strand => $newPred_exon_list->[0]->strand , -tag => { 'ID' => $transcript_id , 'Parent' => $gene_id }) ;
                    push (@$mRNAlistToTakeCare, lc($transcript_id));
                    @level2_list=($new_mRNA_feature);

                    @level3_list=(@$newPred_exon_list, @$new_pred_cds_list, @$new_pred_utr5_list, @$new_pred_utr3_list);

                    #Save the gene (not necesserely new) and mRNA feature (necesseraly new)
                    append_omniscient($tmpOmniscient, \@level1_list, \@level2_list, \@level3_list);

                    #Now we have the new transcript we can test the gene end and start
                    if ( $new_gene ){
                      check_level1_positions( { omniscient => $tmpOmniscient, feature => $level1_list[0] } );
                    }
                    else{
                      check_level1_positions( { omniscient => $tmpOmniscient, feature => $overlaping_gene_ft } );
                    }
                  }
                  else{
                    if($verbose){print "*** Not creating mRNA *** because exon and CDS IDENTIK ! \n";}
                  }

  return $mRNAlistToTakeCare;
}


#create an Uniq gene ID
sub take_care_gene_id{

      my ($gene_id, $tmpOmniscient) = @_;

      #clean geneid if necessary
      $gene_id =~ /^(new[0-9]+_)?(.*)$/;
      my $clean_id=$2;

      #count current gene number - should be one if first analysis
      my $primary_tag_key_general;
      my $numberOfNewGene=1;

      foreach my $primary_tag_key_level1 (keys %{$tmpOmniscient->{'level1'}}){ # primary_tag_key_level1 = gene or repeat etc...
        foreach my $gene_id_from_hash (keys %{$tmpOmniscient->{'level1'}{$primary_tag_key_level1}}){
          if($gene_id_from_hash =~ /(new[1-9]+_)/){
            $numberOfNewGene++;
          }
        }
        #primary tag key containg gene name has been found. No need to see the others.
        $primary_tag_key_general=$primary_tag_key_level1;
        last;
      }

      # From tmpOmniscient: Between ( new1_geneA, new2_geneA, new3_geneA). It happens that new2_geneA has been deleted. In that case we try to create new3_geneA but as already exists we try new2_geneA (--)
      # If new2_geneA also already exist, it means that it exist in hash_omniscient. so we will try decrementing $numberGeneIDToCheck until 1; Then we will try incrementing $numberGeneIDToCheck over new3_geneA  (in other term we try new4_geneA )
      my $testok=undef;
      my $nbToadd=-1;
      my $numberGeneIDToCheck=$numberOfNewGene;
      my $new_id;
      while (! $testok){
        my $newGenePrefix="new".$numberGeneIDToCheck."_";
        $new_id="$newGenePrefix$clean_id";

        if((! defined ($tmpOmniscient->{'level1'}{$primary_tag_key_general}{lc($new_id)})) and (! defined ($hash_omniscient->{'level1'}{$primary_tag_key_general}{lc($new_id)}))){
            $testok=1;
        }
        else{
          if($numberGeneIDToCheck == 1){
            $nbToadd=1;$numberGeneIDToCheck=$numberOfNewGene;
          }
          $numberGeneIDToCheck += $nbToadd;}
      }
      #print "old_gene_id --- $gene_id ***** new_gene_id --- $new_id\n";

  return $new_id;
}

#create an Uniq mRNA ID
sub take_care_mrna_id {

      my ($tmpOmniscient, $mRNA_id) = @_;

      #clean geneid if necessary
      $mRNA_id =~ /^(new[0-9]+_)?(.*)$/;
      my $clean_id=$2;

      #count current gene number - should be one if first analysis
      my %id_to_avoid;

      foreach my $primary_tag_key_level1 (keys %{$tmpOmniscient->{'level1'}}){ # primary_tag_key_level1 = gene or repeat etc...
        foreach my $gene_id_from_hash (keys %{$tmpOmniscient->{'level1'}{$primary_tag_key_level1}}){

          foreach my $primary_tag_key_level2 (keys %{$tmpOmniscient->{'level2'}}){ # primary_tag_key_level1 = gene or repeat etc...
            if( exists_keys($tmpOmniscient, ('level2', $primary_tag_key_level2, $gene_id_from_hash)) ){
              foreach my $featureL2 (@{$tmpOmniscient->{'level2'}{$primary_tag_key_level2}{$gene_id_from_hash}}){
                my $mrna_id_from_hash=$featureL2->_tag_value('ID');
                if($mrna_id_from_hash =~ /(new[1-9]+_)/){
                  $id_to_avoid{lc($mrna_id_from_hash)}++;
                  $PREFIX_CPT_MRNA++;
                }
              }
            }
          }
        }
      }

      my $testok=undef;
      my $nbToadd=-1;
      my $numberMRNA_IDToCheck=$PREFIX_CPT_MRNA;
      my $new_id;
      while (! $testok){
        my $newPrefix="new".$numberMRNA_IDToCheck."_";
        $new_id="$newPrefix$clean_id";

        if( (! defined ($hash_mRNAGeneLink->{lc($new_id)})) and (! defined ($id_to_avoid{lc($new_id)})) ) {
            $testok=1;
        }
        else{
          if($numberMRNA_IDToCheck == 1){
            $nbToadd=1;$numberMRNA_IDToCheck=$PREFIX_CPT_MRNA;
          }
          $numberMRNA_IDToCheck += $nbToadd;}
      }
      #print "old_mrna_id --- $mRNA_id ***** new_mrna_id --- $new_id\n";

  return $new_id;
}

#As based on a Uniq mRNA ID, this will create a Uniq ID;
#PREFIX_CPT_EXON allows to kepp track of name already given during a exon list spliting
sub take_care_level3_id {

      my ($tmpOmniscient, $feature) = @_;

      #clean geneid if necessary
      my $level3_id = $feature->_tag_value('ID');
      $level3_id =~ /^(new[0-9]+_)?(.*)$/;
      my $clean_id=$2;
      my $newPrefix="new".$PREFIX_CPT_EXON."_";
      my $new_id="$newPrefix$clean_id";

      my $primary_tag = lc($feature->primary_tag);

      while(ID_exists_at_level3($tmpOmniscient, $new_id, $primary_tag )){
        $PREFIX_CPT_EXON++;
        $new_id =~ /^(new[0-9]+_)?(.*)$/;
        my $clean_id=$2;
        my $newPrefix="new".$PREFIX_CPT_EXON."_";
        $new_id="$newPrefix$clean_id";
      }
  return $new_id;
}

#return undef if the ID is not existing in tmpOmniscient
sub ID_exists_at_level3{

  my ($tmpOmniscient, $ID, $primary_tag ) = @_;

      foreach my $level2_ID (keys %{$tmpOmniscient->{'level3'}{$primary_tag}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...

        foreach my $feature_l3 ( @{$tmpOmniscient->{'level3'}{$primary_tag}{$level2_ID}}) {
          my $existingID = $feature_l3->_tag_value('ID');
          if ($existingID eq $ID){
             return 1;
          }
        }
      }
  return undef;
}

# Yes if mRNA doesnt overlap an other existing isoform
# mRNA "true" true mean no overlap at CDS level
sub must_be_a_new_gene_new_mrna{
  my ($omniscient, $new_pred_cds_list, $newPred_exon_list)=@_;

  my $overlaping_mrna_ft=undef;
  my $overlaping_gene_ft=undef;
  my $Need_new_gene="true";
  my $Need_new_mRNA="true";
  my $strand=$new_pred_cds_list->[0]->strand;

  foreach my $primary_tag_key_level1 (keys %{$omniscient->{'level1'}}){ # primary_tag_key_level1 = gene or repeat etc...
    foreach my $gene_id_from_hash (keys %{$omniscient->{'level1'}{$primary_tag_key_level1}}){
      my $gene_feature= $omniscient->{'level1'}{$primary_tag_key_level1}{$gene_id_from_hash};

      if($strand eq $gene_feature->strand){
        foreach my $primary_tag_key_level2 (keys %{$omniscient->{'level2'}}){ # primary_tag_key_level1 = gene or repeat etc...
          if( exists_keys($omniscient, ('level2', $primary_tag_key_level2, $gene_id_from_hash)) ){

            foreach my $featureL2 (@{$omniscient->{'level2'}{$primary_tag_key_level2}{$gene_id_from_hash}}){

            # get level2 id
            my $featureL2_id = lc($featureL2->_tag_value('ID'));
            my $featureL2_original_id=lc($newPred_exon_list->[0]->_tag_value('Parent'));

              if($featureL2_id ne $featureL2_original_id){

                #Now check if overlap
                my @cds_feature_list = @{$omniscient->{'level3'}{'cds'}{$featureL2_id}};
                my @exon_feature_list = @{$omniscient->{'level3'}{'exon'}{$featureL2_id}};

                my $overlap_cds = featuresList_overlap(\@cds_feature_list, $new_pred_cds_list);
                if(defined ($overlap_cds)){ #If CDS overlap
                  $Need_new_gene=undef;
                  $overlaping_gene_ft=$gene_feature;
                  #print "CDS Overlap entre $featureL2_id and $featureL2_original_id !\n";
                  if(featuresList_identik(\@cds_feature_list, $new_pred_cds_list)){
                    #print "cds identik !\n";
                    if(featuresList_identik(\@exon_feature_list, $newPred_exon_list)){
                      if($verbose) { print "RNA identik BETWEEN $featureL2_id and $featureL2_original_id \n"; }
                      $Need_new_mRNA=undef;
                      $overlaping_mrna_ft=$featureL2_id;
                      last;
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

  return $Need_new_gene, $Need_new_mRNA, $overlaping_gene_ft, $overlaping_mrna_ft;
}

#remove small remaining pieces in the of UTR in the exon shape
sub shape_exon_extremity{
  #exon_features is a sorted list
  #cds_features is a sorted list
  my $modified=undef;
  my ($exon_features,$cds_features)=@_;

   #test between first exon and first cds
   if( (abs($cds_features->[0]->start - $exon_features->[0]->start) < 3) and (abs($cds_features->[0]->start - $exon_features->[0]->start) > 0) ){ #We have to shape the exon start. We don't want a non multiple of 3 inferior to 3
      $exon_features->[0]->start($cds_features->[0]->start);
      $modified=1;
   }
   #test between last exon and last cds
   if(abs($exon_features->[$#{ $exon_features }]->end - $cds_features->[$#{ $cds_features }]->end ) < 3){  #We have to shape the exon end
      $exon_features->[$#{ $exon_features }]->end($cds_features->[$#{ $cds_features }]->end);
      $modified=1;
   }
   return $modified;
}

sub calcul_real_orf_end_and_start{
  #exons_features is sorted
  my ($orf_cds_region, $exons_features)=@_;

  my $realORFstart;
  my $realORFend;

  my $orf_start=$orf_cds_region->[0]; # get start to begin
  my $orf_length=$orf_cds_region->[2]; # get lentgh to map

  my $first="yes";
  my $total_exon_length=0;
  my $total_exon_length_previous_round=0;
  my $mapped_length=0;
  my $mapped_length_total=0;
  my $the_rest_to_map=0;

  foreach my $exon_feature (@$exons_features){
    # Allows to follow the path on mRNA
    my $exon_length=($exon_feature->end - $exon_feature->start)+1;
    $total_exon_length_previous_round=$total_exon_length;
    $total_exon_length += $exon_length;
    # Allows to follow the path on the CDS
    $mapped_length_total += $mapped_length;
    $the_rest_to_map=$orf_length-$mapped_length_total;
    # exon overlap CDS
    if($total_exon_length >= $orf_start){ #they begin to overlap

      if($first eq "yes"){
        #  $realORFstart=$exon_feature->start+($orf_start - 1);
        $realORFstart=$exon_feature->start+($orf_start - $total_exon_length_previous_round );
        my $end_part_of_exon=$exon_feature->start- $realORFstart + 1;
        if($end_part_of_exon >= $orf_length){           #exon      ============================================
           $realORFend=$realORFstart+$orf_length-1;       #cds              =========================
           last;
         }
        $mapped_length=$exon_feature->end - $realORFstart + 1;
        $first="no";
      }
      else{
        $mapped_length=$exon_feature->end - $exon_feature->start + 1;
      }
    }
    #exon are over the end of cds => we finish at this round
    if($total_exon_length >= ($orf_start+$orf_length) ){        #exon      ============================================
      if($realORFstart > $exon_feature->start){                 #cds       =========================
        $realORFend=$realORFstart+$the_rest_to_map - 1 ;
      last;
      }else{
        $realORFend=$exon_feature->start + $the_rest_to_map - 1 ;
      last;
      }
    }
  }
return $realORFstart, $realORFend;
}

sub change_strand{
  my ($feature)=@_;

  if($feature->strand eq "-" or $feature->strand eq "-1"){
    $feature->strand('+');
  }else{$feature->strand('-');}
}

# The exons containing the original cds keep their parent names. The exon containing the new cds will have a new parent name.
sub create_two_exon_lists {
  # orignalFirst == true if original gene is first on the prediction
  my ($tmpOmniscient, $exons_features, $firstEnd, $secondStart, $orignalFirst, $oppDir)=@_;
  my @list_exon_originalPred;
  my @list_exon_newPred;
  #print "firstEnd $firstEnd, secondStart $secondStart, $orignalFirst, $oppDir\n";
  $PREFIX_CPT_EXON=1;

  my $value = $exons_features->[0]->_tag_value('Parent');
  my $NewParentName = take_care_mrna_id($tmpOmniscient, $value);

  foreach my $exon_feature (@$exons_features){ #for each exon
    if(two_positions_on_feature($exon_feature,$firstEnd,$secondStart)){  # We have to split the exon_feature
      my $duplicated_exon_feature=clone($exon_feature);#create a copy of the feature

      $exon_feature->end($secondStart-1);
      $duplicated_exon_feature->start($secondStart);

      if($orignalFirst eq "true"){

        push( @list_exon_originalPred, $exon_feature);

        my $value = take_care_level3_id($tmpOmniscient,$duplicated_exon_feature);
        create_or_replace_tag($duplicated_exon_feature,'ID', $value);
        create_or_replace_tag($duplicated_exon_feature,'Parent', $NewParentName);
        if($oppDir){
          change_strand($duplicated_exon_feature);
        }
        push( @list_exon_newPred, $duplicated_exon_feature);
        next;
      }else{ #original pred after
        $duplicated_exon_feature->start($secondStart-1);
        push( @list_exon_originalPred, $duplicated_exon_feature);

        my $value = take_care_level3_id($tmpOmniscient, $exon_feature);
        create_or_replace_tag($exon_feature,'ID', $value);

        create_or_replace_tag($exon_feature,'Parent', $NewParentName);
        if($oppDir){
          change_strand($exon_feature);
        }
        push( @list_exon_newPred, $exon_feature);
        next;
      }
    }
    if(! (($exon_feature->end <=  $secondStart) and ($exon_feature->start >=  $firstEnd))){ # avoid exon between CDSs
      if ($exon_feature->end <=  $secondStart) {
        if ($orignalFirst eq "true"){
          push( @list_exon_originalPred, $exon_feature);
        }else{
          my $duplicated_exon_feature=clone($exon_feature);#create a copy of the feature
          my $value = take_care_level3_id($tmpOmniscient, $duplicated_exon_feature);
          create_or_replace_tag($duplicated_exon_feature,'ID', $value);
          create_or_replace_tag($duplicated_exon_feature,'Parent', $NewParentName);
          if($oppDir){
            change_strand($duplicated_exon_feature);
          }
          push( @list_exon_newPred, $duplicated_exon_feature);
        }
      }
      if ($exon_feature->start >=  $firstEnd) {
        if($orignalFirst eq "true"){
          my $duplicated_exon_feature=clone($exon_feature);#create a copy of the feature
          my $value = take_care_level3_id($tmpOmniscient, $duplicated_exon_feature);
          create_or_replace_tag($duplicated_exon_feature,'ID', $value);
          create_or_replace_tag($duplicated_exon_feature,'Parent', $NewParentName);
           if($oppDir){
            change_strand($duplicated_exon_feature);
          }
          push( @list_exon_newPred, $duplicated_exon_feature);
        }
        else{
          push( @list_exon_originalPred, $exon_feature);
        }
      }
    }
    if(($exon_feature->end <=  $secondStart) and ($exon_feature->start >=  $firstEnd)){ # Exon between CDSs
      if ($orignalFirst eq "true"){
        push( @list_exon_originalPred, $exon_feature);
      }else{
        my $duplicated_exon_feature=clone($exon_feature);#create a copy of the feature
        my $value = take_care_level3_id($tmpOmniscient, $duplicated_exon_feature);
        create_or_replace_tag($duplicated_exon_feature,'ID', $value);
        create_or_replace_tag($duplicated_exon_feature,'Parent', $NewParentName);
        if($oppDir){
          change_strand($duplicated_exon_feature);
        }
        push( @list_exon_newPred, $duplicated_exon_feature);
      }
    }
  }
  my @list_exon_originalPred_sorted = sort {$a->start <=> $b->start} @list_exon_originalPred;
  my @list_exon_newPred_sorted = sort {$a->start <=> $b->start} @list_exon_newPred;
  #  print "list1: @list_exon_originalPred_sorted\n";
  #  foreach my $u (@list_exon_originalPred_sorted){
  #    print $u->gff_string."\n";
  #  }
  #  print "list2: @list_exon_newPred_sorted\n";
  #  foreach my $u (@list_exon_newPred_sorted){
  #    print $u->gff_string."\n";
  # }
  return \@list_exon_originalPred_sorted, \@list_exon_newPred_sorted;
}

#Check if feature overlap one position
sub position_on_feature {

  my ($feature,$position)=@_;

  my $isOnSameExon=undef;
  if ( ($position >= $feature->start and $position <= $feature->end)){
    $isOnSameExon="true";
  }
  return $isOnSameExon;
}

#Check if feature overlap two positions (start and stop)
sub two_positions_on_feature {

  my ($feature,$position1,$position2)=@_;

  my $areOnSameExon=undef;
  if ( ($position1 >= $feature->start and $position1 <= $feature->end) and ($position2 >= $feature->start and $position2 <= $feature->end) ){
    $areOnSameExon="true";
  }
  return $areOnSameExon;
}

# We do not use the official translate function from the PrimarySeqI object/lib because we want to keep track to the ORF positions too. So we have modified it consequently.
sub translate_JD {
   my ($self,@args) = @_;
     my ($terminator, $unknown, $frame, $codonTableId, $complete,
     $complete_codons, $throw, $codonTable, $orf, $start_codon, $offset);

   ## new API with named parameters, post 1.5.1
   if ($args[0] && $args[0] =~ /^-[A-Z]+/i) {
         ($terminator, $unknown, $frame, $codonTableId, $complete,
         $complete_codons, $throw,$codonTable, $orf, $start_codon, $offset) =
       $self->_rearrange([qw(TERMINATOR
                                               UNKNOWN
                                               FRAME
                                               CODONTABLE_ID
                                               COMPLETE
                                               COMPLETE_CODONS
                                               THROW
                                               CODONTABLE
                                               ORF
                                               START
                                               OFFSET)], @args);
   ## old API, 1.5.1 and preceding versions
   } else {
     ($terminator, $unknown, $frame, $codonTableId,
      $complete, $throw, $codonTable, $offset) = @args;
   }

    ## Initialize termination codon, unknown codon, codon table id, frame
    $terminator = '*'    unless (defined($terminator) and $terminator ne '');
    $unknown = "X"       unless (defined($unknown) and $unknown ne '');
    $frame = 0           unless (defined($frame) and $frame ne '');
    $codonTableId = 1    unless (defined($codonTableId) and $codonTableId ne '');
    $complete_codons ||= $complete || 0;

    ## Get a CodonTable, error if custom CodonTable is invalid
    if ($codonTable) {
     $self->throw("Need a Bio::Tools::CodonTable object, not ". $codonTable)
      unless $codonTable->isa('Bio::Tools::CodonTable');
    } else {

        # shouldn't this be cached?  Seems wasteful to have a new instance
        # every time...
    $codonTable = Bio::Tools::CodonTable->new( -id => $codonTableId);
   }

    ## Error if alphabet is "protein"
    $self->throw("Can't translate an amino acid sequence.") if
    ($self->alphabet =~ /protein/i);

    ## Error if -start parameter isn't a valid codon
   if ($start_codon) {
     $self->throw("Invalid start codon: $start_codon.") if
      ( $start_codon !~ /^[A-Z]{3}$/i );
   }

   my $seq;

   if ($offset) {
    $self->throw("Offset must be 1, 2, or 3.") if
        ( $offset !~ /^[123]$/ );
    my ($start, $end) = ($offset, $self->length);
    ($seq) = $self->subseq($start, $end);
   } else {
    ($seq) = $self->seq();
   }

         ## ignore frame if an ORF is supposed to be found
   my $orf_region;
   if ( $orf ) {
            ($orf_region) = $self->_find_orfs_nucleotide($seq, $codonTable, $start_codon, $orf eq 'longest' ? 0 : 'first_only' );
            $seq = $self->_orf_sequence( $seq, $orf_region );
   } else {
   ## use frame, error if frame is not 0, 1 or 2
     $self->throw("Valid values for frame are 0, 1, or 2, not $frame.")
      unless ($frame == 0 or $frame == 1 or $frame == 2);
     $seq = substr($seq,$frame);
         }

    ## Translate it
    my $output = $codonTable->translate($seq, $complete_codons);
    # Use user-input terminator/unknown
    $output =~ s/\*/$terminator/g;
    $output =~ s/X/$unknown/g;

    ## Only if we are expecting to translate a complete coding region
    if ($complete) {
     my $id = $self->display_id;
     # remove the terminator character
     if( substr($output,-1,1) eq $terminator ) {
       chop $output;
     } else {
       $throw && $self->throw("Seq [$id]: Not using a valid terminator codon!");
       $self->warn("Seq [$id]: Not using a valid terminator codon!");
     }
     # test if there are terminator characters inside the protein sequence!
     if ($output =~ /\Q$terminator\E/) {
             $id ||= '';
       $throw && $self->throw("Seq [$id]: Terminator codon inside CDS!");
       $self->warn("Seq [$id]: Terminator codon inside CDS!");
     }
     # if the initiator codon is not ATG, the amino acid needs to be changed to M
     if ( substr($output,0,1) ne 'M' ) {
       if ($codonTable->is_start_codon(substr($seq, 0, 3)) ) {
         $output = 'M'. substr($output,1);
       }  elsif ($throw) {
         $self->throw("Seq [$id]: Not using a valid initiator codon!");
       } else {
         $self->warn("Seq [$id]: Not using a valid initiator codon!");
       }
     }
    }

    my $seqclass;
    if ($self->can_call_new()) {
     $seqclass = ref($self);
    } else {
     $seqclass = 'Bio::PrimarySeq';
     $self->_attempt_to_load_Seq();
    }
    my $out = $seqclass->new( '-seq' => $output,
                    '-display_id'  => $self->display_id,
                    '-accession_number' => $self->accession_number,
                    # is there anything wrong with retaining the
                    # description?
                    '-desc' => $self->desc(),
                    '-alphabet' => 'protein',
                    '-verbose' => $self->verbose
            );
    return $out, $orf_region;
}

sub concatenate_feature_list{

  my ($feature_list) = @_;

  my $seq = "";
  my $ExtremStart=1000000000000;
  my $ExtremEnd=0;

  foreach my $feature (@$feature_list) {
    my $start=$feature->start();
    my $end=$feature->end();
    my $seqid=$feature->seq_id();
    $seq .= $db->seq( $seqid, $start, $end );

    if ($start < $ExtremStart){
      $ExtremStart=$start;
    }
    if($end > $ExtremEnd){
              $ExtremEnd=$end;
    }
  }
   return $ExtremStart, $seq, $ExtremEnd;
}

__END__

=head1 NAME

agat_sp_fix_fusion.pl

=head1 DESCRIPTION

The script looks for other ORF in UTRs (UTR3 and UTR5) of each gene model described in the gff file.
Several ouput files will be written if you specify an output.
One will contain the gene not modified (intact), one the gene models fixed.

=head1 SYNOPSIS

    agat_sp_fix_fusion.pl --gff infile.gff --fasta genome.fa [ -o outfile ]
    agat_sp_fix_fusion.pl --help

=head1 OPTIONS

=over 8

=item B<-gff>

Input GTF/GFF file.

=item B<-fa> or B<--fasta>

Input fasta file.

=item B<--ct>, B<--codon> or B<--table>

Codon table to use. [default 1]

=item B<-t> or B<--threshold>

This is the minimum length of new protein predicted that will be taken in account.
By default this value is 100 AA.

=item B<-s> or B<--stranded>

By default we predict protein in UTR3 and UTR5 and in both direction. The fusion assumed can be between gene in same direction and in opposite direction.
If RNAseq data used during the annotation was stranded, only fusion of close genes oriented in same direction are expected. In that case this option should be activated.
When activated, we will try to predict protein in UTR3 and UTR5 only in the same orientation than the gene investigated.

=item B<-v> or B<--verbose>

Output verbose information.

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
