#!/usr/bin/env perl

## if IUPAC:
## We consider a stop only if we are sure it is one
## CDS can contains putative stop codon (but not sure stop one like YAA that can be TAA or CAA).
## We consider a start even if is not sure like AYG that can be ATG or ACG

##TO DO
## Consider longest ORF wihtout checking start (can be incomplete) <= otpion to check start

use strict;
use warnings;
use Carp;
use Clone 'clone';
use File::Basename;
use List::MoreUtils qw(uniq);
use Bio::DB::Fasta;
use Bio::SeqIO;
use AGAT::AGAT;

# avoid case of ambiguous start codon (translated into X) -> we accept if the ORF is SIZE_OPT AA longer.
# Indeed statistically it has more chance to be a real start codon.
my $SIZE_OPT = 21;

my $header = get_agat_header();
my @copyARGV = @ARGV;
my ( $opt, $usage, $config ) = AGAT::AGAT::describe_script_options(
    $header,
    [ 'gff=s',        'Input reference gff file',   { required => 1 } ],
    [ 'fasta|fa|f=s', 'Input reference fasta file', { required => 1 } ],
    [ 'split|s!',     'Split sequences' ],
    [ 'table|codon|ct=i', 'Codon translation table',
        { default => 1, callbacks => { positive => sub { shift > 0 or die 'must be positive' } } } ],
    [ 'model|m=s', 'Model(s) to test',
        { callbacks => { allowed => sub { my $val = shift; $val =~ /^([1-6](,[1-6])*)?$/ or die 'model must be comma-separated list of integers 1-6'; 1; } } } ],
);

my $gff          = $opt->gff;
my $file_fasta   = $opt->fasta;
my $split_opt    = $opt->split;
my $codonTable   = $opt->table;
my $model_to_test = $opt->model;
my $outfile      = $config->{output};
my $verbose      = $config->{verbose};

my $log;
if ( my $log_name = $config->{log_path} ) {
    open( $log, '>', $log_name )
      or die "Can not open $log_name for printing: $!";
    dual_print( $log, $header, 0 );
}

# --- Check codon table
$codonTable = get_proper_codon_table($codonTable, $log, $verbose);

######################
# Manage output file #
my $gffout_file;
my $gffout2_file;
my $gffout3_file;
#my $gffout4_file;
my $report_file;

if ($outfile) {
  my ($path,$ext);
  ($outfile,$path,$ext) = fileparse($outfile,qr/\.[^.]*/);

  $gffout_file  = $path.$outfile."-intact.gff";
  $gffout2_file = $path.$outfile."-only_modified.gff";
  $gffout3_file = $path.$outfile."-all.gff";
  #$gffout4_file = $path.$outfile."-pseudogenes.gff";
  $report_file  = $path.$outfile."-report.txt";
}

my $gffout  = prepare_gffout($config, $gffout_file);
my $gffout2 = prepare_gffout($config, $gffout2_file);
my $gffout3 = prepare_gffout($config, $gffout3_file);
#my $gffout4 = prepare_gffout($config, $gffout4_file);
my $report  = prepare_fileout($report_file);

my %ListModel;
if(!($model_to_test)){
  $ListModel{1}=0;
  $ListModel{2}=0;
  $ListModel{3}=0;
  $ListModel{4}=0;
  $ListModel{5}=0;
  $ListModel{6}=0;
}else{
  my @fields= split(',', $model_to_test);
  foreach my $field (@fields){
      $ListModel{$field}=0;
  }
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>     MAIN     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

######################
### Parse GFF input #
my ($hash_omniscient, $hash_mRNAGeneLink) = slurp_gff3_file_JD({ input => $gff,
                                                                config => $config
                                                              });
dual_print( $log, "GFF3 file parsed\n" );


####################
# index the genome #
my $db = Bio::DB::Fasta->new($file_fasta);
dual_print( $log, "Fasta file parsed\n" );

####################
my $pseudo_threshold=70;
#counters
my $counter_case21=0;
my $geneCounter=0;
my $mRNACounter=0;
my $mRNACounter_fixed=0;
#my $mrna_pseudo_suspected=0;
#my $gene_pseudo_suspected=0;
#my $mrna_pseudo_removed=0;
#my $gene_pseudo_removed=0;

my %omniscient_modified_gene; initialize_omni_from(\%omniscient_modified_gene, $hash_omniscient);
#my %omniscient_pseudogene;
my @modified_gene_list;
my @intact_gene_list;

foreach my $primary_tag_key_level1 (keys %{$hash_omniscient->{'level1'}}){ # primary_tag_key_level1 = gene or repeat etc...
  foreach my $gene_id_tag_key (keys %{$hash_omniscient->{'level1'}{$primary_tag_key_level1}}){
    my $gene_feature=$hash_omniscient->{'level1'}{$primary_tag_key_level1}{$gene_id_tag_key};

    my $one_ORFmodified="no";
    #my $mrna_pseudo=0;
    #my @list_mrna_pseudo;
    my $number_mrna=0;

    foreach my $primary_tag_key_level2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
      if ( exists_keys( $hash_omniscient, ('level2', $primary_tag_key_level2, $gene_id_tag_key) ) ){
        foreach my $level2_feature ( @{$hash_omniscient->{'level2'}{$primary_tag_key_level2}{$gene_id_tag_key}}) {

          my $ORFmodified="no";
          $number_mrna=$#{$hash_omniscient->{'level2'}{$primary_tag_key_level2}{$gene_id_tag_key}}+1;

          # get level2 id
          my $id_level2 = lc($level2_feature->_tag_value('ID'));

          ##############################
          #If it's a mRNA = have CDS. #
          if ( exists ($hash_omniscient->{'level3'}{'cds'}{$id_level2} ) ){

            ##############
            # Manage CDS #
            my @cds_feature_list = sort {$a->start <=> $b->start} @{$hash_omniscient->{'level3'}{'cds'}{$id_level2}}; # be sure that list is sorted
						shrink_cds_offset(\@cds_feature_list);
            my ($cdsExtremStart, $cds_dna_seq, $cdsExtremEnd) = concatenate_feature_list(\@cds_feature_list);
            #create the cds object
            my $cds_obj = Bio::Seq->new(-seq => $cds_dna_seq, -alphabet => 'dna' );
            #Reverse the object depending on strand
            if ($level2_feature->strand == -1 or $level2_feature->strand eq "-"){
              $cds_obj = $cds_obj->revcom();
            }
            #translate cds in protein
            my $original_prot_obj = $cds_obj->translate(-codontable_id => $codonTable) ; #codontable_id by default=0 strict M as start codon. IUPAC => STOP codon even if not sure ...
            my $cds_prot=$original_prot_obj->seq;
            #print $original_prot_obj->seq."\n";
            my $originalProt_size=length($cds_prot);

            ################################################
            # mRNA: extract the concatenated exon sequence #
            my @exons_features = sort {$a->start <=> $b->start} @{$hash_omniscient->{'level3'}{'exon'}{$id_level2}};
            my ($exonExtremStart, $mrna_seq, $exonExtremEnd) = concatenate_feature_list(\@exons_features);
            #create the mrna object
            my $mrna_obj = Bio::Seq->new(-seq => $mrna_seq, -alphabet => 'dna' );

            #Reverse complement according to strand
            if ($level2_feature->strand == -1 or $level2_feature->strand eq "-"){
              $mrna_obj = $mrna_obj->revcom();
            }

            #######################
            # Get the longest ORF ## record ORF = start, end (half-open), length, and frame
            my ($longest_ORF_prot_obj, $orf_cds_region) = translate_JD($mrna_obj,
                                                                        -orf => 'longest',
                                                                        -codontable_id => $codonTable);
      #     print Dumper($orf_cds_region)."\n";
            # set real start and stop to orf
            my $realORFstart;
            my $realORFend;
            # change the start for negative strand
            if ($level2_feature->strand == -1 or $level2_feature->strand eq "-"){
              $orf_cds_region->[0]=(length($mrna_seq) - $orf_cds_region->[1]);
            }
            #calcul the real start end stop of cds in genome
     #       print Dumper($orf_cds_region)."\n".$mrna_obj->seq."\n";
            ($realORFstart, $realORFend) = calcul_real_orf_end_and_start($orf_cds_region, \@exons_features);
    #        print "$id_level2 $realORFstart $realORFend\n";
            #save the real start and stop
            $orf_cds_region->[0]=$realORFstart;
            $orf_cds_region->[1]=$realORFend;

  #############
  # Tests     #
  #############

            ########################
            # prediction is longer #
            dual_print( $log, $id_level2." - size before: ".$originalProt_size." size after: ".$longest_ORF_prot_obj->length()."\n" );
                                            dual_print( $log, "Original: ".$original_prot_obj->seq."\n", 3 );
                                            dual_print( $log, "Prediction: ".$longest_ORF_prot_obj->seq."\n", 3 );
            if($longest_ORF_prot_obj->length() > $originalProt_size){

  #Model1     ###############################################
              # sequence original is part of new prediction #
              if (index($longest_ORF_prot_obj->seq,$cds_prot) != -1){
                if ( exists($ListModel{1}) ){

									if( pass_ambiguous_start($longest_ORF_prot_obj, $originalProt_size) ){

                    $ListModel{1}++;
                    dual_print($log, "Model 1: gene=$gene_id_tag_key mRNA=$id_level2\n");
                    dual_print( $log, "original:$cds_prot\nnew:". $longest_ORF_prot_obj->seq."\n" );
                    modify_gene_model($hash_omniscient, \%omniscient_modified_gene, $gene_feature, $gene_id_tag_key, $level2_feature, $id_level2, \@exons_features, \@cds_feature_list, $cdsExtremStart, $cdsExtremEnd, $realORFstart, $realORFend, 'model1', $gffout);
                    $ORFmodified="yes";

                  }
                }
              }
              #################################################
              # protein original and predicted are different
              else{

                #########################################
  #Model2       # Prediction don't overlap original CDS #
                if( ($realORFend < $cdsExtremStart) or ($realORFstart > $cdsExtremEnd) ){
                  my $model;

                  if( exists($ListModel{2}) ){
                    $ListModel{2}++;
                    dual_print($log, "Model 2: gene=$gene_id_tag_key mRNA=$id_level2\n");
                    $model=1;
                    if($split_opt){
                      split_gene_model(\@intact_gene_list, $hash_omniscient, \%omniscient_modified_gene, $gene_feature, $gene_id_tag_key, $level2_feature, $id_level2, \@exons_features, \@cds_feature_list, $cdsExtremStart, $cdsExtremEnd, $realORFstart, $realORFend, 'model2', $gffout);
                    }
                    else{
                      modify_gene_model($hash_omniscient, \%omniscient_modified_gene, $gene_feature, $gene_id_tag_key, $level2_feature, $id_level2, \@exons_features, \@cds_feature_list, $cdsExtremStart, $cdsExtremEnd, $realORFstart, $realORFend, 'model2', $gffout);
                    }
                    $ORFmodified="yes";
                  }
                } # End they don't overlap

                #########################################
                # Prediction Overlap original CDS #
                else { # They overlap
  #Model3         ###############
                  # original protein and predicted one are different; the predicted one is longest, they overlap each other.
                  if( exists($ListModel{3}) ){
                    $ListModel{3}++;
                    dual_print($log, "Model 3: gene=$gene_id_tag_key mRNA=$id_level2\n");

                    modify_gene_model($hash_omniscient, \%omniscient_modified_gene, $gene_feature, $gene_id_tag_key, $level2_feature, $id_level2, \@exons_features, \@cds_feature_list, $cdsExtremStart, $cdsExtremEnd, $realORFstart, $realORFend, 'model3', $gffout);
                    $ORFmodified="yes";
                  }
                }
              }
            }# End prediction longer

            ###########################
            # The real ORF looks to be shorter than the one originaly described ! Selenocysteine ? pseudogene ? Or just case where prediction begin by L instead of M (correct !) or begin by XXXXXX
            elsif($longest_ORF_prot_obj->length() < $originalProt_size){

              # contains stop codon but not at the last position
              if( (index($original_prot_obj->seq, '*') != -1 ) and (index($original_prot_obj->seq, '*') != length($original_prot_obj->seq)-1) ){
                dual_print( $log, "Original sequence contains premature stop codon.\n" );
                #Model4     ###############
                # /!\ Here we compare the CDS traduction (traduct in IUPAC) against longest CDS in mRNA IUPAC modified to take in account for stops codon only those that are trustable only (TGA, TAR...).
                if( exists($ListModel{4}) ){
                  $ListModel{4}++;
                  dual_print($log, "Model 4: gene=$gene_id_tag_key mRNA=$id_level2\n");
                  dual_print( $log, "Original: ".$original_prot_obj->seq."\n" );
                  dual_print( $log, "longestl: ".$longest_ORF_prot_obj->seq."\n" );
                  #remodelate a shorter gene
                  modify_gene_model($hash_omniscient, \%omniscient_modified_gene, $gene_feature, $gene_id_tag_key, $level2_feature, $id_level2, \@exons_features, \@cds_feature_list, $cdsExtremStart, $cdsExtremEnd, $realORFstart, $realORFend, 'model4', $gffout);
                  $ORFmodified="yes";
             ## Pseudogene THRESHOLD ##
    #              my $threshold_size=(length($original_prot_obj->seq)*$pseudo_threshold)/100; #70% of the original size
    #              if(length($longest_ORF_prot_obj->seq) <  $threshold_size){ # inferior to threshold choosen, we suspect it to be a pseudogene
                    #print Dumper($original_prot_obj);
                    #print Dumper($longest_ORF_prot_obj);
    #                $mrna_pseudo++;
    #                push(@list_mrna_pseudo, $id_level2);
                }
              }
                # no premature stop codons
              else{
                if( exists($ListModel{5}) ){
                  $ListModel{5}++;
                  dual_print( $log, "Model 5: gene=$gene_id_tag_key mRNA=$id_level2\n" );
                  dual_print( $log, "Original: ".$original_prot_obj->seq."\n" );
                  dual_print( $log, "longestl: ".$longest_ORF_prot_obj->seq."\n" );
                  #remodelate a shorter gene
                  modify_gene_model($hash_omniscient, \%omniscient_modified_gene, $gene_feature, $gene_id_tag_key, $level2_feature, $id_level2, \@exons_features, \@cds_feature_list, $cdsExtremStart, $cdsExtremEnd, $realORFstart, $realORFend, 'model4', $gffout);
                  $ORFmodified="yes";
                }
              }
            }
            ###########################
            # The real ORF is same size but check if +1 or +2 bp shit that give same number of AA but give frame shifts
            elsif( (index($original_prot_obj->seq, '*') != -1 ) and (index($original_prot_obj->seq, '*') != length($original_prot_obj->seq)-1) ){
              if( exists($ListModel{6}) ){
                $ListModel{6}++;
                dual_print( $log, "Model 6: gene=$gene_id_tag_key mRNA=$id_level2\n" );
                dual_print( $log, "Original: ".$original_prot_obj->seq."\n" );
                dual_print( $log, "longestl: ".$longest_ORF_prot_obj->seq."\n" );
                modify_gene_model($hash_omniscient, \%omniscient_modified_gene, $gene_feature, $gene_id_tag_key, $level2_feature, $id_level2, \@exons_features, \@cds_feature_list, $cdsExtremStart, $cdsExtremEnd, $realORFstart, $realORFend, 'model3', $gffout);
                $ORFmodified="yes";
              }
            }


          } # End there is a CDS
          if($ORFmodified eq "yes"){
            $one_ORFmodified="yes";
            $mRNACounter_fixed++; # Count only mRNA modified
          }
        } # End foreach mRNA
      }
  #    if($mrna_pseudo > 0){
        # all mRNA are pseudogene, we change the gene status to pseudogenes.
  #      if($mrna_pseudo == $number_mrna){
  #        $mrna_pseudo_suspected += $number_mrna;
  #        $gene_pseudo_suspected++;
  #       $gene_feature->primary_tag('pseudogene');
          #transfert the gene and sub-feature to the omniscient_pseudogene hash
  #        my @level1_list=($gene_id_tag_key);
  #        fill_omniscient_from_other_omniscient_level1_id(\@level1_list, $hash_omniscient, \%omniscient_pseudogene); # If already exists in omniscient_modified_gene, it will be replaced by the modified one
  #      }
        #only some of the isoform are pseudo... we remove them
  #      else{
  #        $mrna_pseudo_removed += $mrna_pseudo;
  #        $gene_pseudo_removed++;
  #        my @tag_list=('all');
  #        my @id_list=($gene_id_tag_key);
  #        remove_element_from_omniscient(\@id_list, \@list_mrna_pseudo, $hash_omniscient, 'level2', 'false', \@tag_list);
  #        remove_tuple_from_omniscient(\@list_mrna_pseudo, $hash_omniscient, 'level3', 'false', \@tag_list);
  #        remove_element_from_omniscient(\@id_list, \@list_mrna_pseudo, \%omniscient_modified_gene, 'level2', 'false', \@tag_list);
  #        remove_tuple_from_omniscient(\@list_mrna_pseudo, \%omniscient_modified_gene, 'level3', 'false', \@tag_list);
  #        print "@list_mrna_pseudo has been removed because are isoform containing stop codon\n";
  #      }
  #    }
    }
    if($one_ORFmodified eq "yes"){
      $geneCounter++;
      $mRNACounter += $number_mrna; #add all the mRNA if at least one modified
      #save remodelate gene name
      push(@modified_gene_list, $gene_id_tag_key);
    }
    else{push(@intact_gene_list, $gene_id_tag_key);}

  }
}

###########
# Fix frame
fil_cds_frame(\%omniscient_modified_gene, $db, $log, $verbose, $codonTable);
#fil_cds_frame(\%omniscient_pseudogene);
fil_cds_frame($hash_omniscient, $db, $log, $verbose, $codonTable);

#Clean omniscient_modified_gene of duplicated/identical genes and isoforms
dual_print( $log, "removing duplicates\n" );
merge_overlap_loci(undef, \%omniscient_modified_gene, undef, $verbose);

########
# Print results
dual_print( $log, "print intact...\n" );
print_omniscient_from_level1_id_list( {omniscient => $hash_omniscient, level_id_list =>\@intact_gene_list, output => $gffout} );

dual_print( $log, "print modified...\n" );
print_omniscient( {omniscient => \%omniscient_modified_gene, output => $gffout2} );

# create a hash containing everything
dual_print( $log, "print all with name of overlapping features resolved...\n" );
my $hash_all = subsample_omniscient_from_level1_id_list_delete($hash_omniscient, \@intact_gene_list);
merge_omniscients( $hash_all, \%omniscient_modified_gene);
merge_overlap_loci(undef, \%omniscient_modified_gene, undef, $verbose);
print_omniscient( {omniscient => $hash_all, output => $gffout3} );

#print_omniscient(\%omniscient_pseudogene, $gffout4); #print putative pseudogene in file

#END
my $string_to_print="usage: $0 @copyARGV\nCodon table used:".$codonTable."\n";
$string_to_print .="Results:\n";
$string_to_print .= "$geneCounter genes have been modified. These genes have  $mRNACounter mRNA, and we fixed the ORF of $mRNACounter_fixed of them.\n";
if (exists ($ListModel{1})){
  $string_to_print .= "$ListModel{1} model1: Prediction contains the original prediction but is longer.\n";
}
if (exists ($ListModel{2})){
  $string_to_print .= "$ListModel{2} model2: Longer prediction found non-overlapping the original one.";
  if ($split_opt){
    $string_to_print .= " Split option activated: the sequence is split in two different genes (Consequently $ListModel{2} new genes has been created)";
    }
     $string_to_print .= "\n";
}
if (exists ($ListModel{3})){
  $string_to_print .= "$ListModel{3} model3: Longer prediction found overlapping the original one but doesn't contain it(frame different).\n";
}
if (exists ($ListModel{4})) {
  $string_to_print .="$ListModel{4} model4: The prediction is shorter due to the presence of stop codon in the original CDS.\n";
}
if (exists ($ListModel{5})){
  $string_to_print .="$ListModel{5} model5: The prediction is shorter but the original CDS sequence has not premature stop codon".
  "\n  The original CDS does not start by a start codon, it is probably incomplete or fragmented".
	"\n  (The adjacent sequence of the start side might be NNNN or XXXX). ".
  "\n  The prediction is probably shorter because it is forced here to use a start codon.\n";
}
 # " The threshold to declare them as a pseudogene (comparing to the original size) is $pseudo_threshold percent.\n".
 # "According to this threshold, we change the gene status (primary_tag) of $gene_pseudo_suspected genes (corresponding to $mrna_pseudo_suspected mRNA) to pseudogene.\n".
 # "According to this threshold, we suspect $gene_pseudo_suspected genes to be pseudogenes (corresponding to $mrna_pseudo_suspected mRNA). So they habe been reported in a secpific output file.\n".
 # "$withStop_butstillgene mRNA(s) containing stop but over this treshold has been re-modelate.\n".
 # "Moreover, $mrna_pseudo_removed putative pseudo mRNA isoforms have been removed because the gene has as well non-pseudo mRNA.\n";

if (exists ($ListModel{6})){
  $string_to_print .= "$ListModel{6} model6: The prediction is the same size (AA size) but the original CDS has premature stop codons".
  " while the prediction not.\n  This is a particular case where a +1 or +2 bp shift at the beginning of the sequence".
  "\n  gives a frame shift in the original sequence but they are removed within the new prediction. \n";
}

if ($codonTable == 1){
	$string_to_print .="\n/!\\ Remind: L and M are AA are possible start codons for standard codon table (table 1).\n".
	"Particular case: If we have a triplet as WTG, AYG, RTG, RTR or ATK it will be seen as a possible start codon (but translated into X)\n";
	#"An arbitrary choisce has been done: The longer translate can begin by a L only if it's longer by 21 AA than the longer translate beginning by M. It's happened $counter_case21 times here.\n";
}
dual_print( $log, $string_to_print );
if($outfile){
  print $report $string_to_print;
}
dual_print( $log, "Bye Bye.\n" );


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

# @Purpose: Shrink CDS start related to offset (phase information)
# @input: 1 =>  list reference of CDS objects
# @output 0 => None
sub shrink_cds_offset{
	my ($sortedList) = @_;

	#set strand, check if need to be reverse complement
  my $minus = undef;
  if($sortedList->[0]->strand eq "-1" or $sortedList->[0]->strand eq "-"){ $minus = 1; }

	my $start = $sortedList->[0]->start;
	my $end = $sortedList->[$#{$sortedList}]->end;

	# in minus strand
	if($minus and $sortedList->[$#{$sortedList}]->frame != 0){
		$sortedList->[$#{$sortedList}]->end -= $sortedList->[$#{$sortedList}]->frame;
	}
	# in plus strand
	elsif (! $minus and $sortedList->[0]->frame != 0){
		$sortedList->[0]->start($sortedList->[0]->start + $sortedList->[0]->frame);
	}
}

# @Purpose: Test for ambiguous start codon (translated into X) -> we accept if the ORF is SIZE_OPT AA longer.
#           Indeed statistically it has more chance to be a real start codon.
# @input: 2 =>  object (predicted ORF), int (original protein size)
# @output 1 => Bolean
sub pass_ambiguous_start{
	my ($longest_ORF_prot_obj, $originalProt_size) = @_;
	my $pass=0;

  if( ! ( ($longest_ORF_prot_obj->seq =~ m/^X/) and ($longest_ORF_prot_obj->length() < $originalProt_size+$SIZE_OPT) ) ) {
		$pass=1;
	}
	return $pass;
}

sub modify_gene_model{

  my ($hash_omniscient, $omniscient_modified_gene, $gene_feature, $gene_id_tag_key, $level2_feature, $id_level2, $exons_features, $cds_feature_list, $cdsExtremStart, $cdsExtremEnd, $realORFstart, $realORFend, $model, $gffout)=@_;

                  ###############################################
                  # create CDS for new prediction #
                  my ($new_pred_utr5_list, $new_pred_cds_list, $new_pred_utr3_list) = modelate_utr_and_cds_features_from_exon_features_and_cds_start_stop($exons_features, $realORFstart, $realORFend);

                  #########
                  #RE-SHAPE last/first exon if less than 3 nucleotides (1  or 2 must be romved) when the CDS finish 1 or 2 nuclotide before... because cannot be defined as UTR
                  shape_exon_extremity($exons_features, $new_pred_cds_list);

                  # Create UTR
                  my $variable_not_needed;
                  ($new_pred_utr5_list, $variable_not_needed, $new_pred_utr3_list) = modelate_utr_and_cds_features_from_exon_features_and_cds_start_stop($exons_features, $realORFstart, $realORFend);

                  #############################################################
                  #  Remove ancient cds
                  my @tag_list=('exon');
                  my @id_list=($id_level2);
                  remove_tuple_from_omniscient(\@id_list, $hash_omniscient, 'level3', 'false', \@tag_list);

                  ####################
                  # Add new CDS/UTRs
                  foreach my $cds_feature (@$new_pred_cds_list){
                    push (@{$hash_omniscient->{'level3'}{'cds'}{$id_level2}}, $cds_feature);
                  }
                  foreach my $utr5_feature (@$new_pred_utr5_list){
                    push (@{$hash_omniscient->{'level3'}{'five_prime_utr'}{$id_level2}}, $utr5_feature);
                  }
                  foreach my $utr3_feature (@$new_pred_utr3_list){
                    push (@{$hash_omniscient->{'level3'}{'three_prime_utr'}{$id_level2}}, $utr3_feature);
                  }
                  $level2_feature->add_tag_value('orfix', $model);

                  check_start_end_of_mrna_feature($level2_feature, $exons_features);
                  check_start_end_of_gene_feature($hash_omniscient, $gene_id_tag_key);

                  #transfert the gene and sub-feature to the omniscient_modified_gene hash
                  my @level1_list=($gene_id_tag_key);
                  fill_omniscient_from_other_omniscient_level1_id(\@level1_list, $hash_omniscient, $omniscient_modified_gene); # If already exists in omniscient_modified_gene, it will be replaced by the modified one

}
############ /!\
# P.S: To be perfect, when a gene is newly created, we should verify if it is not created where another one has already been created. If yes, they should be linked together !!
############
sub split_gene_model{

  my ($intact_gene_list, $hash_omniscient, $omniscient_modified_gene, $gene_feature, $gene_id_tag_key, $level2_feature, $id_level2, $exons_features, $cds_feature_list, $cdsExtremStart, $cdsExtremEnd, $realORFstart, $realORFend, $model, $gffout)=@_;

  my $numberOfNewGene=1;

  my @values=$gene_feature->get_tag_values('ID');
  my $realGeneName=shift(@values);

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
  }
  else{ # ($realORFend < $cdsExtremStart)
    $bolean_original_is_first="false";
    $first_end=$realORFend;
    $second_start=$cdsExtremStart;
  }
  my ($newOrignal_exon_list, $newPred_exon_list) = create_two_exon_lists($exons_features,$first_end,$second_start,$bolean_original_is_first);

  ####################################
  # Remodelate ancient gene
  ####################################

  #############################################################
  #  Remove all level3 feature execept cds
  my @tag_list=('cds');
  my @id_list=($id_level2);
  remove_tuple_from_omniscient(\@id_list, $hash_omniscient, 'level3', 'false', \@tag_list);
  #############
  # Recreate original exon
  @{$hash_omniscient->{'level3'}{'exon'}{$id_level2}}=@$newOrignal_exon_list;

  #########
  #RE-SHAPE last/first exon if less than 3 nucleotides (1  or 2 must be romved) when the CDS finish 1 or 2 nuclotide before... because cannot be defined as UTR
  shape_exon_extremity($newOrignal_exon_list,$cds_feature_list);

  ########
  # calcul utr
  my ($original_utr5_list, $variable_not_needed, $original_utr3_list) = modelate_utr_and_cds_features_from_exon_features_and_cds_start_stop($newOrignal_exon_list, $cdsExtremStart, $cdsExtremEnd);

  #########
  #RE-SHAPE mrna extremities
  check_start_end_of_mrna_feature($level2_feature, $newOrignal_exon_list);
  $level2_feature->add_tag_value('orfix',$model);
  #########
  #RE-SHAPE gene model
  if (must_be_a_new_gene($hash_omniscient, $gene_id_tag_key, $id_level2, $level2_feature)){
    ## create a new gene
    my $new_gene_id="new_".$realGeneName."-".$numberOfNewGene;
    $numberOfNewGene++;
    my $new_gene_feature = Bio::SeqFeature::Generic->new(-seq_id => $level2_feature->seq_id, -source_tag => $level2_feature->source_tag, -primary_tag => 'gene' , -start => $level2_feature->start,  -end => $level2_feature->end, -frame => $level2_feature->frame, -strand => $level2_feature->strand , -tag => { 'ID' => $new_gene_id }) ;
    create_or_replace_tag($level2_feature,'Parent',$new_gene_id);

    # append new gene in omniscient_modified_gene
    my @level1_list=($new_gene_feature);
    my @level2_list=($level2_feature);
    my @level3_list=(@$newOrignal_exon_list, @$cds_feature_list, @$original_utr5_list, @$original_utr3_list);
    append_omniscient($omniscient_modified_gene, \@level1_list, \@level2_list, \@level3_list);
  }
  else{ # keep the original gene model that we modified
    # check shape of original gene
    check_start_end_of_gene_feature($hash_omniscient, $gene_id_tag_key);

    #include UTRS
    if ( @$original_utr5_list ){
      $hash_omniscient->{'level3'}{$original_utr5_list->[0]->primary_tag()}{$id_level2}=[@$original_utr5_list];
    }
    if ( @$original_utr3_list ){
      $hash_omniscient->{'level3'}{$original_utr3_list->[0]->primary_tag()}{$id_level2}=[@$original_utr3_list];
    }

    # append gene modified in omniscient_modified_gene
    my @level1_list=($gene_id_tag_key);
    fill_omniscient_from_other_omniscient_level1_id(\@level1_list, $hash_omniscient, $omniscient_modified_gene); # If already exists in omniscient_modified_gene, it will be replaced by the modified one
  }

  ###################################
  # Remodelate New Prediction
  ###################################

  ###############################################
  # Create CDS #
  my ($new_pred_utr5_list, $new_pred_cds_list, $new_pred_utr3_list) = modelate_utr_and_cds_features_from_exon_features_and_cds_start_stop($newPred_exon_list, $realORFstart, $realORFend);

  ####################################
  #RE-SHAPE last/first exon if less than 3 nucleotides (1  or 2 must be romved) when the CDS finish 1 or 2 nuclotide before... because cannot be defined as UTR
  shape_exon_extremity($newPred_exon_list, $new_pred_cds_list);

  #create UTR
  ($new_pred_utr5_list, $variable_not_needed, $new_pred_utr3_list) = modelate_utr_and_cds_features_from_exon_features_and_cds_start_stop($newPred_exon_list, $realORFstart, $realORFend);

  ######################################################
  # Modelate gene and mRNA features for new prediction #
  @values = $newPred_exon_list->[0]->get_tag_values('Parent');
  my $transcript_id = shift @values;
  my $new_mRNA_feature = Bio::SeqFeature::Generic->new(-seq_id => $newPred_exon_list->[0]->seq_id, -source_tag => $newPred_exon_list->[0]->source_tag, -primary_tag => 'mRNA' , -start => $newPred_exon_list->[0]->start,  -end => $newPred_exon_list->[$#{$newPred_exon_list}]->end, -frame => $newPred_exon_list->[0]->frame, -strand => $newPred_exon_list->[0]->strand , -tag => { 'ID' => $transcript_id , 'Parent' => $realGeneName }) ;

  my @level1_list;
  my @level2_list;
  my @level3_list;

  #$numberOfNewGene == 1 mean we already use the gene container. So in the case where we have oly one mRNA, the split will give 2 mRNA. One is linked to the original gene container (done before)
  # The second must be linked to a new gene container. So, even if must_be_a_new_gene method say no, we must create it because the original one has been already used.
  my $create_a_new_gene=must_be_a_new_gene($hash_omniscient, $gene_id_tag_key, $transcript_id, $new_mRNA_feature);
  if ( ($#{$hash_omniscient->{'level2'}{'mrna'}{$gene_id_tag_key}} == 0) and $numberOfNewGene == 1){ $create_a_new_gene="true";}
  if ( $create_a_new_gene ){
    my $new_gene_id="new_".$realGeneName."-".$numberOfNewGene;
    create_or_replace_tag($new_mRNA_feature, 'Parent', $new_gene_id);
    my $new_gene_feature = Bio::SeqFeature::Generic->new(-seq_id => $newPred_exon_list->[0]->seq_id, -source_tag => $newPred_exon_list->[0]->source_tag, -primary_tag => 'gene' , -start => $newPred_exon_list->[0]->start,  -end => $newPred_exon_list->[$#{$newPred_exon_list}]->end, -frame => $newPred_exon_list->[0]->frame, -strand => $newPred_exon_list->[0]->strand , -tag => { 'ID' => $new_gene_id } , 'orfix' => $model) ;
    @level1_list=($new_gene_feature);
    @level2_list=($new_mRNA_feature);
  }
  else{ #the new mRNA still overlap an isoform. So we keep the link with the original gene
    # append new gene in omniscient_modified_gene
    check_start_end_of_gene_feature($hash_omniscient, $gene_id_tag_key);
    @level1_list=($gene_feature);
    @level2_list=($new_mRNA_feature);
  }
  @level3_list=(@$newPred_exon_list, @$new_pred_cds_list, @$new_pred_utr5_list, @$new_pred_utr3_list);
  append_omniscient($omniscient_modified_gene, \@level1_list, \@level2_list, \@level3_list); # If already exists , no replacement

  if ($numberOfNewGene > 1){
    #remove the mRNA from original omnicient (because the two mRNAs form the splited one are no linked to the original gene
    # but are now linked to newly created gene features). The same for all linked level 3 features, so we remove them.
    my @tag_list=('all');
    my @id_list=($id_level2);

    remove_tuple_from_omniscient(\@id_list, $hash_omniscient, 'level3', 'false', \@tag_list);
    @id_list=($gene_id_tag_key);my @id_list2=($id_level2);
    remove_element_from_omniscient(\@id_list, \@id_list2, $hash_omniscient, 'level2', 'false', \@tag_list);
    #reshape end and start
    check_start_end_of_gene_feature($hash_omniscient, $gene_id_tag_key);
    push(@{$intact_gene_list}, $gene_id_tag_key);
  }
}

# Yes if mRNA doesnt overlap an other existing isoform
sub must_be_a_new_gene{
  my ($hash_omniscient, $gene_id, $id_level2, $level2_feature)=@_;

  my $result="true";
  my @list_mrna=@{$hash_omniscient->{'level2'}{'mrna'}{$gene_id}};

  if($#list_mrna > 0){ #more than only one mrna
    foreach my $mrna (@list_mrna){
      # get level2 id
      my @values = $mrna->get_tag_values('ID');
      my $mrna_id = lc(shift @values);
      if(lc($id_level2) ne lc($mrna_id)){ # we dont check mrna against itself
        #Now check if overlap
        if( ($level2_feature->start <= $mrna->end) and ($level2_feature->end >= $mrna->start) ){ # if it overlaps
          $result=undef;last;
        }
      }
    }
  }
  else{$result=undef;} #only one mRNA

  return $result;
}

sub shape_exon_extremity{
  #exon_features is a sorted list
  #cds_features is a sorted list

  my ($exon_features,$cds_features)=@_;

   #test between first exon and first cds
   if( (abs($cds_features->[0]->start - $exon_features->[0]->start) < 3) and (abs($cds_features->[0]->start - $exon_features->[0]->start) > 0) ){ #We have to shape the exon start. We don't want a non multiple of 3 inferior to 3

      $exon_features->[0]->start($cds_features->[0]->start);
#      print "start reshaped\n";
   }
   #test between last exon and last cds
   if(abs($exon_features->[$#{ $exon_features }]->end - $cds_features->[$#{ $cds_features }]->end ) < 3){  #We have to shape the exon end
      $exon_features->[$#{ $exon_features }]->end($cds_features->[$#{ $cds_features }]->end);
#      print "end reshaped\n";
   }
}

sub calcul_real_orf_end_and_start{
  #exons_features is sorted
  my ($orf_cds_region, $exons_features)=@_;

  my $realORFstart;
  my $realORFend;

  my $orf_start=$orf_cds_region->[0]; # get start to begin
  my $orf_length=$orf_cds_region->[2]; # get start to begin

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

# Check the start and end of mRNA and gene feature;
sub check_start_end_of_mrna_feature{

  my ($mRNA_feature, $exon_list)=@_;

  ######
  #Modify mRNA start-end based on exon features
  my $exonStart=$exon_list->[0]->start;
  my $exonEnd=$exon_list->[$#{$exon_list}]->end;
  if ($mRNA_feature->start != $exonStart){
    $mRNA_feature->start($exonStart);
  }
  elsif($mRNA_feature->end != $exonEnd){
    $mRNA_feature->end($exonEnd);
  }
}

# Check the start and end of gene feature based on its mRNA;
sub check_start_end_of_gene_feature{

  my ($hash_omniscient, $gene_id)=@_;

  #####
  #Modify gene start-end (have to check size of each mRNA)
  my $geneExtremStart=1000000000000;
  my $geneExtremEnd=0;
  foreach my $primary_tag_key_level2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
    foreach my $mrna_feature ( @{$hash_omniscient->{'level2'}{'mrna'}{$gene_id}}) {
      my $start=$mrna_feature->start();
      my $end=$mrna_feature->end();

      if ($start < $geneExtremStart){
        $geneExtremStart=$start;
      }
      if($end > $geneExtremEnd){
        $geneExtremEnd=$end;
      }
    }
  }
  my $gene_feature=$hash_omniscient->{'level1'}{'gene'}{$gene_id};
  if ($gene_feature->start != $geneExtremStart){
      $gene_feature->start($geneExtremStart);
   }
   elsif($gene_feature->end != $geneExtremEnd){
      $gene_feature->end($geneExtremEnd);
    }
}


# The exons containing the original cds keep their parent names. The exon containing the new cds will have a new parent name.
sub create_two_exon_lists {
  # orignalFirst == true if original gene is first on the prediction
  my ($exons_features,$firstEnd, $secondStart, $orignalFirst)=@_;
  my @list_exon_originalPred;
  my @list_exon_newPred;

  foreach my $exon_feature (@$exons_features){ #for each exon
#    print "start:".$exon_feature->start." end:".$exon_feature->end."\n";
    if(two_positions_on_feature($exon_feature,$firstEnd,$secondStart)){ # We have to split the exon_feature P.S: We will loss sequence between the two positions
#      print "both on feature\n";
      my $duplicated_exon_feature=clone($exon_feature);#create a copy of the feature
      #manage original exon
      $exon_feature->end($firstEnd);
      $duplicated_exon_feature->start($secondStart);

      if($orignalFirst eq "true"){
        push( @list_exon_originalPred, $exon_feature);

        my @values = $duplicated_exon_feature->get_tag_values('ID');
        my $value = $values[0];
        create_or_replace_tag($duplicated_exon_feature,'ID', 'new_'.$value);
        @values = $duplicated_exon_feature->get_tag_values('Parent');
        $value = $values[0];
        create_or_replace_tag($duplicated_exon_feature,'Parent', 'new_'.$value);
        push( @list_exon_newPred, $duplicated_exon_feature);
        next;
      }else{ #original pred after
        push( @list_exon_originalPred, $duplicated_exon_feature);

        my @values = $exon_feature->get_tag_values('ID');
        my $value = $values[0];
        create_or_replace_tag($exon_feature,'ID', 'new_'.$value);
        @values = $exon_feature->get_tag_values('Parent');
        $value = $values[0];
        create_or_replace_tag($exon_feature,'Parent', 'new_'.$value);
        push( @list_exon_newPred, $exon_feature);
        next;
      }
    }
    if(! (($exon_feature->end <=  $secondStart) and ($exon_feature->start >=  $firstEnd))){ # We remove it because exon between CDSs
      if ($exon_feature->end <=  $secondStart) {
        if ($orignalFirst eq "true"){
          push( @list_exon_originalPred, $exon_feature);
        }else{
          my $duplicated_exon_feature=clone($exon_feature);#create a copy of the feature
          my @values = $duplicated_exon_feature->get_tag_values('ID');
          my $value = $values[0];
          create_or_replace_tag($duplicated_exon_feature,'ID', 'new_'.$value);
          @values = $duplicated_exon_feature->get_tag_values('Parent');
          $value = $values[0];
          create_or_replace_tag($duplicated_exon_feature,'Parent', 'new_'.$value);
          push( @list_exon_newPred, $duplicated_exon_feature);
        }
      }
      if ($exon_feature->start >=  $firstEnd) {
        if($orignalFirst eq "true"){
            my $duplicated_exon_feature=clone($exon_feature);#create a copy of the feature
            my @values = $duplicated_exon_feature->get_tag_values('ID');
            my $value = $values[0];
            create_or_replace_tag($duplicated_exon_feature,'ID', 'new_'.$value);
            @values = $duplicated_exon_feature->get_tag_values('Parent');
            $value = $values[0];
            create_or_replace_tag($duplicated_exon_feature,'Parent', 'new_'.$value);
            push( @list_exon_newPred, $duplicated_exon_feature);
          }else{
            push( @list_exon_originalPred, $exon_feature);
          }
        }
      }
  }
  my @list_exon_originalPred_sorted = sort {$a->start <=> $b->start} @list_exon_originalPred;
  my @list_exon_newPred_sorted = sort {$a->start <=> $b->start} @list_exon_newPred;

  return \@list_exon_originalPred_sorted, \@list_exon_newPred_sorted;
}

sub position_on_feature {

  my ($feature,$position)=@_;

  my $isOnSameExon=undef;
  if ( ($position >= $feature->start and $position <= $feature->end)){
    $isOnSameExon="true";
  }
  return $isOnSameExon;
}

sub two_positions_on_feature {

  my ($feature,$position1,$position2)=@_;

  my $areOnSameExon=undef;
  if ( ($position1 >= $feature->start and $position1 <= $feature->end) and ($position2 >= $feature->start and $position2 <= $feature->end) ){
    $areOnSameExon="true";
  }
  return $areOnSameExon;
}

sub translate_JD {
   my ($self,@args) = @_;
     my ($terminator, $unknown, $frame, $codonTableId, $complete,
     $complete_codons, $throw, $codonTable, $orf, $start_codon, $no_start_by_aa, $offset);

   ## new API with named parameters, post 1.5.1
   if ($args[0] && $args[0] =~ /^-[A-Z]+/i) {
         ($terminator, $unknown, $frame, $codonTableId, $complete,
         $complete_codons, $throw,$codonTable, $orf, $start_codon, $no_start_by_aa, $offset) =
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
                                               NOSTARTBYAA
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
            ($orf_region) = _find_orfs_nucleotide_JD( $self, $seq, $codonTable, $start_codon, $no_start_by_aa, $orf eq 'longest' ? 0 : 'first_only' );
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

sub _find_orfs_nucleotide_JD {
    my ( $self, $sequence, $codon_table, $start_codon, $no_start_by_aa, $first_only ) = @_;
    $sequence    = uc $sequence;
    $start_codon = uc $start_codon if $start_codon;

    my $is_start = $start_codon
        ? sub { shift eq $start_codon }
        : sub { $codon_table->is_start_codon( shift ) };

    # stores the begin index of the currently-running ORF in each
    # reading frame
    my @current_orf_start = (-1,-1,-1);

    #< stores coordinates of longest observed orf (so far) in each
    #  reading frame
    my @orfs;

    # go through each base of the sequence, and each reading frame for each base
    my $seqlen = CORE::length $sequence;
    for( my $j = 0; $j <= $seqlen-3; $j++ ) {
        my $frame = $j % 3;

        my $this_codon = substr( $sequence, $j, 3 );
        my $AA = $codon_table->translate($this_codon);

        # if in an orf and this is either a stop codon or the last in-frame codon in the string
        if ( $current_orf_start[$frame] >= 0 ) {
            if ( $codon_table->is_ter_codon( $this_codon ) ||( my $is_last_codon_in_frame = ($j >= $seqlen-5)) ) {
                # record ORF start, end (half-open), length, and frame
                my @this_orf = ( $current_orf_start[$frame], $j+3, undef, $frame );
                my $this_orf_length = $this_orf[2] = ( $this_orf[1] - $this_orf[0] );

                $self->warn( "Translating partial ORF "
                                 .$self->_truncate_seq( $self->_orf_sequence( $sequence,\@ this_orf ))
                                 .' from end of nucleotide sequence'
                            )
                    if $first_only && $is_last_codon_in_frame;

                return\@ this_orf if $first_only;
                push @orfs,\@ this_orf;
                $current_orf_start[$frame] = -1;
            }
        }
        # if this is a start codon
        elsif ($is_start->($this_codon)) {
          if($no_start_by_aa){

            if($AA ne $no_start_by_aa){
              $current_orf_start[$frame] = $j;
            }
          }
          else{
            $current_orf_start[$frame] = $j;
          }
        }
    }

    return sort { $b->[2] <=> $a->[2] } @orfs;
}

# Sort by locusID !!!!
# L1 => LocusID->level->typeFeature->ID =[ID,start,end]
# L2 and L3 => LocusID->level->typeFeature->Parent->ID = [ID,start,end]
#
#
sub _sort_by_seq{
  my ($omniscient) = @_;

  my %hash_sortBySeq;

    foreach my $tag_level1 (keys %{$omniscient->{'level1'}}){
      foreach my $level1_id (keys %{$omniscient->{'level1'}{$tag_level1}}){
        my $level1_feature = $omniscient->{'level1'}{$tag_level1}{$level1_id};
        my $ID = $level1_feature->_tag_value('ID');
        my $strand="+";
        if($level1_feature->strand != 1){$strand = "-";}
        my $position_l1=$level1_feature->seq_id."".$strand;

        $hash_sortBySeq{$position_l1}{"level1"}{$tag_level1}{$level1_id} = [$ID, int($level1_feature->start), int($level1_feature->end)];
        }
  }
  return \%hash_sortBySeq;
}

# Previous description with pseudogene =>
#one contains the putative pseudogene detected (As they are just putatuve, they are also present among the intacts ),
#and a last a report of the results.
#Pseudogene particularity: If gene contains mRNA models goods and mRNA that look like a pseudogene, the pseudogene one will be removed.

__END__

=head1 NAME

agat_sp_fix_longest_ORF.pl

=head1 DESCRIPTION

The script aims to fix the ORFs of gene models described in the gff file.
By fixing it means replacing the original ORF (defined by the cds)
when the longest predicted one within the mRNA is different. See the --model parameter
for more details about the different cases. Currently the tool does not perform
incomplete prediction (It always look for a start codon). It is consequently advised
to not use the model5 except if you understand what you do.
Several ouput files will be written if you specify an output.
One will contain the gene not modified (intact), one with the gene models fixed (modified),
one will both together (all).

=head1 SYNOPSIS

    agat_sp_fix_longest_ORF.pl -gff infile.gff --fasta genome.fa [ -o outfile ]
    agat_sp_fix_longest_ORF.pl --help

=head1 OPTIONS

=over 8

=item B<--gff>

Input GTF/GFF file.

=item B<-f>, B<--fa> or B<--fasta>

Imput fasta file.

=item B<--ct>, B<--codon> or B<--table>

Codon table to use. [default 1]

=item B<-m> or B<--model>

Kind of ORF Model you want to fix. By default all are used. To select specific models writte e.g --model 1,4

Model1 = The original ORF is part of the new ORF; the new ORF is longer

Model2 = The original ORF and the new one are different; the new one is longer, they do not overlap each other.

Model3 = The original ORF and the new one are different; the new one is longer, they overlap each other.

Model4 = The new ORF is shorter due to the presence of stop codon in the original ORF.

Model5 = The new ORF is shorter but the original ORF has not premature stop codon.
         The shorter predicted ORF can be due to the fact that the original ORF does not start by a start codon,
				 while we force here the prediction to have a start codon.
				 A ORF wihtout start can be the fact of an incomplete or fragmented ORF:
				 annotation tool didn't predict the start because:
				 * the start region is NNNN
				 * the start region is XXXX
				 * correct nucleotides but prediction tool did not annotate this part (e.g incomplete evidence in evidence-based prediction)

Model6 = The ORF is same size but not correct frame (+1 or +2 bp gives a frame shift).

=item B<-s> or B<--split>

This option is usefull for Model2. Indeed when the prediction is non overlapping the original cds, it is possible to split the gene into two different genes. By default we don't split it.
We keep the longest. If you want to split it type: -s

=item B<-o> , B<--output> , B<--out> or B<--outfile>

Output GFF file.  If no output file is specified, the output will be
written to STDOUT.

=item B<-v>

verbose mode. Default off. -v 1 minimum verbosity, -v 3 maximum verbosity

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


MODEL6 example (extra C at the begining of the CDS):
sequence= CATGTTCAAACGTCTCGAAAATATAGCCGTCCAATCATCCTCCTTCCCCCAGGCAATCTCCTTGATCCAGCAAAACCACCTCTCTCCAAAACTCTTCTTTGATCCCCAGACCTACTCCAAGATCTTCCAGAAACTCTCCCTCAAAGACCAATACCCTCGTTCCTCCCAGTCCCTATGCATCATAGACTACCACGTAGGCTTCACGCCCTTCTCCTACTTCCTCCATAAGGAGCTACACCCTGCACATCACGTCATCTTCCCCGATAGTGTCGCTGCCAACAAGTTCTGGACCCAGATGCTACTGAAGGACCCCGACTGTAAGGACATGGTCATAGACGAGACTCAGGCAAACACAGTGCTAAAGCACAATTTCCTCAATAGATCGCTGGAATTGGGCCACGTCGTTGCAGTAGAACAGACAGACCTAACTAAGGTCAACGACTCGATACTATTGACCGGTAACTTCGTCGATACTTCCGGCGGGGACTCTCTACGGATCTTACTCTTCTTTAATCAGATGAAAACTTCCGTCTTTCAGTATAATAACGTCAAGTTCTTGGCGTGGCTGCCCGCCTCGGAGTCTCTGAAGTTCATAGGACCGATTGGATCGAGGCATAGACGGTCCAATGCGCTGATGACCAACCTATTTGCCAACGTTGACGTGGTAGCGTACTCTAATTATGGCAAGAAGAAGAGCGTTTCCCGAGTCTTGGACGAATATAAGGACGCTGTCAAGCTACCACAGATTCCTGGACAGAAAGACGTATGTTTGATGGAATTTCAGTCGAACTATTCCAAATACGACATTAAATTTCCTGACGAATTGCATTTGATCATACACAAGATGTTGATATCGCCCAGCAATAAGTTGATTGACAATCTTCATTTGCTTGGGCCCGGTGCAAAGGAGTATTTGGAGACCAAGTTGGATCCCGAGCTGTTACAGAAGCCTGCGCCGAACATTACGGAGCAGGAGTTTATAGATATCAGCGAGGCGTATTATTACTGGCCGTTCAAGCCTAACGTTCACTTGGAGACGTATTTAGGAGATCCTCCGGAGGAGGAGTAG
GFF =
y922_scaffold13 . gene  1  1069  . + . ID=DEKNAAG101268;Name=DEKNAAG101268
y922_scaffold13 . mRNA  1  1069  . + . ID=DEKNAAT101273;Parent=DEKNAAG101268;Name=DEKNAAT101273;description=Predicted: mitochondrial rna polymerase specificity factor
y922_scaffold13 . exon  1  1069  . + . ID=DEKNAAE101408;Parent=DEKNAAT101273;Name=DEKNAAE101408
y922_scaffold13 . CDS 1  1069  . + 0 ID=DEKNAAC101407;Parent=DEKNAAT101273;Name=DEKNAAC101407
