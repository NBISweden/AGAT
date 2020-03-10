#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use File::Basename;
use Getopt::Long;
use Statistics::R;
use IO::File;
use Pod::Usage;
use List::MoreUtils qw(uniq);
use Clone 'clone';
use AGAT::Omniscient;

my $header = get_agat_header();

#####
# What we call parial gene (containing "_partial_part-" in the ID) ?
# This gene has been seen as patial: During a lift-over a gene can be detected on 2 several contigs.
# (full kraken file => features with kraken attribute to TRUE are on contig of the target genome (Transfert annotation on), the others (kraken attribute to FALSE) are on the reference genome to liftfover (where annotations are taken to try to liftover) )
#####
#
# TODO: add tag kraken_cn (for copy number) with nb of mapping. 1 if only one liftover.
# Ask Manfred why some region map at different location.

my $outfile = undef;
my $gff = undef;
my $valueK = undef;
my $verbose = undef;
my $kraken_tag = "Kraken_mapped";
my $help= 0;

if ( !GetOptions(
    "help|h" => \$help,
    "gtf=s" => \$gff,
    "threshold|t=i" => \$valueK,
    "verbose|v!" => \$verbose,
    "outfile|output|out|o=s" => \$outfile))

{
    pod2usage( { -message => 'Failed to parse command line',
                 -verbose => 1,
                 -exitval => 1 } );
}

# Print Help and exit
if ($help) {
    pod2usage( { -message => "$header\n",
                 -verbose => 99,
                 -exitval => 0 } );
}

if ( ! (defined($gff)) ){
    pod2usage( {
           -message => "$header\nAt least 1 parameter is mandatory:\nInput reference gtf file (--f)\n\n",
           -verbose => 0,
           -exitval => 1 } );
}

## Manage output file
my $gffout;
my $outReport     = IO::File->new();
if ($outfile) {
  $outfile=~ s/.gff//g;

  open(my $fh, '>', $outfile.".gff") or die "Could not open file '$outfile' $!";
  $gffout= Bio::Tools::GFF->new(-fh => $fh, -gff_version => 3 );

  $outReport->open($outfile."_report.txt", 'w') or die "Could not open file '$outfile'_report.txt $!";
}
else{
  $gffout = Bio::Tools::GFF->new(-fh => \*STDOUT, -gff_version => 3);

  $outReport->fdopen( fileno(STDOUT), 'w' ) or die "Could not open file STDOUT $!";
}

# Message
my $messageValue;
if ( defined($valueK) ){
  $messageValue = "You choose to keep in output only genes mapped over $valueK percent.\n"
}else{
  $messageValue = "We will keep all the mapped features.\n";
  $valueK=0;
}
$messageValue.="The kraken attribute tag that will be used is: ".$kraken_tag."\n";

#print info
if ($outfile) {
  print $outReport $messageValue;
  }else{print $messageValue;}

                #####################
                #     MAIN          #
                #####################


######################
### Parse GFF input #
# checks are deactivated except _remove_orphan_l1
my $verbose_omniscient = -1 if (! $verbose);
my ($hash_omniscient, $hash_mRNAGeneLink) =  slurp_gff3_file_JD({
                                                               input => $gff,
                                                               no_check => 1,
                                                               no_check_skip => ["_remove_orphan_l1"],
																															 verbose =>  2
                                                               });

#track stats
my $nbOriginalGene = nb_feature_level1($hash_omniscient); #total gene at the beginning
my $nbRecordMap=0; #total original gene that map before check the threshold
my $nbOriginalL0Map=0; #total map before check the threshold
my $nbOriginal_multiMap_seqdif=0; #number of gene that have multimap on seq location before check the threshold
my $nbOriginal_total_multiMap_seqdif=0; #total of multimap on different seq before check the threshold
my $nbOriginal_multiMap_sameseq=0;  #number of gene that have multimap on same seq after check the threshold
my $nbOriginal_total_multiMap_sameseq=0; #total of multimap on different seq after check the threshold

#my $nbMapL1=0; #total gene that map after check the threshold
#my $nbMapL1Uniq=0;
my $nbGeneIdUniqMap=0;
my $nb_multiMap_seqdif=0;  #number of gene that have multimap on different seq after check the threshold
my $nb_total_multiMap_seqdif=0; #total of multimap on different seq after check the threshold
my $nb_total_multiMap_seqdif_bothcase=0;
my $nb_total_multiMap_sameseq_bothcase=0;
my $nb_multiMap_sameseq=0;  #number of gene that have multimap on same seq after check the threshold
my $nb_total_multiMap_sameseq=0; #total of multimap on different seq after check the threshold
my $bothCase=0;

# track errors:
my $KrakenFakeGene=0;


my %mappedPercentPerGene; #Keep information for R plot
my %n_omniscient;
my $nb_noCaseL3=0;
my $new_omniscient=\%n_omniscient;
my $list_uID_new_omniscient;
my $loop=0;



# sort by seq id
my $hash_sortBySeq = gather_and_sort_l1_by_seq_id($hash_omniscient);

#################
# == LEVEL level1
################
foreach my $seqid (sort { (($a =~ /(\d+)$/)[0] || 0) <=> (($b =~ /(\d+)$/)[0] || 0) } keys %{$hash_sortBySeq}){ # loop over all the feature level1


		foreach my $primary_tag_key_level1 (sort {$a cmp $b} keys %{$hash_omniscient->{'level1'}}){
			foreach my $gene_feature ( @{$hash_sortBySeq->{$seqid}{$primary_tag_key_level1}} ){
				my $id_tag_key_level1 = lc($gene_feature->_tag_value('ID'));


	    ########################################
	    # Prepare hash in case of muli mapping #
	    ########################################
	    my %listOfProperHash;
	    my %listHashWithTrue;
	    my $l1_original_id = $id_tag_key_level1;
	    #################
	    # == LEVEL 1 == #
	    #################
	    $listOfProperHash{$gene_feature->seq_id()}{'level1'}{$primary_tag_key_level1}{$id_tag_key_level1}=$gene_feature;
	    # write down if kraken_mapped=true
	    if($gene_feature->has_tag($kraken_tag)){
	      if( lc($gene_feature->_tag_value($kraken_tag)) eq "true"){
	        $listHashWithTrue{$gene_feature->seq_id()}++;
	      }
	    }
	    #################
	    # == LEVEL 2 == #
	    #################
	    foreach my $primary_tag_key_level2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
	      if ( exists_keys($hash_omniscient, ('level2',$primary_tag_key_level2,$id_tag_key_level1) ) ){
	        foreach my $feature_level2 ( @{$hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1}}) {
	          push(@{$listOfProperHash{$feature_level2->seq_id()}{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1}},$feature_level2) ;
	          my $level2_ID = lc($feature_level2->_tag_value('ID'));
	          # write down if kraken_mapped=true
	          if($feature_level2->has_tag($kraken_tag)){
	            if( lc($feature_level2->_tag_value($kraken_tag)) eq "true"){
	              $listHashWithTrue{$feature_level2->seq_id()}++;
	            }
	          }
	          #################
	          # == LEVEL 3 == #
	          #################
	          foreach my $primary_tag_l3  (keys %{$hash_omniscient->{'level3'}}){ # primary_tag_l3 = cds or exon or start_codon or utr etc...
	            if ( exists_keys( $hash_omniscient, ('level3', $primary_tag_l3, $level2_ID) ) ){
	              foreach my $feature_level3 ( @{$hash_omniscient->{'level3'}{$primary_tag_l3}{$level2_ID}}) {
	                push(@{$listOfProperHash{$feature_level3->seq_id()}{'level3'}{$primary_tag_l3}{$level2_ID}}, $feature_level3);
	                # write down if kraken_mapped=true
	                if($feature_level3->has_tag($kraken_tag)){
	                  if( lc($feature_level3->_tag_value($kraken_tag)) eq "true"){
	                    $listHashWithTrue{$feature_level3->seq_id()}++;

	                  }
	                }
	              }
	            }
	          }
	        }
	      }
	    }

	    #count multi map
	    my $nbMapTrueHere = keys %listHashWithTrue;
	    if($nbMapTrueHere > 1){
	      $nbOriginal_multiMap_seqdif++;
	      $nbOriginal_total_multiMap_seqdif+=$nbMapTrueHere;
	    }


	    my $sucessMapL0=0;
	    my $sucessMapL1OusideScope=0;
	    my $firstL0Map="yes";


	    # A record contains 1 or several LEVEL 0
	    ################
	    # == LEVEL 0
	    ################
	    foreach my $seqid_key ( sort keys %listOfProperHash){

	      #if it contains a feature mapped we continue
	      if(exists_keys(\%listHashWithTrue,$seqid_key)){

	        #keep track of statistics
	        if($firstL0Map){ # first record of
	          $nbRecordMap++;
	          $firstL0Map=undef;
	        }
	        $nbOriginalL0Map++;

	        #The HAsh is made only by feature with TRUE (because the false one are not collected because different)
	        my $hash = clone($listOfProperHash{$seqid_key}); # We have to clone it otherwise the "merge_omniscients" function can change the id and brakes the original data from hash_omniscient that is used to retrieve the original size of the feature
	        if(! exists_keys($hash,('level3'))){
	          $nb_noCaseL3++;
	        }

					# We skip _check_exons and _check_utrs to not fit the exon to the old mRNA size that was making big last or first exon
	        my ($hash_omniscient_clean, $hash_mRNAGeneLink_clean) = slurp_gff3_file_JD({ input => $hash,
																																											 verbose => 2,
																																											 no_check => 1,
												                                                               no_check_skip => ["_check_sequential",
																																																				 "_check_l2_linked_to_l3",
																																																				 "_check_l1_linked_to_l2",
																																																				 "_remove_orphan_l1",
																																																				 "_check_all_level2_positions",
																																																				 "_check_all_level1_positions"],
	                                                                                   });

	        if($verbose){
	          print "\nA proper hash:\n";
	          print_omniscient($hash_omniscient_clean, $gffout);
	          print "\n";
	        }

	        ###################################################################################
	        # NOW we call deal properly with each proper hash containing only mapped features
	        ## ################################################################################
	        ################
	        # == LEVEL 1
	        ################
	        my $sucessMapL1=0;
	        foreach my $primary_tag_key_level1 (keys %{$hash_omniscient_clean->{'level1'}}){ # primary_tag_key_level1 = gene or repeat etc...
	          foreach my $id_tag_key_level1 (keys %{$hash_omniscient_clean->{'level1'}{$primary_tag_key_level1}}){


	            $gene_feature = $hash_omniscient_clean->{'level1'}{$primary_tag_key_level1}{$id_tag_key_level1};
	            my @ListmrnaNoMatch;
	            print "\n\nlevel1 feature:\n".$gene_feature->gff_string."\n\n" if $verbose;

	            ################
	            # == LEVEL 2
	            ################

	            my $sucessMapL2=0;
	            foreach my $primary_tag_key_level2 (keys %{$hash_omniscient_clean->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...

	              if ( exists_keys($hash_omniscient_clean, ('level2',$primary_tag_key_level2,$id_tag_key_level1) ) ){
	                foreach my $feature_level2 ( @{$hash_omniscient_clean->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1}}) {
	                  print "level2 feature:\n".$feature_level2->gff_string."\n" if $verbose;

	                  my $percentMatch=0;

	                  my $level2_ID = lc($feature_level2->_tag_value('ID'));

	                  ###################################################
	                  # == LEVEL 3
	                  # We will look the size of mapped features
	                  # Feature can be exon ,cds or utr in that order
	                  ###################################################

	                  my $refListFetaureL3 = takeOneListLevel3From1idLevel2($hash_omniscient_clean, $level2_ID);
	                  if(! $refListFetaureL3 ){ # No l3 found we clean this l2 feature
	                    my @listF2; push (@listF2, $feature_level2);
	                    remove_omniscient_elements_from_level2_feature_list ($hash_omniscient_clean, \@listF2);
	                    if (! exists_keys($hash_omniscient_clean,('level1',$primary_tag_key_level1,$id_tag_key_level1))){
	                      $KrakenFakeGene++;
	                    }
	                    next;
	                  }

	                  my $matchSize = 0;
	                  my $matchFeatureExample = undef;

	                  foreach my $feature (@{$refListFetaureL3}){
	                    #print "level3 feature:\n".$feature->gff_string."\n" if $verbose;
	                    my $end=$feature->end();
	                    my $start=$feature->start();

	                    my $mapping_state = undef;
	                    if($feature->has_tag($kraken_tag)){
	                      $mapping_state = lc($feature->_tag_value($kraken_tag));
	                    }
	                    else{ print "error !! No $kraken_tag attribute found for the feature".$feature->gff_string()."\n";}

	                    if( $mapping_state eq "true"){
	                        $matchSize+=($end-$start)+1;
	                        $matchFeatureExample=$feature;
	                    }
	                    elsif(! $mapping_state eq "false"){
	                      print "error !! We don't understand the $kraken_tag attribute value found for the feature".$feature->gff_string()."\n Indeed, we expect false or true.\n";
	                    }
	                  }


	                  #compute the total sie (has to be compute against the original hash)
	                  my $totalSize=0;
	                  $totalSize = compute_total_size($hash_omniscient, $l1_original_id, $matchFeatureExample);

	                  #compute the MATCH. A MATCH can be over 100% because we compute the size of the original feature l3 against the new feature l3. The new feature l3 (i.e exon) could have been strenghten to fit a new size/map of feature l2.
	                  $percentMatch=($matchSize*100)/$totalSize;
	                  print "$id_tag_key_level1 / $level2_ID  maps at ".$percentMatch." percent.\n" if $verbose;
	                  #if($percentMatch > 100){
	                  #  print $id_tag_key_level1."\n";exit;
	                  #}

	                  #######
	                  # Add information to gff
	                  ########
	                  if ($percentMatch >= $valueK) {
	                    $sucessMapL2++;
	                    #We print gene only if a percentage match value is superior to the threshold fixed (No threshold equal everything = 0)
	                    $percentMatch = sprintf('%.2f', $percentMatch);

	                    manage_gene_label($gene_feature, $percentMatch, $kraken_tag); # add info to level1 (gene) feature

	                    create_or_replace_tag($feature_level2,$kraken_tag,$percentMatch."%"); # add info to level2 (mRNA) feature
	                    create_or_replace_tag($feature_level2,'description',"Mapped at ".$percentMatch."%"); # add info to level2 (mRNA) feature

	                    #save best value for gene
	                    if(! exists_keys(\%mappedPercentPerGene, ($l1_original_id) ) ){ # case where it doesn t exist
	                      $mappedPercentPerGene{$l1_original_id}=$percentMatch;
	                      $nbGeneIdUniqMap++;
	                    }
	                    elsif($mappedPercentPerGene{$l1_original_id} < $percentMatch){ # case where it exists but better value to save
	                      $mappedPercentPerGene{$l1_original_id}=$percentMatch;
	                    }
	                  }
	                  # Do not pass the threshold we have to remove this l2 feature
	                  else{
	                    my @listF2; push (@listF2, $feature_level2);
	                    remove_omniscient_elements_from_level2_feature_list ($hash_omniscient_clean, \@listF2);
	                    if (! exists_keys($hash_omniscient_clean,('level1',$primary_tag_key_level1,$id_tag_key_level1))){
	                      $KrakenFakeGene++;
	                    }
	                    next;
	                  }
	                }
	              }
	            }
	            if($sucessMapL2){
	              $sucessMapL1++;
	            }
	          }
	        }
	        if ($sucessMapL1){ # We have a result over the threshold to save
	          $sucessMapL0++;
	          $sucessMapL1OusideScope=$sucessMapL1;
	          #$nbMapL1Uniq++;
	          #$nbMapL1 +=$sucessMapL1;
	          #save the result by appending the result hash and take care of duplicated names

	          if($loop == 0){
	            $new_omniscient = $hash_omniscient_clean;
	            $loop++;
	          }
	          elsif($loop == 1){
	            ($new_omniscient, $list_uID_new_omniscient) = merge_omniscients($new_omniscient, $hash_omniscient_clean);
	            $loop++;
	          }
	          else{
	            ($new_omniscient, $list_uID_new_omniscient) = merge_omniscients($new_omniscient, $hash_omniscient_clean, $list_uID_new_omniscient);
	          }

	          #keep track of successful multimap (same sequences) > Cases saved with different GeneID in new_omniscient ()
	          if($sucessMapL1 > 1){
	            $nb_multiMap_sameseq++;
	            $nb_total_multiMap_sameseq+=$sucessMapL1;
	          }
	        }
	      }
	    }
	    #keep track of successful multimap (different sequences) => Cases saved with different GeneID in new_omniscient
	    if($nbMapTrueHere > 1 and $sucessMapL0 > 1){
	      if( $sucessMapL1OusideScope > 1){
	        $bothCase++; print "Both case:\nNb multi map seq diff=$sucessMapL0\nNb multi map same seq =$sucessMapL1OusideScope\n" if $verbose;
	        $nb_total_multiMap_seqdif_bothcase+=$sucessMapL0;
	        $nb_total_multiMap_sameseq_bothcase+=$sucessMapL1OusideScope;
	        $nb_multiMap_sameseq--;
	        $nb_total_multiMap_sameseq-=$sucessMapL1OusideScope;
	      }
	      else{
	        $nb_multiMap_seqdif++;
	        $nb_total_multiMap_seqdif+=$sucessMapL0;
	      }
	    }
	  }
	}
}
print "Calcul of mapped percentage length finished !\n";


######################
# Check if nothing mapped
my $nbKey = keys %mappedPercentPerGene;
if ($nbKey == 0){
 print "No succefully mapped feature found!\n";
}

########
#print GFF the selected features (over the choosen treshold)
########
print_omniscient($new_omniscient, $gffout);

my $nbEndGene = nb_feature_level1($new_omniscient);
my $total_multi_map = $nbEndGene - $nbOriginalL0Map;
#my $nbGeneMapped = $nbOriginalL0Map - ($nbOriginal_total_multiMap_seqdif - $nbOriginal_multiMap_seqdif); # Total of gene mapped

my $messageEnd;
$messageEnd.= "\nTo resume:\n==========\n\n";
$messageEnd.= "The original file contained $nbOriginalGene genes\n\n";

$messageEnd.= "Before filtering:\nWe have $nbRecordMap mapped genes for a total of $nbOriginalL0Map maps.\n$nbOriginal_multiMap_seqdif genes have several maps on different sequences for a total of $nbOriginal_total_multiMap_seqdif maps.\n";
$messageEnd.= "/!\\ Multi map on same sequence not taken into account\n\n";

$messageEnd.= "After filtering:\nWe have $nbGeneIdUniqMap mapped genes for a total of $nbEndGene maps.(Over the $valueK % match threshold).\n";
my $nbGeneOneMap=$nbGeneIdUniqMap - $nb_multiMap_seqdif - $nb_multiMap_sameseq - $bothCase;
$messageEnd.= "We have $nbGeneOneMap genes that map to only one location.\n";
$messageEnd.= "We have $nb_multiMap_seqdif genes that map on different sequences for a total of $nb_total_multiMap_seqdif maps \n";
$messageEnd.= "We have $nb_multiMap_sameseq genes that map on different location of the same sequences for a total of $nb_total_multiMap_sameseq maps \n";
$messageEnd.= "We have $bothCase genes that map on different location of the same sequences and on different sequences for a total of $nb_total_multiMap_sameseq_bothcase maps on same sequence and $nb_total_multiMap_seqdif_bothcase maps on different sequences\n";
#report error if necessary
if( $nb_noCaseL3 ){
  $messageEnd.= "About potential problem met:\nWe found $nb_noCaseL3 cases where a l1/l2 feature mapped to true but do not have any feature l3. As we are using feature l3 to perform the analysis we have skiped them.\n";
  $messageEnd.= "When it was the case for all the l2 of one recored we have removed the record (l1): number of such case: $KrakenFakeGene.\n";
}
$messageEnd.= "\n";

#print info
if ($outfile) {
print $outReport $messageEnd;
}else{print $messageEnd;}

#############
#PLOT
#############
## -----Manage plot output file-----
# Check R is available. If not we try to load it through Module software
my ($pathPlotFile, $pathOutPlot, $ostreamPlotFile);

if ( system("R --version 1>/dev/null 2>/dev/null") == 0 ) {
	print "R is available. We can also provide a plot as result.\n";

	$ostreamPlotFile = new IO::File;
	$pathPlotFile="geneMapped.txt";
	$pathOutPlot="geneMapped_plot.pdf";
	if ($outfile) {
	  $pathPlotFile=$outfile."-geneMapped.txt";
	  $pathOutPlot=$outfile."-geneMapped_plot.pdf";
	}
	$ostreamPlotFile->open($pathPlotFile, 'w' ) or
        croak(
          sprintf( "Can not open '%s' for writing %s", $pathPlotFile, $! )
        );

	###############
	# print the value per gene in a temporary file for R plot
	foreach my $key (keys %mappedPercentPerGene){
		if ($mappedPercentPerGene{$key} > 100){
				print $ostreamPlotFile "100\n";
				warn "Warning: $key mapped value over 100%: ".$mappedPercentPerGene{$key}."%\n";
		}
		else{
	   print $ostreamPlotFile $mappedPercentPerGene{$key}."\n";
	 	}
	}

	my $messagePlot;
	if ($nbGeneIdUniqMap){
	  # Create the legend
	  my $nbOfGeneSelected = $nbGeneIdUniqMap;
	  # parse file name to remove extension
	  my ($file1,$dir1,$ext1) = fileparse($gff, qr/\.[^.]*/);
	  my $legend=$nbOfGeneSelected." genes selected from ".$file1;

	  my @listTuple=([$pathPlotFile,$legend]);
	  my $R_command=rcc_density_one_row_per_file(\@listTuple,"histogram","Percentage of gene length mapped","10","",$pathOutPlot); # create the R command
	  execute_R_command($R_command);

	  $messagePlot = "Plot done in the pdf file named $pathOutPlot\n";
	}
	else{
	  $messagePlot = "Cannot perform any plot without data.\n";
	}

	#print info
	if ($outfile) {
	  print $outReport $messagePlot;
	}
	else{print $messagePlot;}

	# Delete temporary file
	#unlink "$pathPlotFile";
}
else {
	print "R no available. We cannot perform any plot\n";
}
#END
print "We finished !! Bye Bye.\n";

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


# We have the l1 id from hash omniscient and the omniscient, and we are looking for the proper feature l2 and its subfeture l3.
# When this function is called, a new_omniscient containing only one isoform has been newly created.
# To be sure to retrieve the proper l2 from which the current l2 has been created,
# We look a the transcript_id attribute rather than the ID, because only this attribute is stayed un-modified.
#
sub compute_total_size{
  my ($hash_omniscient, $l1_original_id, $feature_l3)=@_;

		print $l1_original_id." = l1_original_id\n" if ($verbose);

		  my $l2_transcipt_id = lc($feature_l3->_tag_value('transcript_id'));
		  my $total_size=0;


      foreach my $primary_tag_key_level2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
        if ( exists_keys($hash_omniscient, ('level2',$primary_tag_key_level2, $l1_original_id) ) ){

          my $found=undef;
          foreach my $feature_level2 ( @{$hash_omniscient->{'level2'}{$primary_tag_key_level2}{$l1_original_id}}) {
            if($l2_transcipt_id eq lc($feature_level2->_tag_value('transcript_id'))){
							#print $l2_transcipt_id." = l2_transcipt_id\n" if ($verbose);
              my $l2_original_id = lc($feature_level2->_tag_value('ID'));
              foreach my $primary_tag_l3  (keys %{$hash_omniscient->{'level3'}}){ # primary_tag_l3 = cds or exon or start_codon or utr etc...
                if( lc($primary_tag_l3) eq lc( $feature_l3->primary_tag() ) ){
                  if ( exists_keys( $hash_omniscient, ('level3', $primary_tag_l3, $l2_original_id) ) ){
                    foreach my $feature_level3 ( @{$hash_omniscient->{'level3'}{$primary_tag_l3}{$l2_original_id}}){
											my $l3_transcipt_id = lc($feature_level3->_tag_value('transcript_id'));
											#print $l3_transcipt_id." = l3_transcipt_id\n" if ($verbose);
                      $total_size+=($feature_level3->end - $feature_level3->start)+1;
                    }
                  }
                }
              }
            $found =1;
            last; #No need to continue the loop if the proper l2 has been already found
            }
          }
          if(! $found){
            print "l2_transcipt_id $l2_transcipt_id not found in hash_omniscient\n";
          }
        }

  }
  if($total_size == 0){
    print "Something went wrong, total_size is 0 while we expect a positive value.\n";
  }
  return $total_size;
}

# Feature we look at are in the order exon ,cds or utr
#
sub takeOneListLevel3From1idLevel2 {
  my ($hash_omniscient, $level2_ID)=@_;

  my $refListFetaureL3=undef;
  my $refListExon=undef;
  my $refListCDS=undef;
  my $refListUTR=undef;
  my $get_one_true=undef;

  if ( exists_keys($hash_omniscient, ('level3','exon',$level2_ID) ) ){
    $refListFetaureL3 = $hash_omniscient->{'level3'}{'exon'}{$level2_ID};
    $refListExon = $hash_omniscient->{'level3'}{'exon'}{$level2_ID};

    #get if one exon mapped otherwise we have to use CDS instead
    foreach my $feature (@{$refListFetaureL3}){
      if($feature->has_tag($kraken_tag)){
        if (lc($feature->_tag_value($kraken_tag)) eq "true"){
          $get_one_true = 1;
        }
      }
    }
  }
  if ( exists_keys($hash_omniscient,('level3','cds',$level2_ID) ) and ! $get_one_true){
    $refListFetaureL3 = $hash_omniscient->{'level3'}{'cds'}{$level2_ID};
    $refListCDS = $hash_omniscient->{'level3'}{'cds'}{$level2_ID};

    #get if one cds mapped otherwise we have to use UTR instead
    foreach my $feature (@{$refListFetaureL3}){
      if($feature->has_tag($kraken_tag)){
        if (lc($feature->_tag_value($kraken_tag)) eq "true"){
          $get_one_true = 1;
        }
      }
    }
  }
  if(! $get_one_true){
    my $match=undef;
    foreach my $tag (keys %{$hash_omniscient->{'level3'}}){
      if($hash_omniscient->{'level3'}{$tag}{$level2_ID}){
        $match="yes";
        if($tag =~ "utr"){
          $refListFetaureL3 = $hash_omniscient->{'level3'}{$tag}{$level2_ID};
          $refListUTR = $hash_omniscient->{'level3'}{$tag}{$level2_ID};

          #get if one cds mapped otherwise we have to use UTR instead
          foreach my $feature (@{$refListFetaureL3}){
            if($feature->has_tag($kraken_tag)){
              if (lc($feature->_tag_value($kraken_tag)) eq "true"){
                $get_one_true = 1;
              }
            }
          }
        }
      }
    }
  }

  if(! $get_one_true){
    if ($refListExon){
      $refListFetaureL3 = $refListExon;
    }
    elsif($refListCDS){
      $refListFetaureL3 = $refListCDS;
    }
    elsif($refListUTR){
      $refListFetaureL3 = $refListUTR;
    }
    else{
      print "No feature level3 expected found for ".$level2_ID." level2 ! (Probalby an error from kraken that have added a fake l1 and consequently a fake l2. So we will remove the case.)\n";
    }
  }
  return  $refListFetaureL3;
}

sub manage_gene_label{

  my ($gene_feature, $percentMatch, $kraken_tag)=@_;
  if (! $gene_feature->has_tag($kraken_tag)){ # No kraken_mapped attribute
    label_by_value($gene_feature, $percentMatch, $kraken_tag);
  }
  else{ # kraken_mapped tag exists, check if we have to change it
    my @values = $gene_feature->get_tag_values($kraken_tag);
    my $alreadyMap = lc(shift @values) ;
    if ($alreadyMap eq "false" or $alreadyMap eq "true"){
      label_by_value($gene_feature, $percentMatch, $kraken_tag);
    }
    elsif ( ($alreadyMap ne 'full') and (( $percentMatch != 0 ) and ($alreadyMap eq 'none')) ){ # if the existing tag is full or the new tag we want to add is none ($percentMatch == 0), we skip it.
      create_or_replace_tag($gene_feature,$kraken_tag,'partial'); # add info to gene feature
    }
  }
}

sub label_by_value{
  my ($gene_feature, $percentMatch, $kraken_tag)=@_;

  if($percentMatch == 100){
    create_or_replace_tag($gene_feature,$kraken_tag,'full'); # add info to gene feature
  }
  elsif ($percentMatch != 0){
    create_or_replace_tag($gene_feature,$kraken_tag,'partial'); # add info to gene feature
  }
  else{
    create_or_replace_tag($gene_feature,$kraken_tag,'none'); # add info to gene feature
  }
}

__END__

=head1 NAME

agat_sp_kraken_assess_lift_coverage.pl

=head1 DESCRIPTION

The script takes as input gtf produced by Kraken (lift-over tool).
It will analyse the kraken_mapped attributes to calculate the mapped percentage of each mRNA.
According to a threshold (0 by default), gene with a mapping percentage over that value will be reported.
A plot nammed geneMapped_plot.pdf is performed to visualize the result.
/!\ The script handles chimeric files (i.e containg gene part mapped on the template genome and others on the de-novo one)
/!\/!\ If the file is complete (containing kraken_mapped="TRUE" and kraken_mapped="FALSE" attributes),
the script calcul the real percentage lentgh that has been mapped.
Else the calcul is only based on feature with kraken_mapped="TRUE" attributes.
So in this case the result most of time will be 100%.
/!\/!\/!\ We met rare cases where Kraken mapped a feature to several locations of the de-novo genome.
As result we could end up with mapping over > 100%. We report them as 100% mapped in the plot
and a warning is raised to allow to check thoses cases.

=head1 SYNOPSIS

    agat_sp_kraken_assess_lift_coverage --gtf infile.gtf [ -o outfile ]
    agat_sp_kraken_assess_lift_coverage --help

=head1 OPTIONS

=over 8

=item B<-gff>

Input gtf file produced by Kraken.

=item B<--threshold> or B<-t>

Gene mapping percentage over which a gene must be reported. By default the value is 0.

=item B<--verbose> or B<-v>

Verbose information.

=item B<-o> , B<--output> , B<--out> or B<--outfile>

Output GFF file.  If no output file is specified, the output will be
written to STDOUT.

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
