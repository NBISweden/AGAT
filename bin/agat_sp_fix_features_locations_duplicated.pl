#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long;
use POSIX qw(strftime);
use Pod::Usage;
use File::Basename;
use List::MoreUtils qw(uniq);
use AGAT::AGAT;

my $header = get_agat_header();
my $config;
my $threads;
my $model_to_test = undef;
my $outfile = undef;
my $ref = undef;
my $verbose = undef;
my $opt_help= 0;

my @copyARGV=@ARGV;
if ( !GetOptions(
    'c|config=s'             => \$config,
    "h|help"                 => \$opt_help,
    "f|file|gff3|gff=s"      => \$ref,
    "v|verbose!"             => \$verbose,
    "m|model=s"              => \$model_to_test,
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

if ( ! (defined($ref)) ){
    pod2usage( {
           -message => "$header\nAt least 1 parameters is mandatory:\n",
           -verbose => 0,
           -exitval => 2 } );
}

# --- Manage config ---
initialize_agat({ config_file_in => $config, input => $ref });
$CONFIG->{threads} = $threads if defined($threads);

######################
# Manage output file #
my $reportout_file;
if ($outfile) {
  my ($filename,$path,$ext) = fileparse($outfile,qr/\.[^.]*/);
  $reportout_file = $path.$filename."_report.txt" ;
}

my $gffout = prepare_gffout( $outfile );
my $reportout = prepare_fileout( $reportout_file );

# END Manage Ouput Directory / File #
#####################################

my %ListModel;
if(!($model_to_test)){
  $ListModel{1}=0;
  $ListModel{2}=0;
  $ListModel{3}=0;
  $ListModel{4}=0;
  $ListModel{5}=0;
}else{
  my @fields= split(',', $model_to_test);
  foreach my $field (@fields){
    if($field =~ m/^[12345]$/){
      $ListModel{$field}=0;
    }else{
      print "This model $field is not known. Must be an Integer !\n";exit;
    }
  }
}

my $string1 = strftime "%m/%d/%Y at %Hh%Mm%Ss", localtime;
$string1 .= "\n\nusage: $0 @copyARGV\n\n";
print $reportout $string1;
if($outfile){print $string1;}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>     MAIN     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

my $nb_gene_removed=0;

### Parse GFF input #
print ("Parse file $ref\n");
my ($omniscient, $hash_mRNAGeneLink) = slurp_gff3_file_JD({ input => $ref });

# sort by seq id
my $hash_sortBySeq = gather_and_sort_l1_location_by_seq_id_and_strand($omniscient);
my $topfeatures = get_feature_type_by_agat_value($omniscient, 'level1', 'topfeature');

#find overlap
my %checked_l1;
foreach my $seqid (sort keys %{$hash_sortBySeq}){ # loop over all the feature level1

  if( exists_keys($hash_sortBySeq,($seqid ) ) ){
    foreach my $tag (sort keys %{$hash_sortBySeq->{$seqid}}){

      #skip top features
      if(exists_keys($topfeatures,($tag))){ next; }

      foreach my $location ( sort {$a->[1].$a->[0] cmp $b->[1].$b->[0]} @{$hash_sortBySeq->{$seqid}{$tag}}){
        my $gene_feature_id = lc($location->[0]);

        if (! exists_keys($omniscient, ('level1',$tag,$gene_feature_id) ) ){ next;} #feature has been removed from previous check
        $checked_l1{$gene_feature_id}{$gene_feature_id}++; # to not check agaisnt himself
        my $gene_feature = $omniscient->{'level1'}{$tag}{$gene_feature_id};

        #################################################
        # START Take care of isoforms with duplicated location:
        print "START Take care of isoforms with duplicated locations\n" if $verbose;
        my @L2_list_to_remove = ();
        foreach my $l2_type ( sort keys %{$omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...

          if(exists_keys($omniscient,('level2', $l2_type, $gene_feature_id)) and scalar @{$omniscient->{'level2'}{$l2_type}{$gene_feature_id}} > 1){ # more than one l2 feature of that type
            #print "More than 2 mRNA let's check them\n" if $verbose;

            my %checked;
            foreach my $l2_1 (sort {$b->_tag_value('ID') cmp $a->_tag_value('ID')} @{$omniscient->{'level2'}{$l2_type}{$gene_feature_id}}){
              my $id_l2_1 = lc($l2_1->_tag_value('ID'));
              $checked{$id_l2_1}{$id_l2_1}++;

              foreach my $l2_2 (sort {$b cmp $a} @{$omniscient->{'level2'}{$l2_type}{$gene_feature_id}}){
                my $id_l2_2 = lc($l2_2->_tag_value('ID'));

                # If not itself and not already checked (A -> B is the same as B -> A), and A or B already removed and must now be skiped (skipme key)
                if( ! exists_keys(\%checked, ( $id_l2_1 , $id_l2_2 ) ) ){ #
                  $checked{$id_l2_1 }{$id_l2_2}++;
                  $checked{$id_l2_2}{$id_l2_1 }++;

                  #check their position are identical
                  if($l2_2->start().$l2_2->end() eq $l2_1->start().$l2_1->end()){

                    if(exists_keys($omniscient,('level3', 'exon', $id_l2_1))){
                      if(exists_keys($omniscient,('level3', 'exon', $id_l2_2))){
                        if(scalar @{$omniscient->{'level3'}{'exon'}{$id_l2_1}} ==  scalar @{$omniscient->{'level3'}{'exon'}{$id_l2_2}}){

                          #Check their subfeature are  identicals
                          if(featuresList_identik(\@{$omniscient->{'level3'}{'exon'}{$id_l2_1}}, \@{$omniscient->{'level3'}{'exon'}{$id_l2_2}}, $verbose )){
                            print "case1: $id_l2_2 and $id_l2_1 have same exon list\n" if ($verbose);

                            my $size_cds1 =  cds_size($omniscient, $id_l2_1);
                            my $size_cds2 =  cds_size($omniscient, $id_l2_2);
                            if($size_cds1 >= $size_cds2 ){
                              push(@L2_list_to_remove, $id_l2_2);
                              print "case1: push1\n" if $verbose;
                            }
                            elsif($size_cds1 < $size_cds2){
                              push(@L2_list_to_remove, $id_l2_1);
                              print "case1: push2\n" if $verbose;
                            }
                            elsif($size_cds1){
                              push(@L2_list_to_remove, $id_l2_2);
                              print "case1: push3\n" if $verbose;
                            }
                            else{
                              push(@L2_list_to_remove, $id_l2_1);
                              print "case1: push4\n" if $verbose;
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
        }

        if( @L2_list_to_remove ){
          if (exists($ListModel{1})){
            my @L2_list_to_remove_filtered = uniq(@L2_list_to_remove);
            $ListModel{1} += scalar @L2_list_to_remove_filtered;
            print "case1 (removing mRNA isoform identic ): ".join(",", @L2_list_to_remove_filtered)."\n";
            remove_omniscient_elements_from_level2_ID_list($omniscient, \@L2_list_to_remove_filtered);
          }
        }
        # END Take care of isoforms with duplicated location:
        #######################################################



        #######################################################
        # START Take care of other gene with duplicated location
        #
        print "START Take care of gene with duplicated locations\n" if $verbose;
        #foreach my $gene_feature_id2 (@sorted_genefeature_ids){
        foreach my $location2 (sort {$a->[1].$a->[0] cmp $b->[1].$b->[0]}  @{$hash_sortBySeq->{$seqid}{$tag}}){
          my $gene_feature_id2 = lc($location2->[0]);

          if (! exists_keys(\%checked_l1,($gene_feature_id,$gene_feature_id2 ) ) ){
            #print $gene_feature_id.":".$hash_sortBySeq->{$seqid}{'level1'}{$tag}{$gene_feature_id}[1]." $gene_feature_id2:".$hash_sortBySeq->{$seqid}{'level1'}{$tag}{$gene_feature_id2}[1]."\n";
            if ( $location2->[1] > $location->[1] ) { last; } # start of gene2 is over start of gene 1 we could stop to loop... no need to look at following genes in the list
            if (! exists_keys($omniscient, ('level1',$tag,$gene_feature_id2) ) ){ next;} #feature has been removed from previous check

            $checked_l1{$gene_feature_id }{$gene_feature_id2}++;
            $checked_l1{$gene_feature_id2 }{$gene_feature_id}++;
            my @L2_list_to_remove = ();
            my $gene_feature2 = $omniscient->{'level1'}{$tag}{$gene_feature_id2};

            #The two genes overlap
            if( ($gene_feature2->start <= $gene_feature->end() ) and ($gene_feature2->end >= $gene_feature->start) ){

              print "$gene_feature_id and $gene_feature_id2 overlap\n" if $verbose;

              # Loop over the L2 from the first gene feature
              foreach my $l2_type ( sort keys %{$omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
                if ( exists ($omniscient->{'level2'}{$l2_type}{$gene_feature_id} ) ){

                  foreach my $l2_1 (sort {$b->_tag_value('ID') cmp $a->_tag_value('ID')} @{$omniscient->{'level2'}{$l2_type}{$gene_feature_id}}){
                    my $id_l2_1 = lc($l2_1->_tag_value('ID'));

                    # Loop over the L2 from the second gene feature
                    if ( exists ($omniscient->{'level2'}{$l2_type}{$gene_feature_id2} ) ){

                      my $keep = 1;
                      foreach my $l2_2 (sort {$b->_tag_value('ID') cmp $a->_tag_value('ID')} @{$omniscient->{'level2'}{$l2_type}{$gene_feature_id2}}){
                        my $id_l2_2 = lc($l2_2->_tag_value('ID'));

                        #check their position are identical
                        if($l2_2->start().$l2_2->end() eq $l2_1->start().$l2_1->end()){

                          print "$id_l2_2  and $id_l2_1 have same start and stop\n" if $verbose;

                          if(exists_keys($omniscient,('level3', 'exon', $id_l2_1))){
                            if(exists_keys($omniscient,('level3', 'exon', $id_l2_2))){

                              my $resu_overlap = check_feature_overlap_from_l3_to_l1($omniscient, $omniscient , $gene_feature_id, $gene_feature_id2);
                              if ($resu_overlap){

                                print "$id_l2_2  and $id_l2_1 overlap at $resu_overlap\n" if $verbose;

                                #EXON identicals
                                if(featuresList_identik(\@{$omniscient->{'level3'}{'exon'}{$id_l2_1}}, \@{$omniscient->{'level3'}{'exon'}{$id_l2_2}}, $verbose )){

                                  print "$id_l2_2 and $id_l2_1 have same exon list\n" if $verbose;
                                  # NO CDS
                                  if ( ! exists_keys($omniscient, ('level3','cds',$id_l2_1)) and  ! exists_keys($omniscient, ('level3','cds',$id_l2_2) ) ) {
                                    if (exists($ListModel{2})){
                                       print "case2: $id_l2_2 and $id_l2_1 have no CDS\n" if $verbose;
                                       $ListModel{2}++;
                                       push(@L2_list_to_remove, $id_l2_2);
                                    }
                                  }
                                  else{ # WITH CDS
                                    if(featuresList_identik(\@{$omniscient->{'level3'}{'cds'}{$id_l2_1}}, \@{$omniscient->{'level3'}{'cds'}{$id_l2_2}}, $verbose ) ){
                                      if ( exists($ListModel{3}) ){
                                        print "case3: $id_l2_2 and $id_l2_1 have same CDS list\n" if $verbose;
                                        $ListModel{3}++;
                                        #identik because no CDS, we could remove one randomly

                                        my $size_cds1 =  cds_size($omniscient, $id_l2_1);
                                        my $size_cds2 =  cds_size($omniscient, $id_l2_2);
                                        if($size_cds1 >= $size_cds2 ){
                                          push(@L2_list_to_remove, $id_l2_2);
                                          print "case3: push1 $size_cds1 $size_cds2\n" if $verbose;
                                        }
                                        elsif($size_cds1 < $size_cds2){
                                          push(@L2_list_to_remove, $id_l2_1);
                                          print "case3: push2\n" if $verbose;
                                        }
                                        elsif($size_cds1){
                                          push(@L2_list_to_remove, $id_l2_2);
                                          print "case3: push3\n" if $verbose;
                                        }
                                        else{
                                          push(@L2_list_to_remove, $id_l2_1);
                                          print "case3: push4\n" if $verbose;
                                        }
                                      }
                                    }

                                    # CDS are not identic Let's reshape UTRS
                                    elsif ( exists($ListModel{4})){
                                      $ListModel{4}++;
                                      print "case4 (Exon structure identic from different genes, but CDS different, Let's reshape the UTRs to make them different.): $id_l2_1 <=> $id_l2_2\n";
                                      reshape_the_2_l2_models($omniscient, $gene_feature, $l2_1, $gene_feature2, $l2_2, $verbose, 4);
                                    }
                                  }
                                }
                                # Exon structure different inside
                                elsif ( exists($ListModel{5})) {
                                  $ListModel{5}++;
                                  print "case5 (Exons overlap but structure different (Same extremities but different internal locations) Let's reshape the UTRs to make them different.): $id_l2_1 <=> $id_l2_2\n";
                                  reshape_the_2_l2_models($omniscient, $gene_feature, $l2_1, $gene_feature2, $l2_2, $verbose, 5);
                                }
                              }
                              # CDS and Exon does not overlap
                              else{
                                print "CDS and Exon does not overlap\n" if ($verbose);

                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
              if(@L2_list_to_remove){
                print "case2 (removing mRNA identic from different genes: ".join(",", @L2_list_to_remove)."\n";
                remove_omniscient_elements_from_level2_ID_list($omniscient, \@L2_list_to_remove);
                if (! exists_keys($omniscient, ('level1',$tag,$gene_feature_id2) ) or ! exists_keys($omniscient, ('level1',$tag,$gene_feature_id) ) ){ $nb_gene_removed++;}
              }
            }
          }
        }
      }
    }
  }
}

my $string_print = "\nresults:\n";
if (exists($ListModel{1})){
$string_print .= "Case1: AGAT found ".$ListModel{1}." cases where isoforms have identical exon structures (AGAT removed duplicates by keeping the one with longest CDS).\n";
}
if (exists($ListModel{2})){
  $string_print .= "Case2: AGAT found ".$ListModel{2}." cases where l2 from different gene identifier have identical exon but no CDS at all (AGAT removed one duplicate).\n";
}
if (exists($ListModel{3})){
  $string_print .= "Case3: AGAT found ".$ListModel{3}." cases where l2 from different gene identifier have identical exon and CDS structures (AGAT removed duplicates by keeping the one with longest CDS).\n";
}
if (exists($ListModel{4})){
  $string_print .= "Case4: AGAT found ".$ListModel{4}." cases where l2 from different gene identifier have identical exon structures and different CDS structures (AGAT reshaped UTRs to modify gene locations).\n";
}
if (exists($ListModel{5})){
  $string_print .= "Case5: AGAT found ".$ListModel{5}." cases where l2 from different gene identifier overlap but have different exon structure. In that case AGAT modified the gene locations by clipping UTRs.\n";
}

if (exists_keys(\%ListModel,("noclip"))){
  foreach my $case ( sort keys %{$ListModel{"noclip"}} ){
    my $nb = $ListModel{"noclip"}{$case};
    $string_print .= "$nb Case$case problems: No UTRs available to modify the gene model safely.\n";
  }
}

$string_print .= "AGAT removed $nb_gene_removed genes because no more l2 were linked to them.\n";

print_omniscient( {omniscient => $omniscient, output => $gffout} );

print $reportout $string_print;
if($outfile){print $string_print;}

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

# shortened UTR and exon by 1 bp in one extremity
sub reshape_the_2_l2_models{
  my ($omniscient, $gene_feature, $l2_1, $gene_feature2, $l2_2, $verbose, $case)=@_;

  my $id_l2_1 = lc($l2_1->_tag_value('ID'));
  my $parent_l2_1 = lc($l2_1->_tag_value('Parent'));
  my $l2_1_strand = $l2_1->strand;

  my $id_l2_2 = lc($l2_2->_tag_value('ID'));
  my $parent_l2_2 = lc($l2_2->_tag_value('Parent'));
  my $l2_2_strand = $l2_2->strand;

  # get info about UTRs presence
  my $left_UTR1 = get_extremity_feature_l3_from_l2id($omniscient, $gene_feature, $id_l2_1, "UTR", "left");
  my $right_UTR1 = get_extremity_feature_l3_from_l2id($omniscient, $gene_feature, $id_l2_1, "UTR", "right");
  my $left_UTR2 = get_extremity_feature_l3_from_l2id($omniscient, $gene_feature2, $id_l2_2, "UTR", "left");
  my $right_UTR2 = get_extremity_feature_l3_from_l2id($omniscient, $gene_feature2, $id_l2_2, "UTR", "right");

  if ($left_UTR1){
    print "modify $id_l2_1 left\n" if $verbose;
    $left_UTR1->start($left_UTR1->start+1);
    my $left_exon = get_extremity_feature_l3_from_l2id($omniscient, $gene_feature, $id_l2_1, "exon", "left");
    $left_exon->start($left_exon->start+1);
    check_record_positions($omniscient, $parent_l2_1);
  }
  elsif ($right_UTR1){
    print "modify $id_l2_1 right\n" if $verbose;
    $right_UTR1->end($right_UTR1->end-1);
    my $right_exon = get_extremity_feature_l3_from_l2id($omniscient, $gene_feature, $id_l2_1, "exon", "right");
    $right_exon->end($right_exon->end-1);
    check_record_positions($omniscient, $parent_l2_1);
  }
  elsif ($left_UTR2){
    print "modify $id_l2_2 left\n" if $verbose;
    $left_UTR2->start($left_UTR2->start+1);
    my $left_exon = get_extremity_feature_l3_from_l2id($omniscient, $gene_feature2, $id_l2_1, "exon", "left");
    $left_exon->start($left_exon->start-1);
    check_record_positions($omniscient, $parent_l2_2);
  }
  elsif ($right_UTR2){
    print "modify $id_l2_2 right\n" if $verbose;
    $right_UTR2->end($right_UTR2->end-1);
    my $right_exon = get_extremity_feature_l3_from_l2id($omniscient, $gene_feature2, $id_l2_1, "exon", "right");
    $right_exon->end($right_exon->end-1);
    check_record_positions($omniscient, $parent_l2_2);
  }
  else{
    print "$id_l2_1 and $id_l2_2 do not have UTRs, we cannot modify one to make the features different.".
    "You might try EvidenceModeler to choose or modify the gene models automatically,".
    " or you can manually modify them.\n";
    # We might add UTR but in someway we should avoid to goes over extremities
    $ListModel{$case}--;
    $ListModel{"noclip"}{$case}++;
  }
}

sub get_extremity_feature_l3_from_l2id{
  my ($omniscient, $gene_feature, $level2_ID, $tag, $side)=@_;

  my $result = undef;
  my $feature = undef;
  if($side ne "right" and $side ne "left"){print "Error - must be right or left\n";exit 1;}

  if($tag eq "UTR"){
    foreach my $tag ( keys %{$omniscient->{'level3'}} ){
      if(lc($tag) =~ m/utr/){
        if(exists_keys($omniscient,('level3', $tag, $level2_ID))){
          foreach my $feature_level3 ( @{$omniscient->{'level3'}{$tag}{$level2_ID}}) {
            if ( $side eq "left" and ( !$result or $feature_level3->start() < $result) ){
                $result = $feature_level3->start();
                $feature = $feature_level3;
            }
            if ( $side eq "right" and ( !$result or $feature_level3->end() > $result) ){
                $result = $feature_level3->end();
                $feature = $feature_level3;
            }
          }
        }
      }
    }
  }
  else{
    if ( exists ($omniscient->{'level3'}{$tag}{$level2_ID} ) ){
      foreach my $feature_level3 ( @{$omniscient->{'level3'}{$tag}{$level2_ID}}) {

        if ( $side eq "left" and ( !$result or $feature_level3->start() < $result) ){
            $result = $feature_level3->start();
            $feature = $feature_level3;
        }
        if ( $side eq "right" and ( !$result or $feature_level3->end() > $result) ){
            $result = $feature_level3->end();
            $feature = $feature_level3;
        }
      }
    }
  }

  if ( $feature and $side eq "left" and $feature->start() != $gene_feature->start() ){
    $feature = undef;
  }
  if ( $feature and $side eq "right" and $feature->end() != $gene_feature->end() ){
    $feature = undef;
  }
  return $feature;
}


sub cds_size{
  my ($omniscient, $id_l2)=@_;
  my $size=undef;

  if ( exists_keys($omniscient, ('level3', 'cds', lc($id_l2)))){
    foreach my $l3 ( @{$omniscient->{'level3'}{'cds'}{lc($id_l2)}} ){
      $size+= $l3->end - $l3->start +1;
    }
  }
  return $size;
}

__END__

=head1 NAME

agat_sp_fix_features_locations_duplicated.pl

=head1 DESCRIPTION

The script aims to modify/remove feature with duplicated locations. Even if it
not an error by itself in a gtf/gff file, it becomes problematic when submitting
the file to ENA (after convertion).
To modify locations, AGAT modify the UTRs (when available) by shortening them by 1 bp (and consequently the Parent features and the exons accordingly)

* Case1: When isoforms have identical exon structures, AGAT removes duplicates by keeping the one with longest CDS;
* Case2: When l2 (e.g. mRNA) from different gene identifier have identical exon but no CDS at all, AGAT removes one duplicate);
* Case3: When l2 (e.g. mRNA) from different gene identifier have identical exon and CDS structures, AGAT removes duplicates by keeping the one with longest CDS);
* Case4: When l2 (e.g. mRNA) from different gene identifier have identical exon structures and different CDS structures, AGAT reshapes UTRs to modify mRNA and gene locations);
* Case5: When l2 (e.g. mRNA) from different gene identifier overlap but have different exon structure. In that case AGAT modified the gene locations by clipping UTRs;

=head1 SYNOPSIS

    agat_sp_fix_features_locations_duplicated.pl --gff infile  [-o outfile]
    agat_sp_fix_features_locations_duplicated.pl --help

=head1 OPTIONS

=over 8

=item B<-f>, B<--file>, B<--gff3> or B<--gff>

Input GTF/GFF file.

=item B<-m> or B<--model>

To select cases you want to fix. By default all are used.
To select specific cases write e.g. --model 1,4,5

Case1: When isoforms have identical exon structures AGAT removes duplicates by keeping the one with longest CDS;
Case2: When l2 (e.g. mRNA) from different gene identifier have identical exon but no CDS at all (AGAT removes one duplicate);
Case3: When l2 (e.g. mRNA) from different gene identifier have identical exon and CDS structures (AGAT removes duplicates by keeping the one with longest CDS);
Case4: When l2 (e.g. mRNA) from different gene identifier have identical exon structures and different CDS structures (AGAT reshapes UTRs to modify mRNA and gene locations);
Case5: When l2 (e.g. mRNA) from different gene identifier overlap but have different exon structure. In that case AGAT modified the gene locations by clipping UTRs;

=item B<-v> or B<--verbose>

Add verbosity.

=item B<-o>, B<--out>, B<--output> or B<--outfile>

Output file. If none given, will be display in standard output.

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

Test case for first part:
@001900F|arrow|arrow  maker gene  5082  6945  . + . ID=ACAOBTG00000034334;makerName=maker-@001900F|arrow|arrow-exonerate_est2genome-gene-0.0
@001900F|arrow|arrow  maker mRNA  5082  6945  5456  + . ID=ACAOBTM00000062562;Parent=ACAOBTG00000034334;_AED=0.22;_QI=61|1|1|1|0|0|2|575|165;_eAED=0.22;makerName=maker-@001900F|arrow|arrow-exonerate_est2genome-gene-0.0-mRNA-1;product=hypothetical protein
@001900F|arrow|arrow  maker exon  5082  5815  . + . ID=ACAOBTE00000370675;Parent=ACAOBTM00000062562;makerName=maker-@001900F|arrow|arrow-exonerate_est2genome-gene-0.0-mRNA-1:1
@001900F|arrow|arrow  maker exon  6546  6945  . + . ID=ACAOBTE00000370676;Parent=ACAOBTM00000062562;makerName=maker-@001900F|arrow|arrow-exonerate_est2genome-gene-0.0-mRNA-1:2
@001900F|arrow|arrow  maker CDS 5143  5640  . + 0 ID=ACAOBTC00000063258;Parent=ACAOBTM00000062562;makerName=maker-@001900F|arrow|arrow-exonerate_est2genome-gene-0.0-mRNA-1:cds
@001900F|arrow|arrow  maker five_prime_UTR  5082  5142  . + . ID=ACAOBTF00000063257;Parent=ACAOBTM00000062562;makerName=maker-@001900F|arrow|arrow-exonerate_est2genome-gene-0.0-mRNA-1:five_prime_utr
@001900F|arrow|arrow  maker three_prime_UTR 5641  5815  . + . ID=ACAOBTT00000063257;Parent=ACAOBTM00000062562;makerName=maker-@001900F|arrow|arrow-exonerate_est2genome-gene-0.0-mRNA-1:three_prime_utr
@001900F|arrow|arrow  maker three_prime_UTR 6546  6945  . + . ID=ACAOBTT00000063257;Parent=ACAOBTM00000062562;makerName=maker-@001900F|arrow|arrow-exonerate_est2genome-gene-0.0-mRNA-1:three_prime_utr
@001900F|arrow|arrow  maker mRNA  5082  6945  . + . ID=ACAOBTM00000062561;Parent=ACAOBTG00000034334;_AED=0.22;_QI=722|1|1|1|0|0.5|2|28|127;_eAED=0.22;makerName=augustus_masked-@001900F|arrow|arrow-processed-gene-0.22-mRNA-1;product=hypothetical protein
@001900F|arrow|arrow  maker exon  5082  5815  . + . ID=ACAOBTE00000370673;Parent=ACAOBTM00000062561;makerName=augustus_masked-@001900F|arrow|arrow-processed-gene-0.22-mRNA-1:1
@001900F|arrow|arrow  maker exon  6546  6945  . + . ID=ACAOBTE00000370674;Parent=ACAOBTM00000062561;makerName=augustus_masked-@001900F|arrow|arrow-processed-gene-0.22-mRNA-1:2
@001900F|arrow|arrow  maker CDS 5804  5815  . + 0 ID=ACAOBTC00000063257;Parent=ACAOBTM00000062561;makerName=augustus_masked-@001900F|arrow|arrow-processed-gene-0.22-mRNA-1:cds
@001900F|arrow|arrow  maker CDS 6546  6917  . + 0 ID=ACAOBTC00000063257;Parent=ACAOBTM00000062561;makerName=IDmodified-cds-30904
@001900F|arrow|arrow  maker five_prime_UTR  5082  5803  . + . ID=ACAOBTF00000063256;Parent=ACAOBTM00000062561;makerName=augustus_masked-@001900F|arrow|arrow-processed-gene-0.22-mRNA-1:five_prime_utr
@001900F|arrow|arrow  maker three_prime_UTR 6918  6945  . + . ID=ACAOBTT00000063256;Parent=ACAOBTM00000062561;makerName=augustus_masked-@001900F|arrow|arrow-processed-gene-0.22-mRNA-1:three_prime_utr


Test case for second part:
@001900F|arrow|arrow  maker gene  5082  6945  . + . ID=ACAOBTG00000034334;makerName=maker-@001900F|arrow|arrow-exonerate_est2genome-gene-0.0
@001900F|arrow|arrow  maker mRNA  5082  6945  5456  + . ID=ACAOBTM00000062562;Parent=ACAOBTG00000034334;_AED=0.22;_QI=61|1|1|1|0|0|2|575|165;_eAED=0.22;makerName=maker-@001900F|arrow|arrow-exonerate_est2genome-gene-0.0-mRNA-1;product=hypothetical protein
@001900F|arrow|arrow  maker exon  5082  5815  . + . ID=ACAOBTE00000370675;Parent=ACAOBTM00000062562;makerName=maker-@001900F|arrow|arrow-exonerate_est2genome-gene-0.0-mRNA-1:1
@001900F|arrow|arrow  maker exon  6546  6945  . + . ID=ACAOBTE00000370676;Parent=ACAOBTM00000062562;makerName=maker-@001900F|arrow|arrow-exonerate_est2genome-gene-0.0-mRNA-1:2
@001900F|arrow|arrow  maker CDS 5143  5640  . + 0 ID=ACAOBTC00000063258;Parent=ACAOBTM00000062562;makerName=maker-@001900F|arrow|arrow-exonerate_est2genome-gene-0.0-mRNA-1:cds
@001900F|arrow|arrow  maker five_prime_UTR  5082  5142  . + . ID=ACAOBTF00000063257;Parent=ACAOBTM00000062562;makerName=maker-@001900F|arrow|arrow-exonerate_est2genome-gene-0.0-mRNA-1:five_prime_utr
@001900F|arrow|arrow  maker three_prime_UTR 5641  5815  . + . ID=ACAOBTT00000063257;Parent=ACAOBTM00000062562;makerName=maker-@001900F|arrow|arrow-exonerate_est2genome-gene-0.0-mRNA-1:three_prime_utr
@001900F|arrow|arrow  maker three_prime_UTR 6546  6945  . + . ID=ACAOBTT00000063257;Parent=ACAOBTM00000062562;makerName=maker-@001900F|arrow|arrow-exonerate_est2genome-gene-0.0-mRNA-1:three_prime_utr
@001900F|arrow|arrow  maker gene  5082  6945  . + . ID=ACAOBTG00000034333;makerName=augustus_masked-@001900F|arrow|arrow-processed-gene-0.22
@001900F|arrow|arrow  maker mRNA  5082  6945  . + . ID=ACAOBTM00000062561;Parent=ACAOBTG00000034333;_AED=0.22;_QI=722|1|1|1|0|0.5|2|28|127;_eAED=0.22;makerName=augustus_masked-@001900F|arrow|arrow-processed-gene-0.22-mRNA-1;product=hypothetical protein
@001900F|arrow|arrow  maker exon  5082  5815  . + . ID=ACAOBTE00000370673;Parent=ACAOBTM00000062561;makerName=augustus_masked-@001900F|arrow|arrow-processed-gene-0.22-mRNA-1:1
@001900F|arrow|arrow  maker exon  6546  6945  . + . ID=ACAOBTE00000370674;Parent=ACAOBTM00000062561;makerName=augustus_masked-@001900F|arrow|arrow-processed-gene-0.22-mRNA-1:2
@001900F|arrow|arrow  maker CDS 5804  5815  . + 0 ID=ACAOBTC00000063257;Parent=ACAOBTM00000062561;makerName=augustus_masked-@001900F|arrow|arrow-processed-gene-0.22-mRNA-1:cds
@001900F|arrow|arrow  maker CDS 6546  6917  . + 0 ID=ACAOBTC00000063257;Parent=ACAOBTM00000062561;makerName=IDmodified-cds-30904
@001900F|arrow|arrow  maker five_prime_UTR  5082  5803  . + . ID=ACAOBTF00000063256;Parent=ACAOBTM00000062561;makerName=augustus_masked-@001900F|arrow|arrow-processed-gene-0.22-mRNA-1:five_prime_utr
@001900F|arrow|arrow  maker three_prime_UTR 6918  6945  . + . ID=ACAOBTT00000063256;Parent=ACAOBTM00000062561;makerName=augustus_masked-@001900F|arrow|arrow-processed-gene-0.22-mRNA-1:three_prime_utr
