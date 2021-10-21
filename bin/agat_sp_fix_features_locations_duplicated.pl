#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long;
use POSIX qw(strftime);
use Pod::Usage;
use File::Basename;
use List::MoreUtils qw(uniq);
use AGAT::Omniscient;

my $header = get_agat_header();
my $outfile = undef;
my $ref = undef;
my $verbose = undef;
my $opt_help= 0;

my @copyARGV=@ARGV;
if ( !GetOptions(
    "help|h" => \$opt_help,
    "f|file|gff3|gff=s" => \$ref,
    "v" => \$verbose,
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

######################
# Manage output file #
my $gffout;
my $reportout;
if ($outfile) {
  my ($filename,$path,$ext) = fileparse($outfile,qr/\.[^.]*/);
  $reportout=IO::File->new(">".$path.$filename."_report.txt" ) or croak( sprintf( "Can not open '%s' for writing %s", $filename."_report.txt", $! ));

  open(my $fh, '>', $outfile) or die "Could not open file '$filename' $!";
  $gffout= Bio::Tools::GFF->new(-fh => $fh, -gff_version => 3 );
}
else{
  $reportout = \*STDOUT or die ( sprintf( "Can not open '%s' for writing %s", "STDOUT", $! ));
  $gffout = Bio::Tools::GFF->new(-fh => \*STDOUT, -gff_version => 3);
}
# END Manage Ouput Directory / File #
#####################################


my $string1 = strftime "%m/%d/%Y at %Hh%Mm%Ss", localtime;
$string1 .= "\n\nusage: $0 @copyARGV\n\n";
print $reportout $string1;
if($outfile){print $string1;}

                  ######################
                  #        MAIN        #
                  ######################

my $nb_case1=0;
my $nb_case2aa=0;
my $nb_case2a=0;
my $nb_case2b=0;
my $nb_case3=0;
my $nb_gene_removed=0;

### Parse GFF input #
print ("Parse file $ref\n");
my ($omniscient, $hash_mRNAGeneLink) = slurp_gff3_file_JD({ input => $ref
                                                              });
print ("$ref file parsed\n");

# sort by seq id
my $hash_sortBySeq = gather_and_sort_l1_location_by_seq_id_and_strand($omniscient);
my $topfeatures = get_feature_type_by_agat_value($omniscient, 'level1', 'topfeature');

#find overlap
my %checked_l1;
foreach my $seqid (keys %{$hash_sortBySeq}){ # loop over all the feature level1

  if( exists_keys($hash_sortBySeq,($seqid ) ) ){
    foreach my $tag (keys %{$hash_sortBySeq->{$seqid}}){

      #skip top features
      if(exists_keys($topfeatures,($tag))){ next; }

      foreach my $location ( @{$hash_sortBySeq->{$seqid}{$tag}}){
        my $gene_feature_id = lc($location->[0]);

        if (! exists_keys($omniscient, ('level1',$tag,$gene_feature_id) ) ){ next;} #feature has been removed from previous check
        $checked_l1{$gene_feature_id}{$gene_feature_id}++; # to not check agaisnt himself
        my $gene_feature = $omniscient->{'level1'}{$tag}{$gene_feature_id};

        #################################################
        # START Take care of isoforms with duplicated location:
        my @L2_list_to_remove = ();
        foreach my $l2_type (keys %{$omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...

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

        if(@L2_list_to_remove){
          my @L2_list_to_remove_filtered = uniq(@L2_list_to_remove);
          $nb_case1 = $nb_case1 + scalar @L2_list_to_remove_filtered;
          print "case1 (removing mRNA isoform identic ): ".join(",", @L2_list_to_remove_filtered)."\n";
          remove_omniscient_elements_from_level2_ID_list($omniscient, \@L2_list_to_remove_filtered);
        }
        # END Take care of isoforms with duplicated location:
        #######################################################



        #######################################################
        # START Take care of other gene with duplicated location
        #
        #foreach my $gene_feature_id2 (@sorted_genefeature_ids){
        foreach my $location2 (@{$hash_sortBySeq->{$seqid}{$tag}}){
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
              foreach my $l2_type (keys %{$omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
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

                          if(exists_keys($omniscient,('level3', 'exon', $id_l2_1))){
                            if(exists_keys($omniscient,('level3', 'exon', $id_l2_2))){
                              if(scalar @{$omniscient->{'level3'}{'exon'}{$id_l2_1}} ==  scalar @{$omniscient->{'level3'}{'exon'}{$id_l2_2}}){

                                #EXON identicals
                                if(featuresList_identik(\@{$omniscient->{'level3'}{'exon'}{$id_l2_1}}, \@{$omniscient->{'level3'}{'exon'}{$id_l2_2}}, $verbose )){
                                  #EXON CDS
                                  print "case2: $id_l2_2 and $id_l2_1 have same exon list\n" if $verbose;
                                  if ( ! exists_keys($omniscient, ('level3','cds',$id_l2_1)) and  ! exists_keys($omniscient, ('level3','cds',$id_l2_2) ) ) {
                                       print "case2aa: $id_l2_2 and $id_l2_1 have no CDS\n" if $verbose;
                                       $nb_case2aa++;
                                       push(@L2_list_to_remove, $id_l2_2);
                                  }
                                  else{
                                    if(featuresList_identik(\@{$omniscient->{'level3'}{'cds'}{$id_l2_1}}, \@{$omniscient->{'level3'}{'cds'}{$id_l2_2}}, $verbose )){
                                      print "case2: $id_l2_2 and $id_l2_1 have same CDS list\n" if $verbose;
                                      $nb_case2a++;
                                      #identik because no CDS, we could remove one randomly

                                      my $size_cds1 =  cds_size($omniscient, $id_l2_1);
                                      my $size_cds2 =  cds_size($omniscient, $id_l2_2);
                                      if($size_cds1 >= $size_cds2 ){
                                        push(@L2_list_to_remove, $id_l2_2);
                                        print "case2: push1 $size_cds1 $size_cds2\n" if $verbose;
                                      }
                                      elsif($size_cds1 < $size_cds2){
                                        push(@L2_list_to_remove, $id_l2_1);
                                        print "case2: push2\n" if $verbose;
                                      }
                                      elsif($size_cds1){
                                        push(@L2_list_to_remove, $id_l2_2);
                                        print "case2: push3\n" if $verbose;
                                      }
                                      else{
                                        push(@L2_list_to_remove, $id_l2_1);
                                        print "case2: push4\n" if $verbose;
                                      }
                                    }

                                    # CDS are not identic Let's reshape UTRS
                                    else{
                                      $nb_case2b++;
                                      reshape_the_2_gene_models($omniscient, $gene_feature_id, $gene_feature_id2, $verbose);
                                      print "case2-A (Exon structure identic from different genes, but CDS different, Let's reshape the UTRs to make them different.): $id_l2_1 <=> $id_l2_2\n";
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
              }
              if(@L2_list_to_remove){
                print "case2 (removing mRNA identic from different genes: ".join(",", @L2_list_to_remove)."\n";
                remove_omniscient_elements_from_level2_ID_list($omniscient, \@L2_list_to_remove);
                if (! exists_keys($omniscient, ('level1',$tag,$gene_feature_id2) ) or ! exists_keys($omniscient, ('level1',$tag,$gene_feature_id) ) ){ $nb_gene_removed++;}
              }
            }


            # Not identik at exon level but identik a gene level. Whan arriving at this particular case, it means that the CDS do not overlap.
            # We have to shrink the UTR and reshape gene and mRNA extremities.
            if (exists_keys($omniscient, ('level1',$tag,$gene_feature_id2) ) and exists_keys($omniscient, ('level1',$tag,$gene_feature_id) ) ){

              if( ($gene_feature2->start == $gene_feature->start) and ($gene_feature2->end == $gene_feature->end) ){
                print "case3 (reshaping genes): $gene_feature_id2 and $gene_feature_id have same location \n";
                $nb_case3++;

                reshape_the_2_gene_models($omniscient, $gene_feature_id, $gene_feature_id2, $verbose);
              }
            }
          }
        }
      }
    }
  }
}

my $string_print = "\nCase1: We found $nb_case1 cases where isoforms have identical exon structures (we removed duplicates by keeping the one with longest CDS).\n";
$string_print .= "Case2: We found $nb_case2aa cases where l2 from different gene identifier have identical exon but no CDS at all (we removed one duplicate).\n";
$string_print .= "Case3: We found $nb_case2a cases where l2 from different gene identifier have identical exon and CDS structures (we removed duplicates by keeping the one with longest CDS).\n";
$string_print .= "Case4: We found $nb_case2b cases where l2 from different gene identifier have identical exon structures and different CDS structures (we reshaped UTRs to modify gene locations).\n";
$string_print .= "Case5: We found $nb_case3 cases where 2 genes have same location while their exon/CDS locations are differents. In that case we modified the gene locations by clipping UTRs.\n";
$string_print .= "We removed $nb_gene_removed genes because no more l2 were linked to them.\n";
print_omniscient($omniscient, $gffout); #print gene modified

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


sub reshape_the_2_gene_models{
  my ($omniscient, $gene_feature_id, $gene_feature_id2, $verbose)=@_;

  my $extrem_cds_start = get_extrem_cds_start($omniscient, $gene_feature_id);
  my $extrem_cds_start2 = get_extrem_cds_start($omniscient, $gene_feature_id2);
  my $extrem_cds_end = get_extrem_cds_end($omniscient, $gene_feature_id);
  my $extrem_cds_end2 = get_extrem_cds_end($omniscient, $gene_feature_id2);

  if($extrem_cds_start < $extrem_cds_start2){
    print "remove after $gene_feature_id and  before $gene_feature_id2\n" if $verbose;
    #take care of gene_feature_id
    # remodelate exon list
    remodelate_exon_list_right($omniscient, $gene_feature_id, $extrem_cds_end);
    remodelate_exon_list_left($omniscient, $gene_feature_id2, $extrem_cds_end); #+1 to avoid creating overlaping feature
  }
  else{
    print "remove before $gene_feature_id and after $gene_feature_id2\n" if $verbose;
    remodelate_exon_list_right($omniscient, $gene_feature_id2, $extrem_cds_end2);
    remodelate_exon_list_left($omniscient, $gene_feature_id, $extrem_cds_end2);
  }
  handle_l3_features($omniscient, $gene_feature_id2);
  check_record_positions($omniscient, $gene_feature_id2);
  handle_l3_features($omniscient, $gene_feature_id);
  check_record_positions($omniscient, $gene_feature_id);
}

# We will remodelate the l3 features
sub handle_l3_features{
  my ($omniscient, $id_l1)=@_;

  #  Remove all level3 feature execept exon
  my @tag_list=('exon');
  my $l2_id_list= get_all_id_l2($omniscient, $id_l1);
  my %hash_cds_positions; # keep track of start - stop
  foreach my $l2_id_x (@$l2_id_list){
      my ($cds_start, $cds_end) = get_cds_positions($omniscient, $id_l1, $l2_id_x);
      $hash_cds_positions{$l2_id_x} = [$cds_start, $cds_end];
  }
  #remove all l3 feature except exon
  remove_tuple_from_omniscient($l2_id_list, $omniscient, 'level3', 'false', \@tag_list);
  foreach my $l2_id_x (@$l2_id_list){
    my $cds_start = $hash_cds_positions{$l2_id_x}[0];
    my $cds_end = $hash_cds_positions{$l2_id_x}[1];

    my ($utr5_list, $cds_list, $utr3) = modelate_utr_and_cds_features_from_exon_features_and_cds_start_stop($omniscient->{'level3'}{'exon'}{$l2_id_x}, $cds_start, $cds_end);
    my @level1_list;
    my @level2_list;
    my @level3_list=(@$cds_list, @$utr5_list, @$utr3);
    append_omniscient($omniscient, \@level1_list, \@level2_list, \@level3_list);
  }

}

sub get_all_id_l2{
  my ($omniscient, $id_l1)=@_;
  my @result;

      foreach my $l2_type (keys %{$omniscient->{'level2'}}){ # primary_tag_l2 = mrna or mirna or ncrna or trna etc...
        if ( exists ($omniscient->{'level2'}{$l2_type}{$id_l1} ) ){
          foreach my $feature_level2 ( @{$omniscient->{'level2'}{$l2_type}{$id_l1}}) {
            my $level2_ID  = lc($feature_level2->_tag_value('ID'));
            push( @result, $level2_ID);
          }
        }
      }
  return \@result;
}


sub remodelate_exon_list_right{
  my ($omniscient, $id_l1, $limit)=@_;

      foreach my $l2_type (keys %{$omniscient->{'level2'}}){ # primary_tag_l2 = mrna or mirna or ncrna or trna etc...
        if ( exists ($omniscient->{'level2'}{$l2_type}{$id_l1} ) ){
          foreach my $feature_level2 ( @{$omniscient->{'level2'}{$l2_type}{$id_l1}}) {
            my $level2_ID  = lc($feature_level2->_tag_value('ID'));

            if ( exists ($omniscient->{'level3'}{'exon'}{$level2_ID} ) ){
              my $mustModifyList=undef;
              my @listok;
              foreach my $feature_level3 ( @{$omniscient->{'level3'}{'exon'}{$level2_ID}}) {
                if ($feature_level3->start() >= $limit){
                  $mustModifyList="yes";
                }
                elsif ($feature_level3->end() > $limit){
                  $feature_level3->end($limit);
                  push(@listok, $feature_level3);
                }
                else{
                  push(@listok, $feature_level3);
                }
              }
              if($mustModifyList){ # at least one feature has been removed from list. Save the new list
                @{$omniscient->{'level3'}{'exon'}{$level2_ID} } = @listok;
              }
            }
          }
        }
      }
}

sub remodelate_exon_list_left{
  my ($omniscient, $id_l1, $limit)=@_;

  foreach my $l2_type (keys %{$omniscient->{'level2'}}){ # primary_tag_l2 = mrna or mirna or ncrna or trna etc...
    if ( exists ($omniscient->{'level2'}{$l2_type}{$id_l1} ) ){
      foreach my $feature_level2 ( @{$omniscient->{'level2'}{$l2_type}{$id_l1}}) {
        my $level2_ID  = lc($feature_level2->_tag_value('ID'));

        if ( exists ($omniscient->{'level3'}{'exon'}{$level2_ID} ) ){
          my $mustModifyList=undef;
          my @listok;
          foreach my $feature_level3 ( @{$omniscient->{'level3'}{'exon'}{$level2_ID}}) {
            if ($feature_level3->end() <= $limit){
              $mustModifyList="yes";
            }
            elsif ($feature_level3->start() < $limit){
              $feature_level3->start($limit);
               push(@listok, $feature_level3);
            }
            else{
               push(@listok, $feature_level3);
            }
          }
          if($mustModifyList){ # at least one feature has been removed from list. Save the new list
            @{$omniscient->{'level3'}{'exon'}{$level2_ID} } = @listok;
          }
        }
      }
    }
  }
}

sub get_cds_positions{
  my ($omniscient, $id_l1, $id_l2)=@_;
  my $start=0;
  my $end=0;

  foreach my $l2_type (keys %{$omniscient->{'level2'}}){ # primary_tag_l2 = mrna or mirna or ncrna or trna etc...
    if ( exists_keys ($omniscient, ('level2', $l2_type, $id_l1) ) ){

      if ( exists_keys ($omniscient, ('level3', 'cds', $id_l2 ) ) ){
        my  @sorted_cds =  sort { $a->start() <=>  $b->start() } @{$omniscient->{'level3'}{'cds'}{$id_l2}};
        $start = @{$omniscient->{'level3'}{'cds'}{$id_l2}}[0]->start;
        $end = @{$omniscient->{'level3'}{'cds'}{$id_l2}}[$#{$omniscient->{'level3'}{'cds'}{$id_l2}}]->end;

      }
      else{
        print "WARNING $id_l2 do not exists\n";
      }
    }
  }
  return $start,$end;
}

sub get_extrem_cds_start{
  my ($omniscient, $id_l1)=@_;
  my $result=100000000000;

      foreach my $l2_type (keys %{$omniscient->{'level2'}}){ # primary_tag_l2 = mrna or mirna or ncrna or trna etc...
        if ( exists ($omniscient->{'level2'}{$l2_type}{$id_l1} ) ){
          foreach my $feature_level2 ( @{$omniscient->{'level2'}{$l2_type}{$id_l1}}) {
            my $level2_ID  = lc($feature_level2->_tag_value('ID'));
            if ( exists ($omniscient->{'level3'}{'cds'}{$level2_ID} ) ){
              foreach my $feature_level3 ( @{$omniscient->{'level3'}{'cds'}{$level2_ID}}) {
                if ($feature_level3->start() < $result){
                  $result = $feature_level3->start();
                }
              }
            }
          }
        }
      }
  return $result;
}

sub get_extrem_cds_end{
  my ($omniscient, $id_l1)=@_;
  my $result=0;

    foreach my $l2_type (keys %{$omniscient->{'level2'}}){ # primary_tag_l2 = mrna or mirna or ncrna or trna etc...
      if ( exists ($omniscient->{'level2'}{$l2_type}{$id_l1} ) ){
        foreach my $feature_level2 ( @{$omniscient->{'level2'}{$l2_type}{$id_l1}}) {
          my $level2_ID  = lc($feature_level2->_tag_value('ID'));

          if ( exists ($omniscient->{'level3'}{'cds'}{$level2_ID} ) ){
            foreach my $feature_level3 ( @{$omniscient->{'level3'}{'cds'}{$level2_ID}}) {
              if ($feature_level3->end() > $result){
                $result = $feature_level3->end();
              }
            }
          }
        }
      }
    }

  return $result;
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
the file to ena (after convertion).
* Case1: When isoforms have identical exon structures we remove duplicates by keeping the one with longest CDS;
* Case2: When l2 (e.g. mRNA) from different gene identifier have identical exon but no CDS at all (we removed one duplicate);
* Case3: When l2 (e.g. mRNA) from different gene identifier have identical exon and CDS structures (we removed duplicates by keeping the one with longest CDS);
* Case4: When l2 (e.g. mRNA) from different gene identifier have identical exon structures and different CDS structures (we reshape UTRs to modify mRNA and gene locations);
* Case5: When 2 genes have same locations while their exon/CDS locations are differents. In that case we modified the gene locations by clipping UTRs;

=head1 SYNOPSIS

    agat_sp_fix_features_locations_duplicated.pl --gff infile  [-o outfile]
    agat_sp_fix_features_locations_duplicated.pl --help

=head1 OPTIONS

=over 8

=item B<-f>, B<--file>, B<--gff3> or B<--gff>

Input GTF/GFF file.

=item B<-o>, B<--out>, B<--output> or B<--outfile>

Output file. If none given, will be display in standard output.

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
