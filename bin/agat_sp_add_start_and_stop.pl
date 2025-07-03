#!/usr/bin/env perl

use strict;
use warnings;
use POSIX qw(strftime);
use List::MoreUtils  qw(natatime);
use Sort::Naturally;
use Carp;
use Getopt::Long;
use Pod::Usage;
use Bio::DB::Fasta;
use Bio::Tools::CodonTable;
use Clone 'clone';
use AGAT::AGAT;

my $header = get_agat_header();
my $config;
my $cpu;
my $start_id = 1;
my $stop_id = 1;

my $opt_file=undef;
my $file_fasta=undef;
my $codon_table_id=1;
my $opt_output=undef;
my $opt_extend=undef;
my $opt_no_iupac=undef;
my $verbose=undef;
my $opt_help = 0;

my @copyARGV=@ARGV;
if ( !GetOptions( 'i|g|gff=s'        => \$opt_file,
                  "fasta|fa|f=s"     => \$file_fasta,
                  "table|codon|ct=i" => \$codon_table_id,
                  'o|out|output=s'   => \$opt_output,
                  'e|extend!'        => \$opt_extend,
                  'ni|na!'           => \$opt_no_iupac,
                  'v|verbose!'       => \$verbose,
                  'c|config=s'       => \$config,
                  'thread|threads|cpu|cpus|core|cores|job|jobs=i' => \$cpu,
                  'h|help!'          => \$opt_help ) )
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

if(! $opt_file or ! $file_fasta ) {
    pod2usage( {
           -message => "$header\nMust specify at least 2 parameters:\nA gff file (--gff) and a fasta file (--fasta) \n",
           -verbose => 0,
           -exitval => 1 } );
}

# --- Manage config ---
initialize_agat({ config_file_in => $config, input => $opt_file });
$CONFIG->{cpu} = $cpu if defined($cpu);

# #######################
# # START Manage Option #
# #######################
my $gffout = prepare_gffout( $opt_output );

$codon_table_id = get_proper_codon_table($codon_table_id);

my $codon_table = Bio::Tools::CodonTable->new( -id => $codon_table_id, -no_iupac => 0);
# #####################################
# # END Manage OPTION
# #####################################

                                                      #######################
                                                      #        MAIN         #
#                     >>>>>>>>>>>>>>>>>>>>>>>>>       #######################       <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#PART 1
###################################
# Read input gff3 files one by one and save value in hash of list


######################
### Parse GFF input #
my ($hash_omniscient) = slurp_gff3_file_JD({ input => $opt_file });
### END Parse GFF input #
#########################

####################
# index the genome #
my $db = Bio::DB::Fasta->new($file_fasta);
print ("Fasta file parsed\n");

my $counter_start_missing = 0;
my $counter_start_added = 0;
my $counter_end_missing = 0;
my $counter_end_added = 0;

######################
### Parse GFF input #
# get nb of each feature in omniscient;
foreach my $tag_l2 (sort keys %{$hash_omniscient->{'level2'}}){
  foreach my $id_l1 (sort { ncmp ($a, $b) } keys %{$hash_omniscient->{'level2'}{$tag_l2}}){
    foreach my $feature_l2 ( sort { ncmp ($a->start."|".$a->end.$a->_tag_value('ID'), $b->start."|".$b->end.$b->_tag_value('ID') ) } @{$hash_omniscient->{'level2'}{$tag_l2}{$id_l1}} ){

      # get level2 id
      my $id_level2 = lc($feature_l2->_tag_value('ID'));
      my $strand="+";
      if ($feature_l2->strand == -1 or $feature_l2->strand eq "-"){
        $strand="-";
      }
      print "feature strand = $strand\n" if ($verbose); 
      my $seq_id = $feature_l2->seq_id();
      print "sequence length ".$db->length($seq_id)."\n" if ($verbose);
      
      ##############################
      #If it's a mRNA = have CDS. #
      if ( exists ($hash_omniscient->{'level3'}{'cds'}{$id_level2} ) ){

        ##############
        # Manage CDS #
        my @cds_feature_list = sort {$a->start <=> $b->start} @{$hash_omniscient->{'level3'}{'cds'}{$id_level2}}; # be sure that list is sorted
        my $cds_obj = create_cds_object(\@cds_feature_list);


        #-------------------------
        #       START CASE
        #-------------------------
        print "---START CODON TEST---"."\n" if ($verbose);
        if ( exists ($hash_omniscient->{'level3'}{'start_codon'}{$id_level2} ) ){
          print "start_codon already exists for $id_level2\n" if ($verbose);
        }
        else{
          # ----- Find the start codon -----
          my $extension=0;
          my $start_codon = undef;
          if ( !$start_codon ){
            print " Try find a start codon in the CDS (GFF and GTF case) \n" if ($verbose);
            $start_codon = next_codon_is_start(\@cds_feature_list, -3);
          } 
          if ( $opt_extend and !$start_codon ){
            print " Try to extend the sequence to find a start codon further...\n" if ($verbose);
            $extension += 3;  
            # check end of seq
            my $out=undef;
            $out = is_out_of_seq_start(\@cds_feature_list, $extension + 3);# check if next codon is out of seq

            while (! $start_codon and ! $out){
              $start_codon = next_codon_is_start(\@cds_feature_list, $extension);  
              $extension += 3;
              $out = is_out_of_seq_start(\@cds_feature_list, $extension + 3); # check if next codon is out of seq
            }
          }

          # ----- Create the start codon -----
          if ($start_codon) {
            # Modify CDS to reflect result, because CDS will be used as template to create the start_codon 
            if($strand eq "+"){
              $cds_feature_list[0]->start( $cds_feature_list[0]->start() - $extension);
            } else {
              $cds_feature_list[-1]->end( $cds_feature_list[-1]->end() + $extension);
            }

            # Count number of start codon created
            $counter_start_added++;

            # create start feature
            my $ID='start_added-'.$start_id;
            $start_id++;
            my $start_feature = clean_clone( { omniscient => $hash_omniscient,
                                              feature => $cds_feature_list[0],
                                              new_primary_tag => "start_codon",
                                              new_parent => $feature_l2->_tag_value('ID'),
                                              new_id => $ID
                                              } );
            
            if($strand eq "+"){
              #set start position of the start codon
              $start_feature->start($cds_feature_list[0]->start());

              #set stop position of the start codon
              my $step=3;
              my $cpt=0;
              my $size = $cds_feature_list[$cpt]->end()-$cds_feature_list[$cpt]->start()+1;
              while($size < 3){

                my $start_feature_new = clone( $start_feature );
                $start_feature_new->end($cds_feature_list[$cpt]->start()+$size-1);
                my $ID='start_added-'.$start_id;
                $start_id++;
                create_or_replace_tag($start_feature_new,'ID', $ID); #modify ID to replace by parent val
                push @{$hash_omniscient->{'level3'}{'start_codon'}{$id_level2}}, $start_feature_new;

                $cpt++;
                $step-=$size;
                $start_feature->start($cds_feature_list[$cpt]->start());
                $size += $size + $cds_feature_list[$cpt]->end()-$cds_feature_list[$cpt]->start()+1;
              }
              $start_feature->end($cds_feature_list[$cpt]->start()+$step-1);
            }
            else{
              #set start position of the start codon
              $start_feature->end($cds_feature_list[$#cds_feature_list]->end());

              #set stop position of the start codon
              my $step=3;
              my $cpt=$#cds_feature_list;
              my $size=$cds_feature_list[$cpt]->end()-$cds_feature_list[$cpt]->start()+1;
              while($size < 3){

                my $start_feature_new = clone( $start_feature );
                $start_feature_new->start($cds_feature_list[$cpt]->end()-$size+1);
                my $ID='start_added-'.$start_id;
                $start_id++;
                create_or_replace_tag($start_feature_new,'ID', $ID); #modify ID to replace by parent val
                push @{$hash_omniscient->{'level3'}{'start_codon'}{$id_level2}}, $start_feature_new;

                $cpt--;
                $step-=$size;
                $start_feature->end($cds_feature_list[$cpt]->end());
                $size += $size + $cds_feature_list[$cpt]->end()-$cds_feature_list[$cpt]->start()+1;
              }
              $start_feature->start($cds_feature_list[$cpt]->end()-$step+1);
            }
            push @{$hash_omniscient->{'level3'}{'start_codon'}{$id_level2}}, $start_feature;
          }
          else{
            $counter_start_missing++;
          }
        }

        #-------------------------
        #       STOP CASE
        #-------------------------
        print "---STOP CODON TEST---"."\n" if ($verbose);
        if ( exists ($hash_omniscient->{'level3'}{'stop_codon'}{$id_level2} ) ){
          print "stop_codon already exists for $id_level2\n" if ($verbose);
        }
        else{ 
          # ----- Find a stop codon -----
          # Need to try last codon for GFF and codon after CDS in case of GTF
          my $extension = 0;
          my $terminal_codon = undef;
          if ( !$terminal_codon ){
            print " Try find a stop codon in the CDS (GFF case) \n" if ($verbose);
            $terminal_codon = next_codon_is_ter(\@cds_feature_list, -3);
          } 
          if ( !$terminal_codon ){
            print " Try find a stop codon next codon out of the CDS (GTF case) \n" if ($verbose);
            $terminal_codon = next_codon_is_ter(\@cds_feature_list, 0);

            if($strand eq "+"){
              $cds_feature_list[-1]->end( $cds_feature_list[-1]->end() + 3);
            } else {
              $cds_feature_list[0]->start( $cds_feature_list[0]->start() - 3);
            }
          } # with extend option 
          if ($opt_extend and !$terminal_codon){
            print " Try to extend the sequence to find a stop codon further...\n" if ($verbose);
            $extension += 3;  
            # check end of seq
            my $out=undef;
            $out = is_out_of_seq_stop(\@cds_feature_list, $extension + 3);# check if next codon is out of seq

            while (! $terminal_codon and ! $out){
              $terminal_codon = next_codon_is_ter(\@cds_feature_list, $extension);      
              $extension += 3;
              $out = is_out_of_seq_stop(\@cds_feature_list, $extension + 3); # check if next codon is out of seq
            }
          }

          # ----- Create the stop codon -----
          if ( $terminal_codon ){
            # Modify CDS to reflect result, because CDS will be used as template to create the stop_codon 
            if($strand eq "+"){
              $cds_feature_list[-1]->end( $cds_feature_list[-1]->end() + $extension);
            } else {
              $cds_feature_list[0]->start( $cds_feature_list[0]->start() - $extension);
            }

            # Count number of start codon created
            $counter_end_added++;

            # create stop feature
            my $ID='stop_added-'.$stop_id;
            $stop_id++;
            my $stop_feature = clean_clone( { omniscient => $hash_omniscient,
                                              feature => $cds_feature_list[0],
                                              new_primary_tag => "stop_codon",
                                              new_parent => $feature_l2->_tag_value('ID'),
                                              new_id => $ID
                                              } );

            if($strand eq "+"){

              # set start position of the stop codon
              $stop_feature->end($cds_feature_list[$#cds_feature_list]->end());

              #set stop position of the stop codon
              my $step=3;
              my $cpt=$#cds_feature_list;
              my $size=$cds_feature_list[$cpt]->end()-$cds_feature_list[$cpt]->start()+1;
              while($size < 3){

                my $stop_feature_new = clone( $stop_feature );
                $stop_feature_new->start($cds_feature_list[$cpt]->end()-$size+1);
                my $ID='start_added-'.$start_id;
                $start_id++;
                create_or_replace_tag($stop_feature_new,'ID', $ID); #modify ID to replace by parent val
                push @{$hash_omniscient->{'level3'}{'stop_codon'}{$id_level2}}, $stop_feature_new;

                $cpt--;
                $step-=$size;
                $stop_feature->end($cds_feature_list[$cpt]->end());
                $size += $size + $cds_feature_list[$cpt]->end()-$cds_feature_list[$cpt]->start()+1;
              }
              #print $cds_feature_list[$cpt]->end()."\n";
              $stop_feature->start($cds_feature_list[$cpt]->end()-$step+1);
            }
            else{
              #set start position of the stop codon
              $stop_feature->start($cds_feature_list[0]->start());

              #set stop position of the stop codon
              my $step=3;
              my $cpt=0;
              my $size = $cds_feature_list[$cpt]->end()-$cds_feature_list[$cpt]->start()+1;
              while($size < 3){

                my $stop_feature_new = clone( $stop_feature );
                $stop_feature_new->end($cds_feature_list[$cpt]->start()+$size-1);
                my $ID='start_added-'.$start_id;
                $start_id++;
                create_or_replace_tag($stop_feature_new,'ID', $ID); #modify ID to replace by parent val
                push @{$hash_omniscient->{'level3'}{'stop_codon'}{$id_level2}}, $stop_feature_new;

                $cpt++;
                $step-=$size;
                $stop_feature->start($cds_feature_list[$cpt]->start());
                $size += $size + $cds_feature_list[$cpt]->end()-$cds_feature_list[$cpt]->start()+1;
              }
              $stop_feature->end($cds_feature_list[$cpt]->start()+$step-1);
            }
            push @{$hash_omniscient->{'level3'}{'stop_codon'}{$id_level2}}, $stop_feature;
          }
          else{
            $counter_end_missing++;
          }
        }
      }
    }
  }
}

# case we need to check start stop of the features not CDS
if ($opt_extend){
  my ($new_hash_omniscient) = slurp_gff3_file_JD({ input => $hash_omniscient,
                                                   config => $config });
  print_omniscient( {omniscient => $new_hash_omniscient, output => $gffout} );
} else {
  print_omniscient( {omniscient => $hash_omniscient, output => $gffout} );
}

print "$counter_start_added start codon added and $counter_start_missing CDS do not start by a start codon\n";
print "$counter_end_added stop codon added and $counter_end_missing CDS do not end by a stop codon \n";
print "bye bye\n";

      #########################
      ######### END ###########
      #########################

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

sub create_cds_object{

  my ($feature_list) = @_;

  my $seq = "";
  foreach my $feature (@$feature_list) {
    my $start=$feature->start();
    my $end=$feature->end();
    my $seqid=$feature->seq_id();
    $seq .= $db->seq( $seqid, $start, $end );
  }
  print "sequence: $seq\n" if ($verbose);
  my $debut = $db->seq( $feature_list->[0]->seq_id(), $feature_list->[0]->start-1, $feature_list->[0]->start -3 );
  my $fin = $db->seq( $feature_list->[0]->seq_id(), $feature_list->[-1]->end+1, $feature_list->[-1]->end+ 3 );
  print "sequence_extended: $debut$seq$fin\n" if ($verbose);

  #create the cds object
  my $cds_obj = Bio::Seq->new(-seq => $seq, -alphabet => 'dna' );
  
  #Reverse the object depending on strand
  my $strand="+";
  if ($feature_list->[0]->strand == -1 or $feature_list->[0]->strand eq "-"){
      $cds_obj = $cds_obj->revcom();
      $strand = "-";
      print "feature on minus strand\n" if ($verbose);
      print "sequence: ".$cds_obj->seq."\n" if ($verbose);
  }
  
  return $cds_obj;
}

sub next_codon_is_start{
  my ($feature_list, $more) = @_;
  if(! $more){
    $more=0;
  }
  print "next_codon_is_start test: \n" if ($verbose);
  my $cds_obj;
  my $seqid=$feature_list->[0]->seq_id();

    # Minus strand
  if ($feature_list->[0]->strand == -1 or $feature_list->[0]->strand eq "-"){
      my $end=$feature_list->[-1]->end();
      my $seq = $db->seq( $seqid,$end+1+$more, $end+3+$more);
      $cds_obj = Bio::Seq->new(-seq => $seq, -alphabet => 'dna' );
      $cds_obj = $cds_obj->revcom();
      print "  Minus strand - most right side: ".($end+3+$more)."\n" if ($verbose);
  }
  else{ # Plus strand
      my $start=$feature_list->[0]->start();
      my $seq = $db->seq( $seqid, $start-3-$more,  $start-1-$more);
      $cds_obj = Bio::Seq->new(-seq => $seq, -alphabet => 'dna' );
      print "  Plus strand - most right side: ".($start-3-$more)."\n" if ($verbose);
  }

  my $codon = $cds_obj->seq ;
  print "  codon tested is = $codon \n" if ($verbose);

  if ( !is_ambiguous_codon($codon) and $codon_table->is_start_codon( $codon )){
    print "  It is considered as a start codon!\n" if ($verbose);;
    return 1;
  } else{
    return 0;
  }
}

sub next_codon_is_ter{

  my ($feature_list, $more) = @_;
  if(! $more){
    $more=0;
  }
  print "next_codon_is_ter test: \n" if ($verbose);
  my $cds_obj;
  my $seqid=$feature_list->[0]->seq_id();

  # Minus strand
  if ($feature_list->[0]->strand == -1 or $feature_list->[0]->strand eq "-"){
      my $start=$feature_list->[0]->start();
      my $seq = $db->seq( $seqid,$start-3-$more, $start-1-$more);
      $cds_obj = Bio::Seq->new(-seq => $seq, -alphabet => 'dna' );
      $cds_obj = $cds_obj->revcom();
      print "  Minus strand - most right side: ".($start-3-$more)."\n" if ($verbose);
  }
  else{ # Plus strand
      my $end=$feature_list->[-1]->end();
      my $seq = $db->seq( $seqid, $end+1+$more,  $end+3+$more);
      $cds_obj = Bio::Seq->new(-seq => $seq, -alphabet => 'dna' );
      print "  Plus strand - most right side: ".($end+3+$more)."\n" if ($verbose);
  }
  
  my $codon = $cds_obj->seq ;
  print "  codon tested is = $codon \n" if ($verbose);

  if ( !is_ambiguous_codon($codon) and $codon_table->is_ter_codon( $codon )){
    print "  It is considered as a stop codon!\n" if ($verbose);
    return 1;
  } else{
    return 0;
  }

}

sub is_out_of_seq_start{

  my ($feature_list, $more) = @_;
  my $out = undef;

  my $strand="+";
  if ($feature_list->[0]->strand == -1 or $feature_list->[0]->strand eq "-"){
    $strand="-";
  }
  my $length_seqid = $db->length($feature_list->[0]->seq_id())."\n";
  
  if($strand eq "+"){
    if( $feature_list->[0]->start() - $more < 1 ){
      $out = 1;
      print "is_out_of_seq_start!! Plus strand - Most left out of seq: ".($feature_list->[0]->start() - $more)."\n" if ($verbose);
    }
  }
  else{
    if( $feature_list->[-1]->end() + $more > $length_seqid ){
      $out = 1;
      print "is_out_of_seq_start!! Minus strand - Most right out of seq: ".($feature_list->[-1]->end() + $more)."\n" if ($verbose);
    }
  }
  return $out;
}

sub is_out_of_seq_stop{

  my ($feature_list, $more) = @_;
  my $out = undef;

  my $strand="+";
  if ($feature_list->[0]->strand == -1 or $feature_list->[0]->strand eq "-"){
    $strand="-";
  }
  my $length_seqid = $db->length($feature_list->[0]->seq_id())."\n";
  
  if($strand eq "+"){
    if( $feature_list->[-1]->end() + $more > $length_seqid ){
      $out = 1;
      print "is_out_of_seq_stop!! Plus strand - Most right out of seq: ".($feature_list->[-1]->end() + $more)."\n" if ($verbose);
    }
  }
  else{
    if( $feature_list->[0]->start() - $more < 1 ){
      $out = 1;
      print "is_out_of_seq_stop!! Minus strand - Most right out of seq: ".($feature_list->[0]->start() - $more)."\n" if ($verbose);
    }
  }
  return $out;
}

# Test done only if option activated
sub is_ambiguous_codon{

  my ($codon) = @_;

  if ($opt_no_iupac){
    if ( $codon !~ /[ATGC]{3}/) {
      print "$codon is an ambiguous codon we skip it because the no_iupac option is activated!\n" if ($verbose);
      return 1;
    } 
  }
  return 0;
}

__END__
EXAMPLE NORMAL
##gff-version 3
Pcoprophilum_scaf_9 . contig  1 1302582 . . . ID=Pcoprophilum_scaf_9;Name=Pcoprophilum_scaf_9
Pcoprophilum_scaf_9 maker gene  189352  192747  . + . ID=genemark-Pcoprophilum_scaf_9-processed-gene-2.0;Name=genemark-Pcoprophilum_scaf_9-processed-gene-2.0
Pcoprophilum_scaf_9 maker mRNA  189352  192747  . + . ID=genemark-Pcoprophilum_scaf_9-processed-gene-2.0-mRNA-1;Parent=genemark-Pcoprophilum_scaf_9-processed-gene-2.0;Name=genemark-Pcoprophilum_scaf_9-processed-gene-2.0-mRNA-1;_AED=0.00;_eAED=0.00;_QI=398|1|1|1|0.5|0.33|3|343|825
Pcoprophilum_scaf_9 maker exon  189352  189520  . + . ID=genemark-Pcoprophilum_scaf_9-processed-gene-2.0-mRNA-1:exon:96;Parent=genemark-Pcoprophilum_scaf_9-processed-gene-2.0-mRNA-1
Pcoprophilum_scaf_9 maker exon  189643  189922  . + . ID=genemark-Pcoprophilum_scaf_9-processed-gene-2.0-mRNA-1:exon:97;Parent=genemark-Pcoprophilum_scaf_9-processed-gene-2.0-mRNA-1
Pcoprophilum_scaf_9 maker exon  189978  192747  . + . ID=genemark-Pcoprophilum_scaf_9-processed-gene-2.0-mRNA-1:exon:98;Parent=genemark-Pcoprophilum_scaf_9-processed-gene-2.0-mRNA-1
Pcoprophilum_scaf_9 maker five_prime_UTR  189352  189520  . + . ID=genemark-Pcoprophilum_scaf_9-processed-gene-2.0-mRNA-1:five_prime_utr;Parent=genemark-Pcoprophilum_scaf_9-processed-gene-2.0-mRNA-1
Pcoprophilum_scaf_9 maker five_prime_UTR  189643  189871  . + . ID=genemark-Pcoprophilum_scaf_9-processed-gene-2.0-mRNA-1:five_prime_utr;Parent=genemark-Pcoprophilum_scaf_9-processed-gene-2.0-mRNA-1
Pcoprophilum_scaf_9 maker CDS 189872  189922  . + 0 ID=genemark-Pcoprophilum_scaf_9-processed-gene-2.0-mRNA-1:cds;Parent=genemark-Pcoprophilum_scaf_9-processed-gene-2.0-mRNA-1
Pcoprophilum_scaf_9 maker CDS 189978  192404  . + 0 ID=genemark-Pcoprophilum_scaf_9-processed-gene-2.0-mRNA-1:cds;Parent=genemark-Pcoprophilum_scaf_9-processed-gene-2.0-mRNA-1
Pcoprophilum_scaf_9 maker three_prime_UTR 192405  192747  . + . ID=genemark-Pcoprophilum_scaf_9-processed-gene-2.0-mRNA-1:three_prime_utr;Parent=genemark-Pcoprophilum_scaf_9-processed-gene-2.0-mRNA-1
Pcoprophilum_scaf_9 maker gene  197438  198714  . - . ID=genemark-Pcoprophilum_scaf_9-processed-gene-2.3;Name=genemark-Pcoprophilum_scaf_9-processed-gene-2.3
Pcoprophilum_scaf_9 maker mRNA  197438  198714  . - . ID=genemark-Pcoprophilum_scaf_9-processed-gene-2.3-mRNA-1;Parent=genemark-Pcoprophilum_scaf_9-processed-gene-2.3;Name=genemark-Pcoprophilum_scaf_9-processed-gene-2.3-mRNA-1;_AED=0.00;_eAED=0.00;_QI=208|1|1|1|1|1|2|259|211
Pcoprophilum_scaf_9 maker exon  197438  198116  . - . ID=genemark-Pcoprophilum_scaf_9-processed-gene-2.3-mRNA-1:exon:141;Parent=genemark-Pcoprophilum_scaf_9-processed-gene-2.3-mRNA-1
Pcoprophilum_scaf_9 maker exon  198291  198714  . - . ID=genemark-Pcoprophilum_scaf_9-processed-gene-2.3-mRNA-1:exon:140;Parent=genemark-Pcoprophilum_scaf_9-processed-gene-2.3-mRNA-1
Pcoprophilum_scaf_9 maker five_prime_UTR  198507  198714  . - . ID=genemark-Pcoprophilum_scaf_9-processed-gene-2.3-mRNA-1:five_prime_utr;Parent=genemark-Pcoprophilum_scaf_9-processed-gene-2.3-mRNA-1
Pcoprophilum_scaf_9 maker CDS 198291  198506  . - 0 ID=genemark-Pcoprophilum_scaf_9-processed-gene-2.3-mRNA-1:cds;Parent=genemark-Pcoprophilum_scaf_9-processed-gene-2.3-mRNA-1
Pcoprophilum_scaf_9 maker CDS 197697  198116  . - 0 ID=genemark-Pcoprophilum_scaf_9-processed-gene-2.3-mRNA-1:cds;Parent=genemark-Pcoprophilum_scaf_9-processed-gene-2.3-mRNA-1
Pcoprophilum_scaf_9 maker three_prime_UTR 197438  197696  . - . ID=genemark-Pcoprophilum_scaf_9-processed-gene-2.3-mRNA-1:three_prime_utr;Parent=genemark-Pcoprophilum_scaf_9-processed-gene-2.3-mRNA-1
EXAMPLE WITH SPREADED START AND STOP
##gff-version 3
Pcoprophilum_scaf_9 maker gene  189352  192747  . + . ID=genemark-Pcoprophilum_scaf_9-processed-gene-2.0;Name=genemark-Pcoprophilum_scaf_9-processed-gene-2.0
Pcoprophilum_scaf_9 maker mRNA  189352  192747  . + . ID=genemark-Pcoprophilum_scaf_9-processed-gene-2.0-mRNA-1;Parent=genemark-Pcoprophilum_scaf_9-processed-gene-2.0;Name=genemark-Pcoprophilum_scaf_9-processed-gene-2.0-mRNA-1;_AED=0.00;_eAED=0.00;_QI=398|1|1|1|0.5|0.33|3|343|825
Pcoprophilum_scaf_9 maker exon  189352  189520  . + . ID=genemark-Pcoprophilum_scaf_9-processed-gene-2.0-mRNA-1:exon:96;Parent=genemark-Pcoprophilum_scaf_9-processed-gene-2.0-mRNA-1
Pcoprophilum_scaf_9 maker exon  189643  189922  . + . ID=genemark-Pcoprophilum_scaf_9-processed-gene-2.0-mRNA-1:exon:97;Parent=genemark-Pcoprophilum_scaf_9-processed-gene-2.0-mRNA-1
Pcoprophilum_scaf_9 maker exon  189978  192747  . + . ID=genemark-Pcoprophilum_scaf_9-processed-gene-2.0-mRNA-1:exon:98;Parent=genemark-Pcoprophilum_scaf_9-processed-gene-2.0-mRNA-1
Pcoprophilum_scaf_9 maker five_prime_UTR  189352  189520  . + . ID=genemark-Pcoprophilum_scaf_9-processed-gene-2.0-mRNA-1:five_prime_utr;Parent=genemark-Pcoprophilum_scaf_9-processed-gene-2.0-mRNA-1
Pcoprophilum_scaf_9 maker five_prime_UTR  189643  189871  . + . ID=genemark-Pcoprophilum_scaf_9-processed-gene-2.0-mRNA-1:five_prime_utr;Parent=genemark-Pcoprophilum_scaf_9-processed-gene-2.0-mRNA-1
Pcoprophilum_scaf_9 maker CDS 189872  189873  . + 0 ID=genemark-Pcoprophilum_scaf_9-processed-gene-2.0-mRNA-1:cds;Parent=genemark-Pcoprophilum_scaf_9-processed-gene-2.0-mRNA-1
Pcoprophilum_scaf_9 maker CDS 189874  189922  . + 0 ID=genemark-Pcoprophilum_scaf_9-processed-gene-2.0-mRNA-1:cds;Parent=genemark-Pcoprophilum_scaf_9-processed-gene-2.0-mRNA-1
Pcoprophilum_scaf_9 maker CDS 189978  192402  . + 0 ID=genemark-Pcoprophilum_scaf_9-processed-gene-2.0-mRNA-1:cds;Parent=genemark-Pcoprophilum_scaf_9-processed-gene-2.0-mRNA-1
Pcoprophilum_scaf_9 maker CDS 192403  192404  . + 0 ID=genemark-Pcoprophilum_scaf_9-processed-gene-2.0-mRNA-1:cds;Parent=genemark-Pcoprophilum_scaf_9-processed-gene-2.0-mRNA-1
Pcoprophilum_scaf_9 maker three_prime_UTR 192405  192747  . + . ID=genemark-Pcoprophilum_scaf_9-processed-gene-2.0-mRNA-1:three_prime_utr;Parent=genemark-Pcoprophilum_scaf_9-processed-gene-2.0-mRNA-1
Pcoprophilum_scaf_9 maker gene  197438  198714  . - . ID=genemark-Pcoprophilum_scaf_9-processed-gene-2.3;Name=genemark-Pcoprophilum_scaf_9-processed-gene-2.3
Pcoprophilum_scaf_9 maker mRNA  197438  198714  . - . ID=genemark-Pcoprophilum_scaf_9-processed-gene-2.3-mRNA-1;Parent=genemark-Pcoprophilum_scaf_9-processed-gene-2.3;Name=genemark-Pcoprophilum_scaf_9-processed-gene-2.3-mRNA-1;_AED=0.00;_eAED=0.00;_QI=208|1|1|1|1|1|2|259|211
Pcoprophilum_scaf_9 maker exon  197438  198116  . - . ID=genemark-Pcoprophilum_scaf_9-processed-gene-2.3-mRNA-1:exon:141;Parent=genemark-Pcoprophilum_scaf_9-processed-gene-2.3-mRNA-1
Pcoprophilum_scaf_9 maker exon  198291  198714  . - . ID=genemark-Pcoprophilum_scaf_9-processed-gene-2.3-mRNA-1:exon:140;Parent=genemark-Pcoprophilum_scaf_9-processed-gene-2.3-mRNA-1
Pcoprophilum_scaf_9 maker five_prime_UTR  198507  198714  . - . ID=genemark-Pcoprophilum_scaf_9-processed-gene-2.3-mRNA-1:five_prime_utr;Parent=genemark-Pcoprophilum_scaf_9-processed-gene-2.3-mRNA-1
Pcoprophilum_scaf_9 maker CDS 198505  198506  . - 0 ID=genemark-Pcoprophilum_scaf_9-processed-gene-2.3-mRNA-1:cds;Parent=genemark-Pcoprophilum_scaf_9-processed-gene-2.3-mRNA-1
Pcoprophilum_scaf_9 maker CDS 198291  198504  . - 0 ID=genemark-Pcoprophilum_scaf_9-processed-gene-2.3-mRNA-1:cds;Parent=genemark-Pcoprophilum_scaf_9-processed-gene-2.3-mRNA-1
Pcoprophilum_scaf_9 maker CDS 197699  198116  . - 0 ID=genemark-Pcoprophilum_scaf_9-processed-gene-2.3-mRNA-1:cds;Parent=genemark-Pcoprophilum_scaf_9-processed-gene-2.3-mRNA-1
Pcoprophilum_scaf_9 maker CDS 197697  197698  . - 0 ID=genemark-Pcoprophilum_scaf_9-processed-gene-2.3-mRNA-1:cds;Parent=genemark-Pcoprophilum_scaf_9-processed-gene-2.3-mRNA-1
Pcoprophilum_scaf_9 maker three_prime_UTR 197438  197696  . - . ID=genemark-Pcoprophilum_scaf_9-processed-gene-2.3-mRNA-1:three_prime_utr;Parent=genemark-Pcoprophilum_scaf_9-processed-gene-2.3-mRNA-1


=head1 NAME

agat_sp_add_start_and_stop.pl.pl

=head1 DESCRIPTION

The script adds start and stop codons when a CDS feature exists.
The script looks at the nucleotide sequence and checks the presence of start and stop codons.
The script works even if the start or stop codon are split over several CDS features.

=head1 SYNOPSIS

    agat_sp_add_start_and_stop.pl.pl --gff infile.gff --fasta genome.fa --out outfile.gff
    agat_sp_add_start_and_stop.pl.pl --help

=head1 OPTIONS

=over 8

=item B<--gff>, B<-i> or B<-g>

Input GTF/GFF file.

=item B<--fasta>, B<--fa> or B<-f>

Input fasta file. Needed to check that CDS sequences start by start codon and stop by stop codon.

=item B<--ct>, B<--codon> or B<--table>

Codon table to use. [default 1]

=item B<--out>, B<--output> or B<-o>

Output gff file updated

=item B<-e> or B<--extend>

Boolean - When no start/stop codon found, try to extend the CDS to meet the next start/stop codon in the sequence. 

=item B<--ni> or B<--na>

Boolean - no iupac / no ambiguous, avoid usage of IUPAC. By default IUPAC is used that means, NNN is seen as start and/or stop codon.

=item B<-v> or B<--verbose>

Verbose for debugging purpose.

=item B<-thread>, B<threads>, B<cpu>, B<cpus>, B<core>, B<cores>, B<job> or B<jobs>

Integer - Number of parallel processes to use for file input parsing (via forking).

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
