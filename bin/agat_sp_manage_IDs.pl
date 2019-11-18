#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Bio::Tools::GFF;
use Agat::Omniscient;

my $header = get_agat_header();
my $gff = undef;
my $opt_help= 0;
my @opt_tag=();
my $outfile=undef;
my $prefix=undef;
my $nbIDstart=1;

if ( !GetOptions(
    "help|h" => \$opt_help,
    "gff|f=s" => \$gff,
    "nb=i" => \$nbIDstart,
    "prefix=s" => \$prefix,
    "p|t|l=s" => \@opt_tag,
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

if ( ! (defined($gff)) ){
    pod2usage( {
           -message => "$header\nAt least 1 parameter is mandatory:\nInput reference gff file (--gff) \n\n",
           -verbose => 0,
           -exitval => 2 } );
}

my $gffout;
if ($outfile) {
  $outfile=~ s/.gff//g;
  open(my $fh, '>', $outfile.".gff") or die "Could not open file '$outfile' $!";
  $gffout= Bio::Tools::GFF->new(-fh => $fh, -gff_version => 3 );
}
else{
  $gffout = Bio::Tools::GFF->new(-fh => \*STDOUT, -gff_version => 3);
}

# Manage $primaryTag
my %ptagList;
if(! @opt_tag){
  print "We will work on attributes from all features\n";
  $ptagList{'level1'}++;
  $ptagList{'level2'}++;
  $ptagList{'level3'}++;
}
else{
  foreach my $tag (@opt_tag){
    if($tag eq ""){next;}
    if($tag eq "all"){
      print "We will work on attributes from all features\n";
      $ptagList{'level1'}++;
      $ptagList{'level2'}++;
      $ptagList{'level3'}++;
    }
    else{
      print "We will work on attributes from all the $tag features\n";
      $ptagList{lc($tag)}++;
    }
  }
}




                #####################
                #     MAIN          #
                #####################

my %keepTrack;
my %tag_hash;
my @tagLetter_list;

######################
### Parse GFF input #
my ($hash_omniscient, $hash_mRNAGeneLink) = slurp_gff3_file_JD({ input => $gff
                                                              });
print ("GFF3 file parsed\n");

# sort by seq id
my %hash_sortBySeq;
foreach my $tag_level1 ( keys %{$hash_omniscient->{'level1'}}){
  foreach my $level1_id ( keys %{$hash_omniscient->{'level1'}{$tag_level1}}){
    my $position=$hash_omniscient->{'level1'}{$tag_level1}{$level1_id}->seq_id;
    push (@{$hash_sortBySeq{$position}{$tag_level1}}, $hash_omniscient->{'level1'}{$tag_level1}{$level1_id});
  }
}

#Read by seqId to sort properly for ID naming
foreach my $seqid (sort alphaNum keys %hash_sortBySeq){ # loop over all the feature level1

  foreach my $tag_l1 (sort {$a cmp $b} keys %{$hash_sortBySeq{$seqid}}){

    foreach my $feature_l1 ( sort {$a->start <=> $b->start} @{$hash_sortBySeq{$seqid}{$tag_l1}}){
      my $id_l1 = lc($feature_l1->_tag_value('ID'));
      my $l1_ID_modified=undef;

      if(exists ($ptagList{$tag_l1}) or  exists ($ptagList{'level1'}) ){
        if(! exists_keys(\%keepTrack,($tag_l1))){$keepTrack{$tag_l1}=$nbIDstart;}
        manage_attributes($feature_l1,\%keepTrack, $prefix);
        $keepTrack{$tag_l1}++;
        $l1_ID_modified=$feature_l1->_tag_value('ID');
        $hash_omniscient->{'level1'}{$tag_l1}{lc($l1_ID_modified)} = delete $hash_omniscient->{'level1'}{$tag_l1}{$id_l1};
      }

      #################
      # == LEVEL 2 == #
      #################
      foreach my $tag_l2 (sort {$a cmp $b}  keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...

        if ( exists ($hash_omniscient->{'level2'}{$tag_l2}{$id_l1} ) ){
          foreach my $feature_l2 ( sort {$a->start <=> $b->start} @{$hash_omniscient->{'level2'}{$tag_l2}{$id_l1}}) {

            my $l2_ID_modified=undef;
            my $level2_ID = lc($feature_l2->_tag_value('ID'));

            if(exists ($ptagList{$tag_l2}) or  exists ($ptagList{'level2'}) ){
              if(! exists_keys(\%keepTrack,($tag_l2))){$keepTrack{$tag_l2}=$nbIDstart;}
              manage_attributes($feature_l2,\%keepTrack, $prefix);
              $keepTrack{$tag_l2}++;
              $l2_ID_modified=$feature_l2->_tag_value('ID');
            }

            #Modify parent if necessary
            if($l1_ID_modified){
               create_or_replace_tag($feature_l2,'Parent', $l1_ID_modified);
            }

            #################
            # == LEVEL 3 == #
            #################
            foreach my $tag_l3 (sort {$a cmp $b} keys %{$hash_omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...

              if ( exists_keys($hash_omniscient, ('level3', $tag_l3 , $level2_ID) ) ){

                foreach my $feature_l3 ( sort {$a->start <=> $b->start} @{$hash_omniscient->{'level3'}{$tag_l3}{$level2_ID}}) {

                  if(exists ($ptagList{$tag_l3}) or  exists ($ptagList{'level3'}) ){
                    if(! exists_keys(\%keepTrack,($tag_l3))){$keepTrack{$tag_l3}=$nbIDstart;}
                    manage_attributes($feature_l3,\%keepTrack, $prefix);
                    $keepTrack{$tag_l3}++;
                  }

                  #Modify parent if necessary
                  if($l2_ID_modified){
                     create_or_replace_tag($feature_l3,'Parent', $l2_ID_modified);
                  }

                }

                if($l2_ID_modified){
                  $hash_omniscient->{'level3'}{$tag_l3}{lc($l2_ID_modified)} = delete $hash_omniscient->{'level3'}{$tag_l3}{$level2_ID};
                }
              }
            }
          }
          if($l1_ID_modified){
            $hash_omniscient->{'level2'}{$tag_l2}{lc($l1_ID_modified)} = delete $hash_omniscient->{'level2'}{$tag_l2}{$id_l1};
          }
        }
      }
    }
  }
}

# Print results
print_omniscient($hash_omniscient, $gffout); #print gene modified


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
  my  ($feature, $keepTrack, $prefix)=@_;

  my $primary_tag = lc($feature->primary_tag);

  if ($prefix){

    my $nbName = $keepTrack->{$primary_tag};

    my $numberNum=11;
    my $GoodNum="";
    for (my $i=0; $i<$numberNum-length($nbName); $i++){
      $GoodNum.="0";
    }
    $GoodNum.=$nbName;

    my $abb = uc(select_abb($feature));

    my $result="$prefix$abb$GoodNum";
    create_or_replace_tag($feature,'ID', $result);
  }
  else{
    create_or_replace_tag($feature,'ID', $primary_tag."-".$keepTrack->{$primary_tag});
  }
}

#Select the proper abbreviation for the tag
sub  select_abb{
  my  ($feature)=@_;

  # get the tag
  my $primary_tag = lc($feature->primary_tag);

  if(! exists_keys (\%tag_hash,( $primary_tag ))) {

    my $cpt=1;
    my $letter = uc(substr($primary_tag, 0, $cpt));

    while(  grep( /^\Q$letter\E$/, @tagLetter_list) ) { # to avoid duplicate
      $cpt++;
      $letter = uc(substr($primary_tag, 0, $cpt));
    }
    $tag_hash{$primary_tag}=$letter;
    push(@tagLetter_list, $letter)
  }
  return $tag_hash{$primary_tag}
}

#Sorting mixed strings => Sorting alphabetically first, then numerically
# how to use: my @y = sort by_number @x;
sub alphaNum {
    my ( $alet , $anum ) = $a =~ /([^\d]+)(\d+)/;
    my ( $blet , $bnum ) = $b =~ /([^\d]+)(\d+)/;
    ( $alet || "a" ) cmp ( $blet || "a" ) or ( $anum || 0 ) <=> ( $bnum || 0 )
}

__END__

=head1 NAME

agat_sp_manage_IDs.pl

=head1 DESCRIPTION

The script take a gff3 file as input and will go through all feature to overwrite the uniq ID.
By default the ID is build as follow:
  primary_tag(i.e. 3rd column)-Number.
If you provide a specific prefix the ID is build as follow (Ensembl like format ENSG00000000022):
 $prefix.$letterCode.0*.Number where the number of 0 i adapted in order to have 11 digits

By default the numbering start to 1, but you can decide to change this value using the --nb option.
The $letterCode is generated on the fly to be uniq. By defaut it used the first letter of the feature type (3rd colum). If two feature types
start with the same letter, the second one meet will have the two first letter as $letterCode (and so one).

=head1 SYNOPSIS

    agat_sp_manage_IDs.pl -gff file.gff -p level2 -p cds -p exon [ -o outfile ]
    agat_sp_manage_IDs.pl --help

=head1 OPTIONS

=over 8

=item B<--gff> or B<-f>

Input GTF/GFF file.

=item B<--prefix>

String. Add a specific prefix to the ID.

=item B<-p>,  B<-t> or  B<-l>

primary tag option, case insensitive, list. Allow to specied the feature types that will be handled.
You can specified a specific feature by given its primary tag name (column 3) as: cds, Gene, MrNa
You can specify directly all the feature of a particular level:
      level2=mRNA,ncRNA,tRNA,etc
      level3=CDS,exon,UTR,etc
By default all feature are taken into account. fill the option by the value "all" will have the same behaviour.

=item B<--nb>

Integer. Start numbering to this value.

=item B<-o> , B<--output> , B<--out> or B<--outfile>

Output GFF file.  If no output file is specified, the output will be
written to STDOUT.

=item B<-h> or B<--help>

Display this helpful text.

=back

=cut
