#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Sort::Naturally;
use Pod::Usage;
use AGAT::AGAT;

my $header = get_agat_header();
my $config;
my $opt_gff = undef;
my $opt_help= 0;
my $opt_gap=0;
my $opt_tair=undef;
my @opt_tag=();
my $outfile=undef;
my $opt_ensembl=undef;
my $opt_prefix=undef;
my $opt_collective=undef;
my $opt_nbIDstart=1;
my $opt_type_dependent = undef;
my $verbose;

if ( !GetOptions(
    'c|config=s'     => \$config,
    "h|help!"        => \$opt_help,
    "gff|f=s"        => \$opt_gff,
    "nb=i"           => \$opt_nbIDstart,
    "gap=i"          => \$opt_gap,
    "tair!"          => \$opt_tair,
    "ensembl!"       => \$opt_ensembl,
    "prefix=s"       => \$opt_prefix,
    "p|t|l=s"        => \@opt_tag,
    "type_dependent!" => \$opt_type_dependent,
		"collective!"    => \$opt_collective,
    "verbose|v!"     => \$verbose,
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

if ( ! (defined($opt_gff)) ){
    pod2usage( {
           -message => "$header\nAt least 1 parameter is mandatory:\nInput reference gff file (--gff) \n\n",
           -verbose => 0,
           -exitval => 1 } );
}

# --- Manage config ---
$config = get_agat_config({config_file_in => $config});

my $gffout = prepare_gffout($config, $outfile);

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
my @l3_out_priority = ("tss", "exon", "cds", "tts");
######################
### Parse GFF input #
my ($hash_omniscient, $hash_mRNAGeneLink) = slurp_gff3_file_JD({ input => $opt_gff,
                                                                 config => $config
                                                            });
print ("GFF3 file parsed\n");

# get spreadfeatire in case of collective option set
my $spreadfeatures = $hash_omniscient->{'other'}{'level'}{'spread'};

# sort by seq id
my $hash_sortBySeq = gather_and_sort_l1_by_seq_id($hash_omniscient);

my $opt_tair_suffix=0;
#Read by seqId to sort properly for ID naming

foreach my $seqid ( sort { ncmp ($a, $b) } keys %{$hash_sortBySeq}){ # loop over all the feature level1

  foreach my $tag_l1 (sort {$a cmp $b} keys %{$hash_sortBySeq{$seqid}}){

    foreach my $feature_l1 ( @{$hash_sortBySeq{$seqid}{$tag_l1}}){ # feature are alredy sorted by function that made that hash
      my $id_l1 = lc($feature_l1->_tag_value('ID'));
      my $l1_ID_modified=undef;

      if(exists ($ptagList{$tag_l1}) or  exists ($ptagList{'level1'}) ){
        manage_attributes('level1', $feature_l1);
        $l1_ID_modified=$feature_l1->_tag_value('ID');
        $hash_omniscient->{'level1'}{$tag_l1}{lc($l1_ID_modified)} = delete $hash_omniscient->{'level1'}{$tag_l1}{$id_l1};
      }

      #################
      # == LEVEL 2 == #
      #################
      $opt_tair_suffix=0;
      foreach my $tag_l2 (sort {$a cmp $b}  keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...

        if ( exists ($hash_omniscient->{'level2'}{$tag_l2}{$id_l1} ) ){
          foreach my $feature_l2 ( sort { ncmp ($a->start."|".$a->end.$a->_tag_value('ID'), $b->start."|".$b->end.$b->_tag_value('ID') ) } @{$hash_omniscient->{'level2'}{$tag_l2}{$id_l1}}) {
            $opt_tair_suffix++;
            my $l2_ID_modified=undef;
            my $level2_ID = lc($feature_l2->_tag_value('ID'));
            #Modify parent if necessary
            if($l1_ID_modified){
               create_or_replace_tag($feature_l2,'Parent', $l1_ID_modified);
            }

            if(exists ($ptagList{$tag_l2}) or  exists ($ptagList{'level2'}) ){
              manage_attributes('level2', $feature_l2, );
              $l2_ID_modified=$feature_l2->_tag_value('ID');
            }

            #################
            # == LEVEL 3 == #
            #################
            ##########
            # Same order as in OmniscientO
            if ( exists_keys($hash_omniscient,('level3','tss',$level2_ID)) ){
              deal_with_level3(\%ptagList, $level2_ID, 'tss', $l2_ID_modified );
            }
            if ( exists_keys( $hash_omniscient, ('level3', 'exon', $level2_ID) ) ){
              deal_with_level3(\%ptagList, $level2_ID, 'exon', $l2_ID_modified );
            }
            if ( exists_keys( $hash_omniscient, ('level3', 'cds', $level2_ID) ) ){
              deal_with_level3(\%ptagList, $level2_ID, 'cds', $l2_ID_modified );
            }
            if ( exists_keys($hash_omniscient,('level3','tts',$level2_ID)) ){
              deal_with_level3(\%ptagList, $level2_ID, 'tts', $l2_ID_modified );
            }
            foreach my $tag_l3 (sort {$a cmp $b} keys %{$hash_omniscient->{'level3'}}){
              if (! grep { $_ eq $tag_l3 } @l3_out_priority){
                if ( exists_keys($hash_omniscient, ('level3', $tag_l3 , $level2_ID) ) ){
                  deal_with_level3(\%ptagList, $level2_ID, $tag_l3, $l2_ID_modified );
                }
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

# Print results
print_omniscient( {omniscient => $hash_omniscient, output => $gffout} );

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
sub deal_with_level3{
  my  ($ptagList, $level2_ID, $tag_l3, $l2_ID_modified)=@_;

	my $first_round = 1;
	my $previous_id = undef;

  foreach my $feature_l3 ( sort {$a->start <=> $b->start} @{$hash_omniscient->{'level3'}{$tag_l3}{$level2_ID}}) {

    #Modify parent if necessary
    if($l2_ID_modified){
       create_or_replace_tag($feature_l3,'Parent', $l2_ID_modified);
    }
    if(exists ($ptagList->{$tag_l3}) or  exists ($ptagList->{'level3'}) ){

			# Case withcollective ID (need option activated + the feature type is a spreadfeature)
			if($opt_collective and exists_keys($spreadfeatures, ($tag_l3) ) ){
				if(! $previous_id){
					manage_attributes('level3', $feature_l3);
					$previous_id = $feature_l3->_tag_value('ID');
				}
				else{
					create_or_replace_tag($feature_l3,'ID', $previous_id);
				}
			}
			# Case without collective ID
			else{
      	manage_attributes('level3', $feature_l3);
			}
    }
  }

  if($l2_ID_modified){
    $hash_omniscient->{'level3'}{$tag_l3}{lc($l2_ID_modified)} = delete $hash_omniscient->{'level3'}{$tag_l3}{$level2_ID};
  }
}

# $opt_prefix.$letterCode.EnsemblSup.$Number
# $primary_tag
# $opt_tair_suffix $opt_nbIDstart $keepTrack, $opt_type_dependent, $opt_prefix, $opt_ensembl, $opt_tair, $opt_tair_suffix variable are available from everywhere in the script
sub  manage_attributes{
  my  ($level, $feature)=@_;

  my $result;
  my $primary_tag = lc($feature->primary_tag);
  my $prefix = undef;
  my $parent_id = undef;
  if ($level ne 'level1'){ $parent_id = $feature->_tag_value('Parent'); }


  # ---- deal with prefix ----
  if (! $opt_prefix){ # Either the one given or the primary_tag by default
    $prefix = $primary_tag."-";
  }
  else{
    $prefix = $opt_prefix;
  }
  print "prefix $prefix \n" if $verbose;

  #  ----- deal with value independent or not of the feature type--------
  my $tag = $opt_type_dependent ? $primary_tag : 'all';
  print "tag: $tag\n" if $verbose;
  print "primary_tag: $primary_tag\n" if $verbose;

  if ($opt_tair){
    if ($level eq 'level1') {
      my $ID_number = get_id_number($tag, $level);
      $ID_number = add_ensembl_id_number_prefix($feature, $ID_number) if ($opt_ensembl);
      $result="$prefix$ID_number";
      print "tair l1: $result\n" if $verbose;
    }
    if($level eq 'level2'){
      $result = $parent_id.".".$opt_tair_suffix; # add .1, .2,etc to level features
        print "tair l2: $result\n" if $verbose;
    }
    elsif ($level eq 'level3'){
      my $ID_number = get_id_number("$parent_id$tag", $level);
      $result = $parent_id."-"."$primary_tag$ID_number";
      print "tair l3: $result\n" if $verbose;
    }

  }
  else{
    # ---- get id number -----
    my $ID_number = get_id_number($tag, $level);
    # ---- deal with Ensembl - depend of value -----
    $ID_number = add_ensembl_id_number_prefix($feature, $ID_number) if ($opt_ensembl);
    # Finalize
    $result="$prefix$ID_number";
  }

  create_or_replace_tag($feature,'ID', $result);
}

# get ensembl id number prefix according to ID_number
sub add_ensembl_id_number_prefix{
  my  ($feature, $ID_number)=@_;

  my $abb = uc(select_abb($feature));
  my $GoodNum="";
  for (my $i=0; $i<11-length($ID_number); $i++){
    $GoodNum.="0";
  }

  return $abb.$GoodNum.$ID_number; #E00000000001
}

# Get id number according to tag
sub get_id_number{
  my  ($tag, $level)=@_;

  # First time for the tag
  if(! exists_keys(\%keepTrack,($tag))){ $keepTrack{$tag}=$opt_nbIDstart;}
  # Not first time
  else{# if Gap asked we add this value between to level1 feature
    if ($opt_gap and $level eq 'level1'){
      $keepTrack{$tag} += $opt_gap+1;
    } # normal incrementation
    else{
      $keepTrack{$tag}++;
    }
  }

  return $keepTrack{$tag};
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

The script takes a gff3 file as input and will go through all feature to overwrite
the value of the ID attribute.
By default the ID is built as follow: primary_tag(i.e. 3rd column)-Number.
If you provide a specific prefix the ID is built as follow: $opt_prefix.$letterCode.Number.
By default the numbering start at 1, but you can decide to change this value using the --nb option.
The $letterCode is the first letter of the feature type (3rd colum). It is uniq for each feature type,
i.e. when two feature types start with the same letter, the second one met will have the two first letter as $letterCode (and so one).

=head1 SYNOPSIS

    agat_sp_manage_IDs.pl --gff file.gff -p level2 -p cds -p exon [ -o outfile ]
    agat_sp_manage_IDs.pl --help

=head1 OPTIONS

=over 8

=item B<--gff> or B<-f>

Input GTF/GFF file.

=item B<--gap>

Integer - Increment the next gene (level1 feature) suffix with this value. Defauft 0.

=item B<--ensembl>

Boolean - For an ID Ensembl like (e.g PREFIXG00000000022). The ID is built as follow:
$opt_prefix.$letterCode.0*.Number where the number of 0 is adapted in order to have 11 digits.

=item B<--prefix>

String - Add a specific prefix to the ID. By defaut if will be the feature type (3rd column).

=item B<--type_dependent>

Boolean - Activate type_dependent numbering. The number is depedendent of the feature type.
i.e instead of:
NbV1Ch01        AUGUSTUS        gene    97932   99714   0.06    -       .       ID=gene1
NbV1Ch01        AUGUSTUS        mRNA    97932   99714   0.06    -       .       ID=mRNA2
NbV1Ch01        AUGUSTUS        exon    97932   98571   .       -       .       ID=exon3
NbV1Ch01        AUGUSTUS        exon    98679   98844   .       -       .       ID=exon4
You will get:
NbV1Ch01        AUGUSTUS        gene    97932   99714   0.06    -       .       ID=gene1
NbV1Ch01        AUGUSTUS        mRNA    97932   99714   0.06    -       .       ID=mRNA1
NbV1Ch01        AUGUSTUS        exon    97932   98571   .       -       .       ID=exon1
NbV1Ch01        AUGUSTUS        exon    98679   98844   .       -       .       ID=exon2

=item B<--collective>

Boolean - In the case of discontinuous features (i.e. a single feature that
exists over multiple genomic locations like CDS, UTR) we set a uniq ID by default.
If you wish to set the a collective ID for those feature, please activate this option.

=item B<--tair>

Boolean - Tair like Output:

NbV1Ch01    TAIR10  gene    5928    8737    .       -       .       ID=AT1G01020
NbV1Ch01    TAIR10  mRNA    5928    8737    .       -       .       ID=AT1G01020.1
NbV1Ch01    TAIR10  exon    5928    8737   .       -       .        ID=AT1G01020.1-exon1

=item B<--nb>

Integer - Start numbering to this value. Default 1.

=item B<-p>,  B<-t> or  B<-l>

primary tag option, case insensitive, list. Allow to specied the feature types that will be handled.
You can specified a specific feature by given its primary tag name (column 3) as: cds, Gene, MrNa
You can specify directly all the feature of a particular level:
      level2=mRNA,ncRNA,tRNA,etc
      level3=CDS,exon,UTR,etc
By default all feature are taken into account. fill the option by the value "all" will have the same behaviour.

=item B<-o> , B<--output> , B<--out> or B<--outfile>

String - Output GFF file. If no output file is specified, the output will be
written to STDOUT.

=item B<-c> or B<--config>

String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any, 
otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

=item B<-h> or B<--help>

Boolean - Display this helpful text.

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
