#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use Carp;
use POSIX qw(strftime);
use Getopt::Long;
use IO::File;
use Pod::Usage;
use List::MoreUtils qw(uniq);
use Bio::DB::Fasta;
use AGAT::AGAT;

#use Data::Dumper; # JN: for dedug printing
my $DEBUG = 0;    # JN: for dedug printing

my $header = get_agat_header();
my $config;
my $cpu;

# PARAMETERS - OPTION
my $opt_reffile;
my $opt_output;
my $opt_BlastFile;
my $opt_CleanNameAttribute; # Should we remove the Name attribute value if already exists - bolean
my $opt_CleanProductAttribute; # Should we remove the product attribute value if already exists - bolean
my $opt_CleanOntology_termAttribute; # Should we remove the Ontology_term attribute value if already exists - bolean
my $opt_CleanDbxrefAttribute; # Should we remove the Dbxref attribute value if already exists - bolean
my $opt_InterproFile;
my $opt_name = undef;
my $opt_nameU;
my $opt_populate_cds = undef;
my $opt_verbose = undef;
my $opt_help = 0;
my $opt_blastEvalue = 1e-6;
my $opt_dataBase = undef;
my $opt_pe = 5;
my $opt_addGnPresentTag = 0; # JN: Optionally add the 'gn_present=yes|no|NA' tag in gff
my %numbering;
my $nbIDstart = 1;
my $prefixName = undef;
my %tag_hash;
my @tag_list;
# END PARAMETERS - OPTION

# FOR FUNCTIONS BLAST#
my %nameBlast;
my %geneNameBlast;
my %mRNANameBlast;
my %mRNAproduct;
my %mRNAUniprotIDFromBlast; # Uniprot ID of the best candidate (from balst file)
my %blast_evalue;           # evalue of the best candidate (from balst file)
my %blast_organism;         # organism (OS) of the best candidate (from balst file)
my %blast_protEvidence;     # Protein Evidence (PR) of the best candidate (from balst file)
my %blast_seqVersion;       # Sequence Version (SV) of the best candidate (from balst file)
my %geneNameGiven;          # Gene Name (GN) of the best candidate (from balst file)
my %duplicateNameGiven;
my %fasta_id_gn_hash = ();     # JN: key: 'sp|a6w1c3|hem1_marms' , value: 'hema' etc or undef
my %l2_gn_present_hash = ();   # JN: Key: 'maker-bi03_p1mp_001088f-est_gff_stringtie-gene-0.2-mrna-1', value: 'yes' or 'no'
my $nbGnNotPresentInDb = 0;    # JN: Count entries without GN in db
my $nbGnNotPresentForMrna = 0; # JN: Count mRNAs without GN in db
my $nbDuplicateNameGiven = 0;
my $nbNamedGene = 0;
# END FOR FUNCTION BLAST#

# FOR FUNCTIONS INTERPRO#
my %TotalTerm;
my %finalID;
my %GeneAssociatedToTerm;
my %mRNAAssociatedToTerm;
my %functionData;
my %functionDataAdded;
my %functionOutput;
my %functionStreamOutput;
my %geneWithoutFunction;
my %geneWithFunction;
my $nbmRNAwithoutFunction = 0;
my $nbmRNAwithFunction = 0;
my $nbGeneWithGOterm = 0;
my $nbTotalGOterm = 0;
# END FOR FUNCTION INTERPRO#

# OPTION MANAGMENT
my @copyARGV = @ARGV;
GetOptions(
 'f|ref|reffile|gff|gff3=s' => \$opt_reffile,
 'b|blast=s'                => \$opt_BlastFile,
 'clean_name!'              => \$opt_CleanNameAttribute,
 'clean_product!'           => \$opt_CleanProductAttribute,
 'clean_dbxref!'            => \$opt_CleanDbxrefAttribute,
 'clean_ontology!'          => \$opt_CleanOntology_termAttribute,
 'd|db=s'                   => \$opt_dataBase,
 'be|blast_evalue=f'        => \$opt_blastEvalue,
 'pe=i'                     => \$opt_pe,
 'pcds!'                    => \$opt_populate_cds,
 'i|interpro=s'             => \$opt_InterproFile,
 'id=s'                     => \$opt_name,
 'idau=s'                   => \$opt_nameU,
 'nb=i'                     => \$nbIDstart,
 'o|output=s'               => \$opt_output,
 'a|addgntag'               => \$opt_addGnPresentTag,
 'v'                        => \$opt_verbose,
 'c|config=s'               => \$config,
 'thread|threads|cpu|cpus|core|cores|job|jobs=i' => \$cpu,
 'h|help!'                  => \$opt_help
)
or pod2usage( {
  -message => 'Failed to parse command line',
  -verbose => 1,
  -exitval => 1
});

# Print Help and exit
if ($opt_help) {
  pod2usage(
    {
      -verbose => 99,
      -exitval => 0,
      -message => "$header\n"
    }
  );
}

if ( !( defined($opt_reffile) ) ) {
  pod2usage(
    {
      -message => "$header\nAt least 1 parameter is mandatory:\nInput reference gff file (--f)\n\n"
        . "Many optional parameters are available. Look at the help documentation to know more.\n",
      -verbose => 0,
      -exitval => 1
    }
  );
}

# --- Manage config ---
initialize_agat({ config_file_in => $config, input => $opt_reffile });
$CONFIG->{cpu} = $cpu if defined($cpu);

#################################################
####### START Manage files (input output) #######
#################################################

if ( ($opt_pe > 5) or ($opt_pe < 1) ) {
  print "Error the Protein Existence (PE) value must be between 1 and 5\n";
  exit;
}

my $streamBlast = IO::File->new();
my $streamInter = IO::File->new();

# Manage Blast File
if (defined $opt_BlastFile) {
  if (! $opt_dataBase) {
    print "To use the blast output we also need the fasta of the database used for the blast (--db)\n";
    exit;
  }
  $streamBlast->open( $opt_BlastFile, 'r' ) or
    croak( sprintf( "Can not open '%s' for reading: %s", $opt_BlastFile, $! ) );
}

# Manage Interpro file
if (defined $opt_InterproFile) {
  $streamInter->open( $opt_InterproFile, 'r' ) or
    croak( sprintf( "Can not open '%s' for reading: %s", $opt_InterproFile, $! ) );
}

##########################
##### Manage Output ######

my $ostreamGFF_file;
my $ostreamLog_file;
my $ostreamReport_file;
if (defined($opt_output)) {
  if (-f $opt_output) {
    print "Cannot create a directory with the name $opt_output because a file with this name already exists.\n";
    exit();
  }
  if (-d $opt_output) {
    print "The output directory choosen already exists. Please give me another Name.\n";
    exit();
  }
  mkdir $opt_output;

  my $file_out_name = fileparse($opt_reffile);
  
  $ostreamGFF_file    = "$opt_output/$file_out_name";
  $ostreamLog_file    = $opt_output."/error.txt";
  $ostreamReport_file = $opt_output."/report.txt";
}

my $ostreamGFF    = prepare_gffout( $ostreamGFF_file);
my $ostreamLog    = prepare_fileout($ostreamLog_file);
my $ostreamReport = prepare_fileout($ostreamReport_file);



###############################################
####### END Manage files (input output) #######
###############################################
#my $stringPrint = strftime "%m/%d/%Y at %Hh%Mm%Ss", localtime;
my $stringPrint = strftime "%m/%d/%Y", localtime;
$stringPrint .= "\nusage: $0 @copyARGV\n";
if ($opt_name) {
  $prefixName = $opt_name;
  $stringPrint .= "->IDs are changed using <$opt_name> as prefix.\nIn the case of discontinuous features (i.e. a single feature that exists over multiple genomic locations) the same ID may appear on multiple lines.".
  " All lines that share an ID collectively represent a signle feature.\n";
}
if ($opt_nameU) {
  $stringPrint .= "->IDs will be changed using <$opt_nameU> as prefix. Features that shared an ID collectively (e.g. CDS, UTRs, etc...) will now have each an uniq ID.\n";
  $prefixName = $opt_nameU;
}

# Display
$ostreamReport->print($stringPrint);
if ($opt_output) {
  print_time("$stringPrint");
} # When ostreamReport is a file we have to also display on screen

                  #          +------------------------------------------------------+
                  #          |+----------------------------------------------------+|
                  #          ||                       MAIN                         ||
                  #          |+----------------------------------------------------+|
                  #          +------------------------------------------------------+

######################
### Parse GFF input #
my ($hash_omniscient) = slurp_gff3_file_JD({ input => $opt_reffile });

### END Parse GFF input #
#########################

#Print directly what has been read
print_time("Compute statistics");
print_omniscient_statistics(
  {
    input => $hash_omniscient,
    output => $ostreamReport
  }
);

################################
# MANAGE FUNCTIONAL INPUT FILE #

#####################
# Manage Blast File #
my $db;
my %allIDs;

if (defined $opt_BlastFile) {
  # read fasta file and save info in memory
  print_time("Look at the fasta database");
  $db = Bio::DB::Fasta->new($opt_dataBase);

  # JN: Begin parse fasta
  # JN: Alternative parsing of fasta. Picking up GNs as we go.
  # Save ID in lower case to avoid cast problems
  my $dbstream = $db->get_PrimarySeq_stream;
  while (my $seqobj = $dbstream->next_seq) {
    my $display_id = $seqobj->display_id;
    my $lc_display_id = lc($display_id);
    $allIDs{$lc_display_id} = $display_id;
    my $desc = $seqobj->desc;
    if ($desc =~ /GN=(\S+)/) {
        my $GN = $1;
        my $lc_GN = lc($GN);
        $fasta_id_gn_hash{$lc_display_id} = $lc_GN;
    }
    else {
      $nbGnNotPresentInDb++;
      $fasta_id_gn_hash{$lc_display_id} = undef;
    }
  } # JN: End parse fasta

  print_time("Parsing Finished");

  # parse blast output
  print( "Reading features from $opt_BlastFile...\n");
  parse_blast($streamBlast, $opt_blastEvalue, $hash_omniscient);
}

########################
# Manage Interpro File #
if (defined $opt_InterproFile) {
  parse_interpro_tsv($streamInter, $opt_InterproFile);

  # create streamOutput
  if ($opt_output) {
    foreach my $type (keys %functionData) {
      my $ostreamFunct = IO::File->new();
      $ostreamFunct->open( $opt_output."/$type.txt", 'w' ) or
        croak(
          sprintf( "Can not open '%s' for writing %s", $opt_output."/$type.txt", $! )
        );
      $functionStreamOutput{$type} = $ostreamFunct;
    }
  }
}
# END MANAGE FUNCTIONAL INPUT FILE #
####################################

###########################
# change FUNCTIONAL information if asked for
if ($opt_BlastFile || $opt_InterproFile ) {
  print_time( "load FUNCTIONAL information" );

  #################
  # == LEVEL 1 == #
  #################
  foreach my $primary_tag_level1 (keys %{$hash_omniscient ->{'level1'}}) { # primary_tag_level1 = gene or repeat etc...
    foreach my $id_level1 (keys %{$hash_omniscient ->{'level1'}{$primary_tag_level1}}) {
      my $feature_level1 = $hash_omniscient->{'level1'}{$primary_tag_level1}{$id_level1};

      #Manage Name 
      clean_attribute($feature_level1, "Name"); # Clean NAME attribute
      if ( $opt_BlastFile ) {

        if (exists ($geneNameBlast{$id_level1})) {
          my @list_names = @{$geneNameBlast{$id_level1}};
          create_or_append_tag($feature_level1, 'Name', \@list_names);
          $nbNamedGene++;

          # Keep track of duplicated gene names <= Find another way
          foreach my $name (@list_names){

            if (exists ($geneNameGiven{$name})) {
              $nbDuplicateNameGiven++; # track total
              $duplicateNameGiven{$name}++; # track diversity
            }
            else { # first time we have given this name
              $geneNameGiven{$name}++;
            }
          }
        }
      }

      #################
      # == LEVEL 2 == #
      #################
      foreach my $primary_tag_key_level2 (keys %{$hash_omniscient->{'level2'}}) { # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...

        if ( exists_keys ($hash_omniscient, ('level2', $primary_tag_key_level2, $id_level1) ) ) {
          foreach my $feature_level2 ( @{$hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_level1}} ) {

            my $level2_ID = lc($feature_level2->_tag_value('ID'));

            # Manage Name 
            clean_attribute($feature_level2, "Name"); # Clean NAME attribute
            if ($opt_BlastFile) {
              # add gene Name
              if (exists ($mRNANameBlast{$level2_ID})) {
                create_or_append_tag($feature_level2, 'Name', $mRNANameBlast{$level2_ID});
                add_attribute_to_cds($hash_omniscient, $level2_ID, 'Name', $mRNANameBlast{$level2_ID});
              }
              # add OS attribute
              if (exists ($blast_organism{$level2_ID})) {
                create_or_append_tag($feature_level2, 'os', $blast_organism{$level2_ID});
                add_attribute_to_cds($hash_omniscient, $level2_ID, 'os', $blast_organism{$level2_ID});
              }
              # add PE attribute
              if (exists ($blast_protEvidence{$level2_ID})) {
                create_or_append_tag($feature_level2, 'pe', $blast_protEvidence{$level2_ID});
                add_attribute_to_cds($hash_omniscient, $level2_ID, 'pe', $blast_protEvidence{$level2_ID});
              }
              # add SV attribute
              if (exists ($blast_seqVersion{$level2_ID})) {
                create_or_append_tag($feature_level2, 'sv', $blast_seqVersion{$level2_ID});
                add_attribute_to_cds($hash_omniscient, $level2_ID, 'sv', $blast_seqVersion{$level2_ID});
              }

              #add UniprotID attribute
              if (exists ($mRNAUniprotIDFromBlast{$level2_ID})) {
                my $mRNAUniprotID = $mRNAUniprotIDFromBlast{$level2_ID};
                create_or_replace_tag($feature_level2, 'uniprot_id', $mRNAUniprotID);
                add_attribute_to_cds($hash_omniscient, $level2_ID, 'uniprot_id', $mRNAUniprotID);
              }

              #add evalue of the blast
              if (exists ($blast_evalue{$level2_ID})) {
                create_or_replace_tag($feature_level2, 'evalue', $blast_evalue{$level2_ID});
                add_attribute_to_cds($hash_omniscient, $level2_ID, 'evalue', $blast_evalue{$level2_ID});
              }

              # JN: Add info on existence of GN= tag in fasta header in blast db file: gn_present=yes|no|NA
              if (exists($l2_gn_present_hash{$level2_ID})) {
                my $gn_status = $l2_gn_present_hash{$level2_ID};
                if ($gn_status eq 'no' ) {
                  $nbGnNotPresentForMrna++;
                }
                create_or_replace_tag($feature_level2, 'gn_present', $gn_status) if ($opt_addGnPresentTag);
              }
              else {
                create_or_replace_tag($feature_level2, 'gn_present', 'NA') if ($opt_addGnPresentTag);
              }

              my $productData = printProductFunct($level2_ID);

              #add product attribute
              clean_attribute($feature_level2, "product"); # Clean product attribute
              if ($productData ne "") {
                add_attribute_to_cds($hash_omniscient, $level2_ID, 'product', $productData);
                if ($feature_level2->has_tag('pseudo')) {
                  create_or_replace_tag($feature_level2, 'Note', "product:$productData");
                }
                else {
                  create_or_append_tag($feature_level2, 'product', $productData);
                }
              }
              else {
                add_attribute_to_cds($hash_omniscient, $level2_ID, 'product', "hypothetical protein");
                if ($feature_level2->has_tag('pseudo')) {
                  create_or_replace_tag($feature_level2, 'Note', "product:hypothetical protein");
                }
                else {
                  create_or_append_tag($feature_level2, 'product', "hypothetical protein");
                }
              } #Case where the protein is not known
            }

            # print function if option
            if ($opt_InterproFile) {
              my $parentID = $feature_level2->_tag_value('Parent');

              if (addFunctions($hash_omniscient, $feature_level2, $opt_output)) {
                $nbmRNAwithFunction++;
                $geneWithFunction{$parentID}++;
                if (exists ($geneWithoutFunction{$parentID})) {
                  delete $geneWithoutFunction{$parentID};
                }
              }
              else {
                $nbmRNAwithoutFunction++;
                if (! exists ($geneWithFunction{$parentID})) {
                  $geneWithoutFunction{$parentID}++;
                }
              }
            }
          }
        }
      }
    }
  }
}

###########################
# change names if asked for
if ($opt_nameU || $opt_name ) { #|| $opt_BlastFile || $opt_InterproFile) {
  print_time("load new IDs");

  my %hash_sortBySeq;
  foreach my $tag_level1 ( keys %{$hash_omniscient->{'level1'}}) {
    foreach my $level1_id ( keys %{$hash_omniscient->{'level1'}{$tag_level1}}) {
      my $position = $hash_omniscient->{'level1'}{$tag_level1}{$level1_id}->seq_id;
      push (@{$hash_sortBySeq{$position}{$tag_level1}}, $hash_omniscient->{'level1'}{$tag_level1}{$level1_id});
    }
  }

  #################
  # == LEVEL 1 == #
  #################
  #Read by seqId to sort properly the output by seq ID
  foreach my $seqid (sort alphaNum keys %hash_sortBySeq) { # loop over all the feature level1

    foreach my $primary_tag_level1 (sort {$a cmp $b} keys %{$hash_sortBySeq{$seqid}}) {

      foreach my $feature_level1 ( sort {$a->start <=> $b->start} @{$hash_sortBySeq{$seqid}{$primary_tag_level1}}) {
        my $level1_ID = $feature_level1->_tag_value('ID');
        my $id_level1 = lc($level1_ID);
        my $newID_level1 = undef;
        #print_time( "Next gene $id_level1\n");

        #keep track of Maker ID
        if ($opt_BlastFile) { #In that case the name given by Maker is removed from ID and from Name. We have to keep track
          create_or_replace_tag($feature_level1, 'makerName', $level1_ID);
        }

        my $letter_tag = get_letter_tag($primary_tag_level1);

        if (! exists_keys(\%numbering, ($letter_tag))) {
          $numbering{$letter_tag} = $nbIDstart;
        }
        $newID_level1 = manageID($prefixName, $numbering{$letter_tag}, $letter_tag );
        $numbering{$letter_tag}++;
        create_or_replace_tag($feature_level1, 'ID', $newID_level1);

        $finalID{$feature_level1->_tag_value('ID')} = $newID_level1;

        #################
        # == LEVEL 2 == #
        #################
        foreach my $primary_tag_level2 (keys %{$hash_omniscient->{'level2'}}) { # primary_tag_level2 = mrna or mirna or ncrna or trna etc...

          if ( exists_keys ($hash_omniscient, ('level2', $primary_tag_level2, $id_level1) ) ) {
            foreach my $feature_level2 ( @{$hash_omniscient->{'level2'}{$primary_tag_level2}{$id_level1}}) {

              my $level2_ID = $feature_level2->_tag_value('ID');
              my $newID_level2 = undef;

              #keep track of Maker ID
              if ($opt_InterproFile) { #In that case the name given by Maker is removed from ID and from Name. We have to keep track
                create_or_replace_tag($feature_level2, 'makerName', $level2_ID);
              }

              my $letter_tag = get_letter_tag($primary_tag_level2);
              if (! exists_keys(\%numbering, ($letter_tag))) {
                $numbering{$letter_tag} = $nbIDstart;
              }
              $newID_level2 = manageID($prefixName, $numbering{$letter_tag}, $letter_tag);
              $numbering{$letter_tag}++;
              create_or_replace_tag($feature_level2, 'ID', $newID_level2);
              create_or_replace_tag($feature_level2, 'Parent', $newID_level1);

              $finalID{$level2_ID} = $newID_level2;

              #################
              # == LEVEL 3 == #
              #################
              foreach my $primary_tag_level3 (keys %{$hash_omniscient->{'level3'}}) { # primary_tag_key_level3 = cds or exon or start_codon or utr etc...
                if ( exists_keys ($hash_omniscient, ('level3', $primary_tag_level3, lc($level2_ID)) ) ) {
                  foreach my $feature_level3 ( @{$hash_omniscient->{'level3'}{$primary_tag_level3}{lc($level2_ID)}}) {
                    #keep track of Maker ID
                    my $level3_ID = $feature_level3->_tag_value('ID');
                    if ($opt_InterproFile) { #In that case the name given by Maker is removed from ID and from Name. We have to kee a track
                      create_or_replace_tag($feature_level3, 'makerName', $level3_ID);
                    }

                    my $letter_tag = get_letter_tag($primary_tag_level3);
                    if (! exists_keys(\%numbering, ($letter_tag))) {
                      $numbering{$letter_tag} = $nbIDstart;
                    }

                    my $newID_level3 = manageID($prefixName, $numbering{$letter_tag}, $letter_tag);
                    if ( $primary_tag_level3 =~ /cds/ or $primary_tag_level3 =~ /utr/ ) {
                      if ($opt_nameU) {
                        $numbering{$letter_tag}++;
                      }
                    }
                    else {
                      $numbering{$letter_tag}++;
                    }
                    create_or_replace_tag($feature_level3, 'ID', $newID_level3);
                    create_or_replace_tag($feature_level3, 'Parent', $newID_level2);

                    $finalID{$level3_ID} = $newID_level3;
                  }
                  #save the new l3 into the new l2 id name
                  $hash_omniscient->{'level3'}{$primary_tag_level3}{lc($newID_level2)} = delete $hash_omniscient->{'level3'}{$primary_tag_level3}{lc($level2_ID)} # delete command return the value before deleting it, so we just transfer the value
                }
                if ( $opt_name and ($primary_tag_level3 =~ /cds/ or $primary_tag_level3 =~ /utr/ ) ) {
                  my $letter_tag = get_letter_tag($primary_tag_level3);
                  $numbering{$letter_tag}++;
                } # with this option we increment UTR name only for each UTR (cds also)
              }
            }
            if ($newID_level1) {
              $hash_omniscient->{'level2'}{$primary_tag_level2}{lc($newID_level1)} = delete $hash_omniscient->{'level2'}{$primary_tag_level2}{$id_level1}; # modify the id key of the hash. The delete command return the value before deleting it, so we just transfer the value
            }
          }
        }

        if ($newID_level1) {
          $hash_omniscient->{'level1'}{$primary_tag_level1}{lc($newID_level1)} = delete $hash_omniscient->{'level1'}{$primary_tag_level1}{$id_level1}; # modify the id key of the hash. The delete command return the value before deleting it, so we just transfer the value
        }
      }
    }
  }
}

###########################
# RESULT PRINTING
###########################

##############################
# print FUNCTIONAL INFORMATION

# first table name\tfunction
if ($opt_output) {
  foreach my $function_type (keys %functionOutput) {
    my $streamOutput = $functionStreamOutput{$function_type};
    foreach my $ID (keys %{$functionOutput{$function_type}}) {

      if ($opt_nameU || $opt_name ) {
        print $streamOutput $finalID{$ID}."\t".$functionOutput{$function_type}{$ID}."\n";
      }
      else {
        print $streamOutput $ID."\t".$functionOutput{$function_type}{$ID}."\n";
      }
    }
  }
}

# NOW summarize
$stringPrint = ""; # reinitialise (use at the beginning)
if ($opt_InterproFile) {
  #print INFO
  my $lineB =      "_________________________________________________________________________________________________________________________________";
  $stringPrint .= " ".$lineB."\n";
  $stringPrint .= "|                         | Nb Total term           | Nb mRNA with term       | Nb mRNA updated by term | Nb gene updated by term |\n";
  $stringPrint .= "|                         | in raw File             | in raw File             | in our annotation file  | in our annotation file  |\n";
  $stringPrint .= "|".$lineB."|\n";

  foreach my $type (sort keys %functionData) {
    my $total_type = $TotalTerm{$type};
    my $mRNA_type_raw = $functionDataAdded{$type};
    my $mRNA_type = keys %{$mRNAAssociatedToTerm{$type}};
    my $gene_type = keys %{$GeneAssociatedToTerm{$type}};
    $stringPrint .= "|".sizedPrint(" $type", 25)."|".sizedPrint($total_type, 25)."|".sizedPrint($mRNA_type_raw, 25)."|".sizedPrint($mRNA_type, 25)."|".sizedPrint($gene_type, 25)."|\n|".$lineB."|\n";
  }

  #RESUME TOTAL OF FUNCTION ATTACHED
  my $listOfFunction;
  foreach my $funct (sort keys %functionData) {
    $listOfFunction .= "$funct,";
  }
  chop $listOfFunction;
  my $nbGeneWithoutFunction = keys %geneWithoutFunction;
  my $nbGeneWithFunction = keys %geneWithFunction;
  $stringPrint .= "\n".
                  "nb mRNA without Functional annotation ($listOfFunction) = $nbmRNAwithoutFunction\n".
                  "nb mRNA with Functional annotation ($listOfFunction) = $nbmRNAwithFunction\n".
                  "nb gene without Functional annotation ($listOfFunction) = $nbGeneWithoutFunction\n".
                  "nb gene with Functional annotation ($listOfFunction) = $nbGeneWithFunction\n";
}

if ($opt_BlastFile) {
  my $nbGeneDuplicated = keys %duplicateNameGiven;
  $nbDuplicateNameGiven += $nbGeneDuplicated; # Until now we have counted only name in more, now we add the original name.
  my $nbGeneNameInBlast =  keys %geneNameBlast;
  $stringPrint .= "\n$nbGeneNameInBlast gene names have been retrieved in the blast file. $nbNamedGene gene names have been successfully inferred.\n".
  "Among them there are $nbGeneDuplicated names that are shared at least per two genes for a total of $nbDuplicateNameGiven genes.\n";

  # JN: Begin summary
  # JN: Report number of entries in $opt_dataBase without GN. Note: tentative output format
  if ($opt_dataBase) {
      $stringPrint .= "\n$nbGnNotPresentInDb entries in $opt_dataBase have no GN\n";
      $stringPrint .= "\n$nbGnNotPresentForMrna mRNA entries in gff output have no GN\n";
  } # JN: End summary

  #Lets keep track the duplicated names
  if ($opt_output) {
    my $duplicatedNameOut = IO::File->new(">".$opt_output."/duplicatedNameFromBlast.txt");
    foreach my $name (sort { $duplicateNameGiven{$b} <=> $duplicateNameGiven{$a} } keys %duplicateNameGiven) {
      print $duplicatedNameOut "$name\t".($duplicateNameGiven{$name} + 1)."\n";
    }
  }
}

if ($opt_name or $opt_nameU) {
  $stringPrint .= "\nList of Letter use to create the uniq ID:\n";
  foreach my $tag (keys %tag_hash) {
    $stringPrint .= "$tag => $tag_hash{$tag}\n";
  }
  $stringPrint .= "\n";
}

# Display
$ostreamReport->print("$stringPrint");

####################
# PRINT IN FILES
####################
print_time("Writing result...");
print_omniscient( {omniscient => $hash_omniscient, output => $ostreamGFF} );
print_time("End of script.");

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

# remove the attribute provided
sub clean_attribute {
  my ($feature, $tag) = @_;
       
  if ($opt_CleanNameAttribute and $tag eq "Name"){
    if ($feature->has_tag('Name')) {
      $feature->remove_tag('Name');
    }
  }
  if ($opt_CleanProductAttribute and $tag eq "product"){
    if ($feature->has_tag('product')) {
      $feature->remove_tag('product');
    }
  }
  if ($opt_CleanDbxrefAttribute and $tag eq "Dbxref"){
    if ($feature->has_tag('Dbxref')) {
      $feature->remove_tag('Dbxref');
    }
  }
  if ($opt_CleanOntology_termAttribute and $tag eq "Ontology_term"){
    if ($feature->has_tag('Ontology_term')) {
      $feature->remove_tag('Ontology_term');
    }
  }
}

sub add_attribute_to_cds {
  my ($hash_omniscient, $level2_ID, $tag, $value) = @_;

  if($opt_populate_cds){
    if ( exists_keys ($hash_omniscient, ('level3', 'cds', lc($level2_ID)) ) ) {
      foreach my $feature_level3 ( @{$hash_omniscient->{'level3'}{'cds'}{lc($level2_ID)}}) {
        clean_attribute($feature_level3, $tag);
        $feature_level3->add_tag_value($tag, $value);
      }
    }
  }
}

#create or take the unique letter TAG
sub get_letter_tag {
  my ($tag) = @_;

  $tag = lc($tag);
  if (! exists_keys (\%tag_hash, ( $tag ))) {

    my $substringLength = 1;
    my $letter = uc(substr($tag, 0, $substringLength));

    while ( grep( /^\Q$letter\E$/, @tag_list) ) { # to avoid duplicate
      $substringLength++;
      $letter = uc(substr($tag, 0, $substringLength));
    }
    $tag_hash{ $tag } = uc($letter);
    push(@tag_list, $letter)
  }
  return $tag_hash{ $tag };
}

# creates gene ID correctly formated (PREFIX,TYPE,NUMBER) like HOMSAPG00000000001 for a Homo sapiens gene.
sub manageID {
  my ($prefix, $nbName, $type) = @_;
  my $result = "";
  my $numberNum = 11;
  my $GoodNum = "";
  for (my $i=0; $i<$numberNum-length($nbName); $i++) {
    $GoodNum .= "0";
  }
  $GoodNum .= $nbName;
  $result = "$prefix$type$GoodNum";

  return $result;
}

# Create String containing the product information associated to the mRNA
sub printProductFunct {
  my ($refname) = @_;
  my $String = "";
  my $first = "yes";
  if (exists $mRNAproduct{$refname}) {
    foreach my $element (sort @{$mRNAproduct{$refname}}) {
      if ($first eq "yes") {
        $String .= "$element";
        $first = "no";
      }
      else {
        $String .= ",$element";
      }
    }
  }
  return $String;
}

sub addFunctions {
  my ($hash_omniscient, $feature, $opt_output) = @_;

  my $functionAdded = undef;
  my $ID = lc($feature->_tag_value('ID'));
  foreach my $function_type ( sort keys %functionData) {

    if (exists ($functionData{$function_type}{$ID})) {
      $functionAdded = "true";

      my $data_list;

      if (lc($function_type) eq "go") {
        clean_attribute($feature, "Ontology_term"); # Clean Ontology_term attribute
        foreach my $data (@{$functionData{$function_type}{$ID}}) {
          $feature->add_tag_value('Ontology_term', $data);
          $data_list .= "$data,";
          $functionDataAdded{$function_type}++;
          add_attribute_to_cds($hash_omniscient, $ID, 'Ontology_term', $data);
        }
      }
      else {
        clean_attribute($feature, "Dbxref"); # Clean Dbxref attribute
        foreach my $data (@{$functionData{$function_type}{$ID}}) {
          $feature->add_tag_value('Dbxref', $data);
          $data_list .= "$data,";
          $functionDataAdded{$function_type}++;
          add_attribute_to_cds($hash_omniscient, $ID, 'Dbxref', $data);
        }
      }

      if ($opt_output) {
          my $ID = $feature->_tag_value('ID');
          chop $data_list;
          $functionOutput{$function_type}{$ID} = $data_list;
        }
    }
  }
  return $functionAdded;
}

# method to parse blast file
sub parse_blast {
  my($file_in, $opt_blastEvalue, $hash_omniscient) = @_;

#################################################################################
####### Step 1 : CATCH all candidates (the better candidate for each mRNA)####### (with a gene name)

  my %candidates;
  my %gene_name_HoH = (); # JN: key: maker-Bi03_p1mp_000319F-est_gff_StringTie-gene-5.8-mRNA-1, value: {'genename' => 1}

  while(my $line = <$file_in>) {
    my @values = split(/\t/, $line);
    my $l2_name = lc($values[0]);  # JN: maker-Bi03_p1mp_000319F-est_gff_StringTie-gene-5.8-mRNA-1
    my $prot_name = $values[1]; # JN: sp|Q4FZT2|PPME1_RAT
    my @prot_name_sliced = split(/\|/, $values[1]);
    my $uniprot_id = $prot_name_sliced[1];
    #print "uniprot_id: ".$uniprot_id."\n" if ($opt_verbose);
    my $evalue = $values[10];
    #print "Evalue: ".$evalue."\n" if ($opt_verbose);

    #if does not exist fill it if over the minimum evalue
    if (! exists_keys(\%candidates, ($l2_name)) or @{$candidates{$l2_name}} > 3 ) { # the second one means we saved an error message as candidates we still have to try to find a proper one

      if ( $evalue <= $opt_blastEvalue ) {
        # JN: begin gene_name_Debug HoH
        my $lc_prot_name = lc($prot_name);
        if (exists($fasta_id_gn_hash{$lc_prot_name})) {    # JN: Key exists if gene name or undef
          if (defined($fasta_id_gn_hash{$lc_prot_name})) { # JN: Only defined if gene name
            my $gn = $fasta_id_gn_hash{$lc_prot_name};     # JN: Get the gene name
            $gene_name_HoH{$l2_name}{$gn}++;               # JN: Count the gene name
          }
          else {                                           # JN: If not defined, the 'GN=' is missing
            undef($gene_name_HoH{$l2_name});
          }
        }
        # JN: End Debug gene_name_HoH

        my $protID_correct = undef;

        if ( exists $allIDs{lc($prot_name)}) {
          $protID_correct = $allIDs{lc($prot_name)};
          my $header = $db->header( $protID_correct );

          if (! $header =~ m/GN=/) {
            $ostreamLog->print( "No gene name (GN=) in this header $header\n") if ($opt_verbose or $opt_output);
            $candidates{$l2_name} = ["error", $evalue, $prot_name."-".$l2_name];
          }

          if ($header =~ /PE=([1-5])/) {
            if ($1 <= $opt_pe) {
              $candidates{$l2_name} = [$header, $evalue, $uniprot_id];
            }
          }
          else {
            $ostreamLog->print("No Protein Existence (PE) information in this header: $header\n") if ($opt_verbose or $opt_output);
          }
        }
        else {
          $ostreamLog->print( "ERROR $prot_name not found among the db! You probably didn't give to me the same fasta file than the one used for the blast. (l2=$l2_name)\n" ) if ($opt_verbose or $opt_output);
          $candidates{$l2_name} = ["error", $evalue, $prot_name."-".$l2_name];
        }
      }
    }
    elsif ( $evalue < $candidates{$l2_name}[1] ) { # better evalue for this record

      # JN: begin gene_name_Debug HoH
      my $lc_prot_name = lc($prot_name);
      if (exists($fasta_id_gn_hash{$lc_prot_name})) {    # JN: Key exists if gene name or undef
        if (defined($fasta_id_gn_hash{$lc_prot_name})) { # JN: Only defined if gene name
          my $gn = $fasta_id_gn_hash{$lc_prot_name};     # JN: Get the gene name
          $gene_name_HoH{$l2_name}{$gn}++;               # JN: Count the gene name
        }
        else {                                           # JN: If not defined, the 'GN=' is missing
          undef($gene_name_HoH{$l2_name});
        }
      }
      # JN: End Debug gene_name_HoH

      my $protID_correct = undef;

      if ( exists $allIDs{lc($prot_name)}) {

        $protID_correct = $allIDs{lc($prot_name)};
        my $header = $db->header( $protID_correct );

        if (! $header =~ m/GN=/) {
          $ostreamLog->print("No gene name (GN=) in this header $header\n") if ($opt_verbose or $opt_output);
        }

        if ($header =~ /PE=([1-5])/) {
          if ($1 <= $opt_pe) {
            $candidates{$l2_name} = [$header, $evalue, $uniprot_id];
          }
        }
        else {
          $ostreamLog->print( "No Protein Existence (PE) information in this header: $header\n") if ($opt_verbose or $opt_output);
        }
      }
      else {
        $ostreamLog->print( "ERROR $prot_name not found among the db! You probably didn't give to me the same fasta file than the one used for the blast. (l2=$l2_name)\n") if ($opt_verbose or $opt_output);
      }
    }
  }

  my $nb_desc = keys %candidates;
  $ostreamLog->print( "We have $nb_desc description candidates.\n") if ($opt_verbose or $opt_output);

  ##################################################
  ####### Step 2 : go through all candidates ####### report gene name for each mRNA
  my %name_checker;
  foreach my $l2 (sort keys %candidates) {
    # JN: Here we need to not(?) return error above to be able to differentiate the cases without GN?
    if ( $candidates{$l2}[0] eq "error" ) {
      $ostreamLog->print("error nothing found for $candidates{$l2}[2]\n") if ($opt_verbose or $opt_output);
      next;
    }

    #Save uniprot id of the best match
   
    $mRNAUniprotIDFromBlast{$l2} = $candidates{$l2}[2];
    print "save protein ID for $l2 : ".$candidates{$l2}[2]."\n" if ($opt_verbose);
    
    #Save evalu
    $blast_evalue{$l2} = $candidates{$l2}[1];
    print "save blast evalue for $l2 : ".$candidates{$l2}[1]."\n" if ($opt_verbose);

    # Parse header
    my $header = $candidates{$l2}[0];
    print "header: ".$header."\n" if ($opt_verbose);

    if ($header =~ m/(^[^\s]+)\s(.+?(?= \w{2}=))(.+)/) {
      my $protID = $1;
      my $description = $2;
      my $theRest = $3;
      $theRest =~ s/\n//g;
      $theRest =~ s/\r//g;
      my $nameGene = undef;
      print "description: ".$description."\n" if ($opt_verbose);
      push ( @{ $mRNAproduct{$l2} }, $description );

      #deal with the rest
      my %hash_rest;
      my $tuple = undef;
      while ($theRest) {
        ($theRest, $tuple) = stringCatcher($theRest);
        my ($type, $value) = split /=/, $tuple;
        print "$protID: type:$type --- value:$value\n" if ($opt_verbose);
        $hash_rest{lc($type)} = $value;
      }

      if (exists($hash_rest{"gn"})) {
        $nameGene = $hash_rest{"gn"};

        if (exists_keys($hash_omniscient, ('other', 'l2tol1', $l2) ) ){
          # Save mRNA name into mRNA features
          $mRNANameBlast{$l2} = $nameGene;

          my $geneID = $hash_omniscient->{'other'}{'l2tol1'}{$l2};

          # Save mRNA names into gene features (a list because we can have several gene names if mRNA isoforms were refering to different gene names)
          if (! exists_keys(\%name_checker,(lc($geneID),lc($nameGene) ))){
            push ( @{ $geneNameBlast{lc($geneID)} }, $nameGene );
          }
          $name_checker{lc($geneID)}{lc($nameGene)}++;

        }
        else {
          $ostreamLog->print( "No parent found for $l2 (defined in the blast file) in hash_omniscient{l2tol1} (created by the gff file).\n") if ($opt_verbose or $opt_output);
        }
      }
      else {
        $ostreamLog->print( "No gene name (GN) tag found in the header: $header\n") if ($opt_verbose or $opt_output);
      }

      # catch organism
      if (exists($hash_rest{"os"})) {
        $blast_organism{$l2} = $hash_rest{"os"};
      }

      # catch Protein evidence
      if (exists($hash_rest{"pe"})) {
        $blast_protEvidence{$l2} = $hash_rest{"pe"};
      }

      # catch Sequence Version
      if (exists($hash_rest{"sv"})) {
        $blast_seqVersion{$l2} = $hash_rest{"sv"};
      }
    }
    else {
      $ostreamLog->print( "Header from the db fasta file doesn't match the regular expression: $header\n") if ($opt_verbose or $opt_output);
    }
  }

  # JN: Begin traversing gene_name_HoH, and populate global hash l2_gn_present_hash
  while ( my ($l2_key, $values) = each %gene_name_HoH ) {
    my $size = 0;
    if (defined($values)) { # JN: If defined, we have at least one GN
      $size = scalar(%{$values});
      $l2_gn_present_hash{$l2_key} = "yes"; # JN: gn_present=yes
      # JN: TODO: need to check and handle(?) cases where we have several different hits
      if ($size > 1) {
        my (@vals) = keys (%{$values});
        #$ostreamLog->print( "DEBUG JN: level 2 label \'$l2_key\' have several GN values: @vals\n") if ($opt_verbose or $opt_output);
        $ostreamLog->print( "DEBUG JN: level 2 label \'$l2_key\' have several GN values: @vals\n") if ($DEBUG); # JN: Debug printing
      }
    }
    else {
      $l2_gn_present_hash{$l2_key} = "no"; # JN: gn_present=no
    }
  }
  # JN: End traversing gene_name_HoH
}

#uniprotHeader string splitter
sub stringCatcher {
  my($String) = @_;
  my $newString = undef;

  if ( $String =~ m/(\w{2}=.+?(?= \w{2}=))(.+)/ ) {
    $newString = substr $String, length($1) + 1;
    return ($newString, $1);
  }
  else {
    $String =~ s/^\s+|\s+$//g; # Removes leading and trailing spaces
    return (undef, $String);
  }
}

# method to parse Interpro file
sub parse_interpro_tsv {
  my($file_in, $fileName) = @_;
  print("Reading features from $fileName...\n");

  while( my $line = <$file_in>) {

    my @values = split(/\t/, $line);
    my $sizeList = @values;
    my $mRNAID = lc($values[0]);

    #Check for the specific DB
    my $db_name = $values[3];
    my $db_value = $values[4];
    my $db_tuple = $db_name.":".$db_value;
    print "Specific dB: ".$db_tuple."\n" if ($opt_verbose);

    if (! grep( /^\Q$db_tuple\E$/, @{$functionData{$db_name}{$mRNAID}} )) { #to avoid duplicate
      $TotalTerm{$db_name}++;
      push ( @{$functionData{$db_name}{$mRNAID}}, $db_tuple );
      if (exists_keys($hash_omniscient, ('other', 'l2tol1', $mRNAID) ) ){ ## check if exists among our current gff annotation file analyzed
        $mRNAAssociatedToTerm{$db_name}{$mRNAID}++;
        $GeneAssociatedToTerm{$db_name}{$hash_omniscient->{'other'}{'l2tol1'}{$mRNAID}}++;
      }
    }

    #check for interpro
    if ( $sizeList > 11 ) {
      my $db_name = "InterPro";
      my $interpro_value = $values[11];
      $interpro_value =~ s/\n//g;
      my $interpro_tuple = "InterPro:".$interpro_value;
      print "interpro dB: ".$interpro_tuple."\n" if ($opt_verbose);
      next if $interpro_value eq "-"; #fix 147

      if (! grep( /^\Q$interpro_tuple\E$/, @{$functionData{$db_name}{$mRNAID}} )) { #to avoid duplicate
        $TotalTerm{$db_name}++;
        push ( @{$functionData{$db_name}{$mRNAID}}, $interpro_tuple );
        if (exists_keys($hash_omniscient, ('other', 'l2tol1', $mRNAID) ) ){ ## check if exists among our current gff annotation file analyzed
          $mRNAAssociatedToTerm{$db_name}{$mRNAID}++;
          $GeneAssociatedToTerm{$db_name}{$hash_omniscient->{'other'}{'l2tol1'}{$mRNAID}}++;
        }
      }
    }

    #check for GO
    if ( $sizeList > 13 ) {
      my $db_name = "GO";
      my $go_flat_list = $values[13];
      $go_flat_list =~ s/\n//g;
      my @go_list = split(/\|/, $go_flat_list); #cut at character |
      foreach my $go_tuple (@go_list) {
        print "GO term: ".$go_tuple."\n" if ($opt_verbose);
        next if $go_tuple eq "-"; #fix kira
        
        if (! grep( /^\Q$go_tuple\E$/, @{$functionData{$db_name}{$mRNAID}} )) { #to avoid duplicate
          $TotalTerm{$db_name}++;
          push ( @{$functionData{$db_name}{$mRNAID}}, $go_tuple );
          if (exists_keys($hash_omniscient, ('other', 'l2tol1', $mRNAID) ) ){ ## check if exists among our current gff annotation file analyzed
            $mRNAAssociatedToTerm{$db_name}{$mRNAID}++;
            $GeneAssociatedToTerm{$db_name}{$hash_omniscient->{'other'}{'l2tol1'}{$mRNAID}}++;
          }
        }
      }
    }

    #check for pathway
    if ( $sizeList > 14 ) {
      my $pathway_flat_list = $values[14];
      $pathway_flat_list =~ s/\n//g;
      $pathway_flat_list =~ s/ //g;
      my @pathway_list = split(/\|/, $pathway_flat_list); #cut at character |
      foreach my $pathway_tuple (@pathway_list) {
        my @tuple = split(/:/, $pathway_tuple); #cut at character :
        my $db_name = $tuple[0];
        print "pathway info: ".$pathway_tuple."\n" if ($opt_verbose);
        next if ($pathway_tuple eq "-"); # avoid empty pathway tuple
        if (! grep( /^\Q$pathway_tuple\E$/, @{$functionData{$db_name}{$mRNAID}} ) ) { # to avoid duplicate
          $TotalTerm{$db_name}++;
          push ( @{$functionData{$db_name}{$mRNAID}} , $pathway_tuple );
          if (exists_keys($hash_omniscient, ('other', 'l2tol1', $mRNAID) ) ){ ## check if exists among our current gff annotation file analyzed
            $mRNAAssociatedToTerm{$db_name}{$mRNAID}++;
            $GeneAssociatedToTerm{$db_name}{$hash_omniscient->{'other'}{'l2tol1'}{$mRNAID}}++;
          }
        }
      }
    }
  }
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

agat_sp_manage_functional_annotation.pl

=head1 DESCRIPTION

The script take a gff3 file as input and blast and/or interpro output in order
to attach functional annotation to corresponding features within the gff file.

>The blast against Protein Database (outfmt 6) allows to fill the field/attribute
NAME for gene and PRODUCT for mRNA.

>The Interpro result (.tsv) file allows to fill the DBXREF field/attribute with
pfam, tigr, interpro, GO, KEGG, etc... terms data.

With the <id> option the script will change all the ID field by an Uniq ID
created from the given prefix, a letter to specify the kind of feature (G,T,C,E,U),
and the feature number.

The result is written to the specified output file, or to STDOUT.

About the TSV format from interproscan:
=======================================

The TSV format presents the match data in columns as follows:

 1.Protein Accession (e.g. P51587)
 2.Sequence MD5 digest (e.g. 14086411a2cdf1c4cba63020e1622579)
 3.Sequence Length (e.g. 3418)
 4.Analysis (e.g. Pfam / PRINTS / Gene3D)
 5.Signature Accession (e.g. PF09103 / G3DSA:2.40.50.140)
 6.Signature Description (e.g. BRCA2 repeat profile)
 7.Start location
 8.Stop location
 9.Score - is the e-value (or score) of the match reported by member database method (e.g. 3.1E-52)
 10.Status - is the status of the match (T: true)
 11.Date - is the date of the run
 12.(InterPro annotations - accession (e.g. IPR002093) - optional column; only displayed if -iprlookup option is switched on)
 13.(InterPro annotations - description (e.g. BRCA2 repeat) - optional column; only displayed if -iprlookup option is switched on)
 14.(GO annotations (e.g. GO:0005515) - optional column; only displayed if --goterms option is switched on)
 15.(Pathways annotations (e.g. REACT_71) - optional column; only displayed if --pathways option is switched on)

P.S: The 9th column contains most of time e-value, but can contain also score (e.g Prosite). To understand the difference: https://myhits.isb-sib.ch/cgi-bin/help?doc=scores.html

About the outfmt 6 from blast:
==============================

 1.  qseqid  query (e.g., gene) sequence id
 2.  sseqid  subject (e.g., reference genome) sequence id
 3.  pident  percentage of identical matches
 4.  length  alignment length
 5.  mismatch  number of mismatches
 6.  gapopen   number of gap openings
 7.  qstart  start of alignment in query
 8.  qend  end of alignment in query
 9.  sstart  start of alignment in subject
 10.   send  end of alignment in subject
 11.   evalue  expect value
 12.   bitscore  bit score

Currently the best e-value win... That means another hit with a lower e-value
(but still over the defined threshold anyway) even if it has a better PE value
will not be reported.

=head1 SYNOPSIS

    agat_sp_manage_functional_annotation.pl -f infile.gff [-b blast_infile][-d uniprot.fasta][-i interpro_infile.tsv][--id ABCDEF][-a][-o output]
    agat_sp_manage_functional_annotation.pl --help

=head1 OPTIONS

=over 8

=item B<-f>, B<--reffile>,B<-ref> , B<--gff> or B<--gff3>

String - Input GTF/GFF file.

=item B<-b> or B<--blast>

String - Input blast ( outfmt 6 = tabular ) usually made by blasting the proteins resulting from the GFF/GTF file provided as input
and a confident protein database (e.g. Swissprot/Uniprot). The file makse a bridge between the feature ID from the GFF/GTF and the 
best protein ID matched in the used database. Thanks to that link the Name and products (sometimes called descriptions) information 
will be extracted from the database fasta file and added in the GFF file. You must provide the same database via --db as the one used 
to create this blast output file.

=item B<--clean_name> 

Bolean - When activated, if the Name attribute already exists, it we be cleaned. Otherwise Name retrieved by --blast + --db options 
will be appended. Default False (Name attribute not cleaned).

=item B<--clean_product> 

Bolean - When activated, if the product attribute already exists, it we be cleaned. Otherwise product retrieved by --blast + --db options 
will be appended. Default False (product attribute not cleaned).

=item B<--clean_dbxref> 

Bolean - When activated, if the Dbxref attribute already exists, it we be cleaned. Otherwise Dbxref retrieved by --interpro option 
will be appended. Default False (Dbxref attribute not cleaned).

=item B<--clean_ontology> 

Bolean - When activated, if the Ontology_term attribute already exists, it we be cleaned. Otherwise Ontology_term retrieved by --interpro option 
will be appended. Default False (Ontology_term attribute not cleaned).

=item B<-d> or B<--db>

String - The fasta file that has been used as DB for the blast. Gene names and products  (sometimes called descriptions) will be fished from this file.

=item B<--be> or B<--blast_evalue>

Float - Maximum e-value to keep the annotation from the blast file. By default 1e-6.

=item B<--pe>

Integer - The PE (protein existence) in the uniprot header indicates the type of evidence that supports the existence of the protein.
You can decide until which protein existence level you want to consider to lift the functional information. Default 5.

1. Experimental evidence at protein level
2. Experimental evidence at transcript level
3. Protein inferred from homology
4. Protein predicted
5. Protein uncertain

=item B<-i> or B<--interpro>

String - Input interpro file (.tsv) that will be used to complement the features read from
the first file (specified with B<--ref>).

=item B<-id>

String - This option will changed the id name. It will create from id prefix (usually 6 letters) given as input, uniq IDs like prefixE00000000001. Where E mean exon. Instead E we can have C for CDS, G for gene, T for mRNA, U for Utr.
In the case of discontinuous features (i.e. a single feature that exists over multiple genomic locations) the same ID may appear on multiple lines. All lines that share an ID collectively represent a signle feature.

=item B<-idau>

Boolean - This option (id all uniq) is similar to -id option but Id of features that share an ID collectively will be change by different and uniq ID.

=item B<-nb>

Integer - Usefull only if -id is used.
This option is used to define the number that will be used to begin the numbering. By default begin by 1.

=item B<-a> or B<--addgntag>

Add information in ouptut gff about if gene-name tag ('GN=') is present in blast db fasta ('gn_present=yes')
or not ('gn_present=no'). Blast hits without an entry in the blast db will receive 'gn_present=NA'.

=item B<-o> or B<--output>

String - Output folder name with summary files. If no output file is specified, the output will be
written to STDOUT.

=item B<--pcds>

Boolean - pcds stands for populate cds. It copies the Name, product, Ontology_term, Dbxref and uniprot_id attributes from mRNA to the CDS.

=item B<-v>

Boolean - Verbose, for debug purpose.

=item B<-thread>, B<threads>, B<cpu>, B<cpus>, B<core>, B<cores>, B<job> or B<jobs>

Integer - Number of parallel processes to use for file input parsing (via forking).

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
