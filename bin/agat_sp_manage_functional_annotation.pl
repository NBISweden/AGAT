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
use Bio::Tools::GFF;
use AGAT::Omniscient;

use Data::Dumper; # JN: dedug
my $DEBUG = 1; # JN: debug

my $header = get_agat_header();
# PARAMETERS - OPTION
my $opt_reffile;
my $opt_output;
my $opt_BlastFile;
my $opt_InterproFile;
my $opt_name = undef;
my $opt_nameU;
my $opt_verbose = undef;
my $opt_help = 0;
my $opt_blastEvalue = 1e-6;
my $opt_dataBase = undef;
my $opt_pe = 5;
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
my %mRNAUniprotIDFromBlast;
my %mRNAproduct;
my %geneNameGiven;
my %duplicateNameGiven;
my $nbDuplicateNameGiven = 0;
my $nbDuplicateName = 0;
my $nbNamedGene = 0;
my $nbGeneNameInBlast = 0;
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
 'd|db=s'                   => \$opt_dataBase,
 'be|blast_evalue=i'        => \$opt_blastEvalue,
 'pe=i'                     => \$opt_pe,
 'i|interpro=s'             => \$opt_InterproFile,
 'id=s'                     => \$opt_name,
 'idau=s'                   => \$opt_nameU,
 'nb=i'                     => \$nbIDstart,
 'o|output=s'               => \$opt_output,
 'v'                        => \$opt_verbose,
 'h|help!'                  => \$opt_help
)
or pod2usage( {
  -message => 'Failed to parse command line',
  -verbose => 1,
  -exitval => 1
});

# Print Help and exit
if ($opt_help) {
  pod2usage( { -verbose => 99,
               -exitval => 0,
               -message => "$header\n" } );
}

if ( ! (defined($opt_reffile)) ) {
  pod2usage( {
         -message => "$header\nAt least 1 parameter is mandatory:\nInput reference gff file (--f)\n\n".
         "Many optional parameters are available. Look at the help documentation to know more.\n",
         -verbose => 0,
         -exitval => 1 } );
}

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
# JN: example file maker_evidence_appendedByAbinitio_blast.out
if (defined $opt_BlastFile) {
  if (! $opt_dataBase) {
    print "To use the blast output we also need the fasta of the database used for the blast (--db)\n";
    exit;
  }
  $streamBlast->open( $opt_BlastFile, 'r' ) or
    croak( sprintf( "Can not open '%s' for reading: %s", $opt_BlastFile, $! ) );
}

# Manage Interpro file
# JN: example file maker_evidence_appendedByAbinitio_interpro.tsv
if (defined $opt_InterproFile) {
  $streamInter->open( $opt_InterproFile, 'r' ) or
    croak( sprintf( "Can not open '%s' for reading: %s", $opt_InterproFile, $! ) );
}

##########################
##### Manage Output ######
my $ostreamReport;
my $ostreamGFF;
my $ostreamLog;
if (defined($opt_output)) {
  if (-f $opt_output) {
    print "Cannot create a directory with the name $opt_output because a file with this name already exists.\n";
    exit();
  }
  if (-d $opt_output) {
    print "The output directory choosen already exists. Please geve me another Name.\n";
    exit();
  }
  mkdir $opt_output;

  $ostreamReport = IO::File->new(">".$opt_output."/report.txt") or
    croak( sprintf( "Can not open '%s' for writing %s", $opt_output."/report.txt", $! ));

  my $file_out_name = fileparse($opt_reffile);
  $ostreamGFF = Bio::Tools::GFF->new(-file => ">$opt_output/$file_out_name", -gff_version => 3 ) or
    croak(sprintf( "Can not open '%s' for writing %s", $opt_output."/".$opt_reffile, $! ));

  $ostreamLog = IO::File->new(">".$opt_output."/error.txt") or
    croak( sprintf( "Can not open '%s' for writing %s", $opt_output."/log.txt", $! ));
}
else {
  $ostreamReport = \*STDOUT or die ( sprintf( "Can not open '%s' for writing %s", "STDOUT", $! ));
  $ostreamLog = \*STDOUT or die ( sprintf( "Can not open '%s' for writing %s", "STDOUT", $! ));
  $ostreamGFF = Bio::Tools::GFF->new( -fh => \*STDOUT, -gff_version => 3) or croak( sprintf( "Can not open STDOUT for writing: %s", $! ) );
}

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
my ($hash_omniscient, $hash_mRNAGeneLink) = slurp_gff3_file_JD({ input => $opt_reffile});
print_time("Parsing Finished");
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
my %fasta_id_gn_hash = (); # JN: key: lower case display_id, value: lower case GN
my $missing_gn_in_fasta_counter = 0; # JN: Count entries in db with no GN

if (defined $opt_BlastFile) {
  # read fasta file and save info in memory
  print ("look at the fasta database\n");
  $db = Bio::DB::Fasta->new($opt_dataBase); # JN: example file uniprot_sprot.fasta
  # save ID in lower case to avoid cast problems
  #my @ids = $db->get_all_primary_ids; # JN: here we parse the fasta file
  #foreach my $id (@ids) {
  #  $allIDs{lc($id)} = $id;
  #}
  # JN: Alternative parsing of fasta. Picking up GNs as we go
  my $dbstream = $db->get_PrimarySeq_stream;
  while (my $seqobj = $dbstream->next_seq) {
    my $display_id = $seqobj->display_id;
    my $lc_display_id = lc($display_id);
    $allIDs{$lc_display_id} = $display_id; # JN: example in %allIDs: 'sp|a0le17|rbfa_magmm' => 'sp|A0LE17|RBFA_MAGMM'
    my $desc = $seqobj->desc;
    if ($desc =~ /GN=(\S+)/) {
        my $GN = $1;
        my $lc_GN = lc($GN);
        $fasta_id_gn_hash{$lc_display_id} = $lc_GN;
    }
    else {
      $missing_gn_in_fasta_counter++;
      $fasta_id_gn_hash{$lc_display_id} = "missing_gn";
    }
  }
  print_time("Parsing Finished\n\n");

  # parse blast output
  # JN: example file maker_evidence_appendedByAbinitio_blast.out
  print( "Reading features from $opt_BlastFile...\n");
  parse_blast($streamBlast, $opt_blastEvalue, $hash_mRNAGeneLink); # JN: Need to see how we parse here
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
  print_time( "load FUNCTIONAL information\n" );

  #################
  # == LEVEL 1 == #
  #################
  my $missing_name_counter = 0; # JN: DEBUG
  my $missing_name_in_blast_counter = 0; # JN: DEBUG

  foreach my $primary_tag_level1 (keys %{$hash_omniscient ->{'level1'}}) { # primary_tag_level1 = gene or repeat etc...
    foreach my $id_level1 (keys %{$hash_omniscient ->{'level1'}{$primary_tag_level1}}) {

      my $feature_level1 = $hash_omniscient->{'level1'}{$primary_tag_level1}{$id_level1};

      # Clean NAME attribute
      # JN: Why do we need to remove the tag?
      # JN: Note: all entries have a Name tag for the debug example I'm using
      if ($feature_level1->has_tag('Name')) {
        $feature_level1->remove_tag('Name');
      }

      #Manage Name if option setting
      if ( $opt_BlastFile ) {
        if (exists ($geneNameBlast{$id_level1})) { # JN: Does the Name exists in the geneNameBlast hash? If not, no name is stored!
          create_or_replace_tag($feature_level1, 'Name', $geneNameBlast{$id_level1});
          $nbNamedGene++;

          # Check name duplicated given
          my $nameClean = $geneNameBlast{$id_level1};
          $nameClean =~ s/_([2-9]{1}[0-9]*|[0-9]{2,})*$//;

          my $nameToCompare;
          if (exists ($nameBlast{$nameClean})) { # We check that is really a name where we added the suffix _1
            $nameToCompare = $nameClean;
          }
          else {
            $nameToCompare = $geneNameBlast{$id_level1};
          } # it was already a gene_name like BLABLA_12

          if (exists ($geneNameGiven{$nameToCompare})) {
            $nbDuplicateNameGiven++; # track total
            $duplicateNameGiven{$nameToCompare}++; # track diversity
          }
          else {
            $geneNameGiven{$nameToCompare}++;
          } # first time we have given this name
        }
        else { # JN: Start DEBUG
          create_or_replace_tag($feature_level1, 'Name', 'DEBUG_unnamed_gene_level1'); # JN: Debug output
        } # End DEBUG
      }

      #################
      # == LEVEL 2 == #
      #################
      foreach my $primary_tag_key_level2 (keys %{$hash_omniscient->{'level2'}}) { # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...

        if ( exists_keys ($hash_omniscient, ('level2', $primary_tag_key_level2, $id_level1) ) ) {
          foreach my $feature_level2 ( @{$hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_level1}} ) {

            my $level2_ID = lc($feature_level2->_tag_value('ID'));
            # Clean NAME attribute
            if ($feature_level2->has_tag('Name')) {
              $feature_level2->remove_tag('Name');
            }

            #Manage Name if option set
            if ($opt_BlastFile) {
              # add gene Name
              if (exists ($mRNANameBlast{$level2_ID})) {
                my $mRNABlastName = $mRNANameBlast{$level2_ID};
                create_or_replace_tag($feature_level2, 'Name', $mRNABlastName);
              }
              else { # JN: Start DEBUG
                create_or_replace_tag($feature_level2, 'Name', 'DEBUG_unnamed_gene_level2'); # JN: Debug output
              } # End DEBUG

              my $productData = printProductFunct($level2_ID);

              #add UniprotID attribute
              if (exists ($mRNAUniprotIDFromBlast{$level2_ID})) {
                my $mRNAUniprotID = $mRNAUniprotIDFromBlast{$level2_ID};
                create_or_replace_tag($feature_level2, 'uniprot_id', $mRNAUniprotID);
              }

              #add product attribute
              if ($productData ne "") {
                if ($feature_level2->has_tag('pseudo')) {
                  create_or_replace_tag($feature_level2, 'Note', "product:$productData");
                }
                else {
                  create_or_replace_tag($feature_level2, 'product', $productData);
                }
              }
              else {
                if ($feature_level2->has_tag('pseudo')) {
                  create_or_replace_tag($feature_level2, 'Note', "product:hypothetical protein");
                }
                else {
                  create_or_replace_tag($feature_level2, 'product', "hypothetical protein"); # JN: check Luciles case here?
                }
              } #Case where the protein is not known
            }

            # print function if option
            if ($opt_InterproFile) {
              my $parentID = $feature_level2->_tag_value('Parent');

              if (addFunctions($feature_level2, $opt_output)) {
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
  $stringPrint .= "|                         | in raw File             |   in raw File           | in our annotation file  | in our annotation file  |\n";
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
  $stringPrint .= "nb mRNA without Functional annotation ($listOfFunction) = $nbmRNAwithoutFunction\n".
                  "nb mRNA with Functional annotation ($listOfFunction) = $nbmRNAwithFunction\n".
                  "nb gene without Functional annotation ($listOfFunction) = $nbGeneWithoutFunction\n".
                  "nb gene with Functional annotation ($listOfFunction) = $nbGeneWithFunction\n";
}

if ($opt_BlastFile) {
  my $nbGeneDuplicated = keys %duplicateNameGiven;
  $nbDuplicateNameGiven = $nbDuplicateNameGiven + $nbGeneDuplicated; # Until now we have counted only name in more, now we add the original name.
  $stringPrint .= "$nbGeneNameInBlast gene names have been retrieved in the blast file. $nbNamedGene gene names have been successfully inferred.\n".
  "Among them there are $nbGeneDuplicated names that are shared at least per two genes for a total of $nbDuplicateNameGiven genes.\n";
  # "We have $nbDuplicateName gene names duplicated ($nbDuplicateNameGiven - $nbGeneDuplicated).";
  $stringPrint .= "$missing_gn_in_fasta_counter entries in db have no GN\n"; # JN: Tentative output

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
print_omniscient($hash_omniscient, $ostreamGFF); # JN: check content of $hash_omniscient

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

# Each mRNA of a gene has its proper gene name. Most often is the same, and annie added a number at the end.
# To provide only one gene name, we remove this number and then remove duplicate name (case insensitive).
# If it stay at the end of the process more than one name, they will be concatenated together.
# It removes redundancy intra name.
sub manageGeneNameBlast {

  my ($geneName) = @_;
  foreach my $element (keys %$geneName) {
    my @tab = @{$geneName->{$element}};

    my %seen;
    my @unique;
    for my $w (@tab) { # remove duplicate in list case insensitive
      $w =~ s/_[0-9]+$// ;
      next if $seen{lc($w)}++;
      push(@unique, $w);
    }

    my $finalName = "";
    my $cpt = 0;
    foreach my $name (@unique) { #if several names we will concatenate them together

        if ($cpt == 0) {
          $finalName .= "$name";
          $cpt++;
        }
        else {
          $finalName .= "_$name";
        }

    }
    $geneName->{$element} = $finalName;
    $nameBlast{lc($finalName)}++;
  }
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
    foreach my $element (@{$mRNAproduct{$refname}}) {
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
  my ($feature, $opt_output) = @_;

  my $functionAdded = undef;
  my $ID = lc($feature->_tag_value('ID'));
  foreach my $function_type (keys %functionData) {

    if (exists ($functionData{$function_type}{$ID})) {
      $functionAdded = "true";

      my $data_list;

      if (lc($function_type) eq "go") {
        foreach my $data (@{$functionData{$function_type}{$ID}}) {
          $feature->add_tag_value('Ontology_term', $data);
          $data_list .= "$data,";
          $functionDataAdded{$function_type}++;
        }
      }
      else {
        foreach my $data (@{$functionData{$function_type}{$ID}}) {
          $feature->add_tag_value('Dbxref', $data);
          $data_list .= "$data,";
          $functionDataAdded{$function_type}++;
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
  my($file_in, $opt_blastEvalue, $hash_mRNAGeneLink) = @_;

#################################################################################
####### Step 1 : CATCH all candidates (the better candidate for each mRNA)####### (with a gene name)

  my %candidates;

  while(my $line = <$file_in>) {
    my @values = split(/\t/, $line);
    my $l2_name = lc($values[0]);
    my $prot_name = $values[1];
    my @prot_name_sliced = split(/\|/, $values[1]);
    my $uniprot_id = $prot_name_sliced[1];
    print "uniprot_id: ".$uniprot_id."\n" if ($opt_verbose);
    my $evalue = $values[10];
    print "Evalue: ".$evalue."\n" if ($opt_verbose);

    #if does not exist fill it if over the minimum evalue
    if (! exists_keys(\%candidates, ($l2_name)) or @{$candidates{$l2_name}} > 3 ) { # the second one means we saved an error message as candidates we still have to try to find a proper one
      if ( $evalue <= $opt_blastEvalue ) {
        my $protID_correct = undef;

        if ( exists $allIDs{lc($prot_name)}) { # JN: Look for same entry as in fasta file
          $protID_correct = $allIDs{lc($prot_name)};
          my $header = $db->header( $protID_correct ); # JN: Get header from fasta db
          if (! $header =~ m/GN=/) {
            # JN: No gene name
            $ostreamLog->print( "No gene name (GN=) in this header $header\n") if ($opt_verbose or $opt_output);
            $candidates{$l2_name} = ["error", $evalue, $prot_name."-".$l2_name];
          }
          if ($header =~ /PE=([1-5])\s/) {
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
      my $protID_correct = undef;

      if ( exists $allIDs{lc($prot_name)}) {
        $protID_correct = $allIDs{lc($prot_name)};
        my $header = $db->header( $protID_correct );
        if (! $header =~ m/GN=/) {
          # JN: No gene name
          $ostreamLog->print("No gene name (GN=) in this header $header\n") if ($opt_verbose or $opt_output);
        }
        if ($header =~ /PE=([1-5])\s/) {
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

  my %geneName;
  my %linkBmRNAandGene;

  foreach my $l2 (keys %candidates) {
    if ( $candidates{$l2}[0] eq "error" ) {
      $ostreamLog->print("error nothing found for $candidates{$l2}[2]\n") if ($opt_verbose or $opt_output);
      next;
    }

    #Save uniprot id of the best match
    print "save for $l2  ".$candidates{$l2}[2]."\n" if ($opt_verbose);
    $mRNAUniprotIDFromBlast{$l2} = $candidates{$l2}[2];
    print "save for $l2  ".$candidates{$l2}[2]."\n" if ($opt_verbose);
    my $header = $candidates{$l2}[0];
    print "header: ".$header."\n" if ($opt_verbose);

    if ($header =~ m/(^[^\s]+)\s(.+?(?= \w{2}=))(.+)/) {
      my $protID = $1;
      my $description = $2;
      my $theRest = $3;
      $theRest =~ s/\n//g;
      $theRest =~ s/\r//g;
      my $nameGene = undef;
      push ( @{ $mRNAproduct{$l2} }, $description );

      #deal with the rest
      my %hash_rest;
      my $tuple = undef;
      while ($theRest) {
        ($theRest, $tuple) = stringCatcher($theRest);
        my ($type, $value) = split /=/, $tuple;
        #print "$protID: type:$type --- value:$value\n";
        $hash_rest{lc($type)} = $value;
      }

      if (exists($hash_rest{"gn"})) { # JN: Check for Gene name?
        $nameGene = $hash_rest{"gn"};

        if (exists_keys ($hash_mRNAGeneLink, ($l2)) ) { # JN: Gene name is only captured if key exists here 
          my $geneID = $hash_mRNAGeneLink->{$l2};
          #print "push $geneID $nameGene\n";
          push ( @{ $geneName{lc($geneID)} }, lc($nameGene) );
          push(@{ $linkBmRNAandGene{lc($geneID)}}, lc($l2)); # save mRNA name for each gene name
        }
        else {
          $ostreamLog->print( "No parent found for $l2 (defined in the blast file) in hash_mRNAGeneLink (created by the gff file).\n") if ($opt_verbose or $opt_output);
        }
      }
      else {
        # JN: No gene name
        $ostreamLog->print( "Header from the db fasta file doesn't match the regular expression: $header\n") if ($opt_verbose or $opt_output);
      }
    }
  }

  ####################################################
  ####### Step 3 : Manage NAME final gene name ####### several isoforms could have different gene name reported. So we have to keep that information in some way to report only one STRING to gene name attribute of the gene feature.
  ################# Remove redundancy to have only one name for each gene

  manageGeneNameBlast(\%geneName); # JN: check what manageGeneNameBlast will set and report


  ##########################################################
  ####### Step 4 : CLEAN NAMES REDUNDANCY inter gene #######

  my %geneNewNameUsed;
  foreach my $geneID (keys %geneName) {
    $nbGeneNameInBlast++;
    my @mRNAList = @{$linkBmRNAandGene{$geneID}};
    my $String = $geneName{$geneID};
    #print "$String\n";
    if (! exists( $geneNewNameUsed{$String})) {
      $geneNewNameUsed{$String}++;
      $geneNameBlast{$geneID} = $String;
      # link name to mRNA and and isoform name _1 _2 _3 if several mRNA
      my $cptmRNA = 1;
      if ($#mRNAList != 0) {
        foreach my $mRNA (@mRNAList) {
          $mRNANameBlast{$mRNA} = $String."_iso".$cptmRNA;
          $cptmRNA++;
        }
      }
      else {
        $mRNANameBlast{$mRNAList[0]} = $String;
      }
    }
    else { #in case where name was already used, we will modify it by adding a number like "_2"
      $nbDuplicateName++;
      $geneNewNameUsed{$String}++;
      my $nbFound = $geneNewNameUsed{$String};
      $String .= "_$nbFound";
      $geneNewNameUsed{$String}++;
      $geneNameBlast{$geneID} = $String;
      # link name to mRNA and and isoform name _1 _2 _3 if several mRNA
      my $cptmRNA = 1;
      if ($#mRNAList != 0) {
        foreach my $mRNA (@mRNAList) {
          $mRNANameBlast{$mRNA} = $String."_iso".$cptmRNA;
          $cptmRNA++;
        }
      }
      else {
        $mRNANameBlast{$mRNAList[0]} = $String;
      }
    }
  }
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
      if ( exists $hash_mRNAGeneLink->{$mRNAID}) { ## check if exists among our current gff annotation file analyzed
        $mRNAAssociatedToTerm{$db_name}{$mRNAID}++;
        $GeneAssociatedToTerm{$db_name}{$hash_mRNAGeneLink->{$mRNAID}}++;
      }
    }

    #check for interpro
    if ( $sizeList > 11 ) {
      my $db_name = "InterPro";
      my $interpro_value = $values[11];
      $interpro_value =~ s/\n//g;
      my $interpro_tuple = "InterPro:".$interpro_value;
      print "interpro dB: ".$interpro_tuple."\n" if ($opt_verbose);

      if (! grep( /^\Q$interpro_tuple\E$/, @{$functionData{$db_name}{$mRNAID}} )) { #to avoid duplicate
        $TotalTerm{$db_name}++;
        push ( @{$functionData{$db_name}{$mRNAID}}, $interpro_tuple );
        if ( exists $hash_mRNAGeneLink->{$mRNAID}) { ## check if exists among our current gff annotation file analyzed
          $mRNAAssociatedToTerm{$db_name}{$mRNAID}++;
          $GeneAssociatedToTerm{$db_name}{$hash_mRNAGeneLink->{$mRNAID}}++;
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

        if (! grep( /^\Q$go_tuple\E$/, @{$functionData{$db_name}{$mRNAID}} )) { #to avoid duplicate
          $TotalTerm{$db_name}++;
          push ( @{$functionData{$db_name}{$mRNAID}}, $go_tuple );
          if ( exists $hash_mRNAGeneLink->{$mRNAID}) { ## check if exists among our current gff annotation file analyzed
            $mRNAAssociatedToTerm{$db_name}{$mRNAID}++;
            $GeneAssociatedToTerm{$db_name}{$hash_mRNAGeneLink->{$mRNAID}}++;
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

        if (! grep( /^\Q$pathway_tuple\E$/, @{$functionData{$db_name}{$mRNAID}} ) ) { # to avoid duplicate
          $TotalTerm{$db_name}++;
          push ( @{$functionData{$db_name}{$mRNAID}} , $pathway_tuple );
          if ( exists $hash_mRNAGeneLink->{$mRNAID}) { ## check if exists among our current gff annotation file analyzed
            $mRNAAssociatedToTerm{$db_name}{$mRNAID}++;
            $GeneAssociatedToTerm{$db_name}{$hash_mRNAGeneLink->{$mRNAID}}++;
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

    agat_sp_manage_functional_annotation.pl -f infile.gff [-b blast_infile][-d uniprot.fasta][-i interpro_infile.tsv][-id ABCDEF][-o output]
    agat_sp_manage_functional_annotation.pl --help

=head1 OPTIONS

=over 8

=item B<-f>, B<--reffile>,B<-ref> , B<--gff> or B<--gff3>

String - Input GTF/GFF file.

=item B<-b> or B<--blast>

String - Input blast ( outfmt 6 = tabular ) file that will be used to complement the features
read from the first file (specified with --ref).

=item B<-d> or B<--db>

String - The fasta file that has been used as DB for the blast. Gene names and products/descriptions will be fished from this file.

=item B<--be> or B<--blast_evalue>

Integer - Maximum e-value to keep the annotation from the blast file. By default 1e-6.

=item B<--pe>

Integer - The PE (protein existence) in the uniprot header indicates the type of evidence that supports the existence of the protein.
You can decide until which protein existence level you want to consider to lift the finctional information. Default 5.

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

=item B<-o> or B<--output>

String - Output folder name with summary files. If no output file is specified, the output will be
written to STDOUT.

=item B<-v>

Boolean - Verbose, for debug purpose.

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
