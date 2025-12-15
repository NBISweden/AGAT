#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Try::Tiny;
use File::Basename;
use Bio::DB::Taxonomy;
use Getopt::Long;
use Bio::DB::Fasta;
use Pod::Usage;
use List::MoreUtils qw(uniq);
use AGAT::AGAT;


start_script();
my $header = get_agat_header();
# -------------------------------- LOAD OPTIONS --------------------------------
#The cases are exclusive, one result could not be part of several cases.
my %cases_explanation = (
  -1  => "No protein alignement overlap the gene model.",
  0   => "There is protein alignment that overlap the gene model but none goes over the threshold defined.",
  1   => "There is protein alignment that overlap the gene model, the overlap is over the threshold defined, the protein is from one of the species defined, and the PE value is as defined. P.S: Case only possible when both options sp and pe are activated.",
  2.1  => "There is protein alignment that overlap the gene model, the overlap is over the threshold defined, the PE value is as defined.",
  2.2  => "There is protein alignment that overlap the gene model, the overlap is over the threshold defined, the protein is from one of the species defined.",
  3   => "There is protein alignment that overlap the gene model, the overlap is over the threshold defined."
);

my $opt_output             = undef;
my $annotation_gff         = undef;
my $protein_gff            = undef;
my $protein_fasta          = undef;
my $valueK                 = 50;
my $opt_test               = "=";
my $attributes             = undef ;
my $whole_sequence_opt     = undef ;
my $priority_opt           = "pe";
my $sort_method_by_species = undef ;
my $sort_method_by_pe      = undef ;
my $method_opt             = "replace";
my $opt_help               = 0;

my @copyARGV=@ARGV;

# Split shared vs script options and parse script options
my ($shared_argv, $script_argv) = split_argv_shared_vs_script(\@ARGV);
# Parse script-specific options
my $parser = Getopt::Long::Parser->new();
$parser->configure('bundling','no_auto_abbrev');
if ( !$parser->getoptionsfromarray(
  $script_argv,
  "h|help"                 => \$opt_help,
  "annotation|a=s"         => \$annotation_gff,
  "pgff=s"                 => \$protein_gff,
  "sp:s"                   => \$sort_method_by_species,
  'test=s'                 => \$opt_test,
  'pe:i'                   => \$sort_method_by_pe,
  'priority|p=s'           => \$priority_opt,
  "fasta|pfasta=s"         => \$protein_fasta,
  "w"                      => \$whole_sequence_opt,
  "value|threshold=i"      => \$valueK,
  'method|m:s'             => \$method_opt,
  "output|out|o=s"         => \$opt_output ) )
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

if ( ! ($annotation_gff and $protein_gff and $protein_fasta) ){
    pod2usage( {
           -message => "$header\nAt least 3 parameters are mandatory:\nAnnotation gff file (-a), Protein gff file (--pgff), Protein fasta file (--pfasta)\n",
           -verbose => 0,
           -exitval => 1 } );
}

# --- Manage config ---
my $shared_opts = parse_shared_options($shared_argv);
initialize_agat({ config_file_in => ( $shared_opts->{config} ), input => $annotation_gff, shared_opts => $shared_opts });

#               +------------------------------------------------------+
#               |+----------------------------------------------------+|
#               ||                     Manage OPTIONS                 ||
#               |+----------------------------------------------------+|
#               +------------------------------------------------------+



#Manage test option
if($opt_test ne "<" and $opt_test ne ">" and $opt_test ne "<=" and $opt_test ne ">=" and $opt_test ne "="){
  die "The test to apply is Wrong: $opt_test.\nWe want something among this list: <,>,<=,>= or =.";
}
if(defined ($sort_method_by_pe)){
  if ($sort_method_by_pe){
    if(! (5 >= $sort_method_by_pe and  $sort_method_by_pe >= 1) ){
      die "The value of the Protein Existence value is Wrong: $sort_method_by_pe.\n We want a value between 1 and 5.";
    }
  }
  else{
    $sort_method_by_pe=1;
  }
}

##########################
##### Manage Output ######

our @outputTab;

if (defined($opt_output) ) {
  if (-f $opt_output){
    die "Cannot create a directory with the name $opt_output because a file with this name already exists.\n";
  }
  if (-d $opt_output){
    die "The output directory choosen already exists. Please geve me another Name.\n";
  }

  #create the folder
  mkdir $opt_output;

  my $outfile_gff = $annotation_gff;
  my ($file1,$dir1,$ext1) = fileparse($outfile_gff, qr/\.[^.]*/);
  $outfile_gff = $file1."_updated.gff";

  #0 txt
  my $ostreamReport = IO::File->new(">".$opt_output."/report.txt" ) or
  croak( sprintf( "Can not open '%s' for writing %s", $opt_output."/report.txt", $! ));
  push (@outputTab, $ostreamReport);
  #1 txt
  my $ostreamFAadded = IO::File->new(">".$opt_output."/function_added.txt" ) or
  croak( sprintf( "Can not open '%s' for writing %s", $opt_output."/function_added.txt", $! ));
  push (@outputTab, $ostreamFAadded);
  _print1("Gene ID\tmRNA ID\tGene name\tProduct\n", 1);
  #2 gff
  my $ostreamCoding = Bio::Tools::GFF->new(-file => ">".$opt_output."/".$outfile_gff, -gff_version => $CONFIG->{gff_output_version} ) or
  croak(sprintf( "Can not open '%s' for writing %s", $opt_output."/".$outfile_gff, $! ));
  push (@outputTab, $ostreamCoding);

}
else{
  my $ostreamReport = \*STDOUT ;
  push (@outputTab, $ostreamReport);
  my $ostreamFAadded = \*STDOUT ;
  push (@outputTab, $ostreamFAadded);
   my $ostream  = IO::File->new();
  $ostream->fdopen( fileno(STDOUT), 'w' ) or croak( sprintf( "Can not open STDOUT for writing: %s", $! ) );
  my $outputGFF = Bio::Tools::GFF->new( -fh => $ostream, -gff_version => $CONFIG->{gff_output_version} ) or croak( sprintf( "Can not open STDOUT for writing: %s", $! ) );
  push (@outputTab, $outputGFF);
}

# Manage species names
my $db = Bio::DB::Taxonomy->new(-source => 'entrez');

if(defined($sort_method_by_species) ){
  if($sort_method_by_species){
    my %hash;
    my @listValue= split /,/, $sort_method_by_species;
    my $counter=1;
    foreach my $element (@listValue ){
      $element =~ s/^\s+|\s+$//g; #removing leading and trailing white spaces
      my @taxonids = $db->get_taxonids($element);
      if(@taxonids){
        $hash{$counter++}=$element;
      }
      else{
        _print1("/!\\ species <$element> is unknown in the NCBI taxonomy DB, we skip it !\n",0);
      }
    }
    $sort_method_by_species = \%hash;
  }
  else{
    $sort_method_by_species = _taxid_ref_sorted();
  }
  _print1 ("Priority in this order will be used for selecting the referential protein form matching proteins:\n");
  foreach my $priority (sort {$a <=> $b} keys %{$sort_method_by_species}){
    _print1( $priority." - ".$sort_method_by_species->{$priority}."\n",0);
  }
}

#Manage priority
if ($priority_opt ne "pe" and $priority_opt ne "sp"){
  die "Priority option with value $priority_opt doesn't exist. Please select pe or sp. (i.e help)\n";
}

#Manage method
if ($method_opt eq "replace"){
  _print1 ("We will add or replace the product and Name values when a protein maps properly.\n");
}
elsif($method_opt eq "add"){
  _print1 ("We will add the lfp_product and lfp_name values when a protein maps properly.\n");
}
elsif($method_opt eq "complete"){
  _print1 ("We will add the product and Name values when a protein maps properly and no product and/or Name value exists.\n");
}
else{ die "Method option must be replace, add or complete. Please check the help for more information. (replace by default)\n";}


#               +------------------------------------------------------+
#               |+----------------------------------------------------+|
#               ||                          MAIN                      ||
#               |+----------------------------------------------------+|
#               +------------------------------------------------------+

######################
### Parse GFF input #
_print1( "Parsing file $annotation_gff\n",0);
my ($hash_omniscient) = slurp_gff3_file_JD({ input => $annotation_gff });
_print1( "Done\nParsing file $protein_gff\n",0);
my ($prot_omniscient) = slurp_gff3_file_JD({ input => $protein_gff });
_print1( "Done\n",0);

###########################
#### Parse protein fasta #
_print1( "Parsing the protein fasta file $protein_fasta\n",0);
my $nbFastaSeq=0;
my $db_prot = Bio::DB::Fasta->new($protein_fasta);
my @long_ids_prot      = $db_prot->get_all_primary_ids;
my %allIDs_prot; # save ID in lower case to avoid cast problems
foreach my $long_id (@long_ids_prot){
  #uniprot case long_id=>sp|Q5R8S7|PPIA_PONPY short_id=Q5R8S7
  my $short_id = _take_clean_id($long_id);
  $allIDs_prot{lc($short_id)}=$long_id;
}
_print1( "Done\n",0);

my $omniscient1_sorted = gather_and_sort_l1_location_by_seq_id_and_strand($hash_omniscient);
my $omniscient2_sorted = gather_and_sort_l1_location_by_seq_id_and_strand($prot_omniscient);
my $topfeatures = get_feature_type_by_agat_value($hash_omniscient, 'level1', 'topfeature');

my %cases;

foreach my $locusID ( keys %{$omniscient1_sorted}){ # tag_l1 = protein_match match ....
  foreach my $tag_l1 ( keys %{$omniscient1_sorted->{$locusID}} ) {

    #skip top features
    if(exists_keys($topfeatures,($tag_l1))){ next; }

    # Go through location from left to right ### !!
    my @aligns;
    my $selected=undef;
    while ( @{$omniscient1_sorted->{$locusID}{$tag_l1}} ){

      my $location = shift  @{$omniscient1_sorted->{$locusID}{$tag_l1}};# This location will be updated on the fly
      my $id1_l1 = lc($location->[0]);

      if( exists_keys($omniscient2_sorted, ($locusID ) ) ) {

        foreach my $tag2_l1 ( keys %{$omniscient2_sorted->{$locusID}} ) {

          #skip top features
          if(exists_keys($topfeatures,($tag2_l1))){ next; }

          while ( @{$omniscient2_sorted->{$locusID}{$tag2_l1}} ){

            my $location2 = shift  @{$omniscient2_sorted->{$locusID}{$tag2_l1}};# This location will be updated on the fly
            my $id2_l1 = lc($location2->[0]);

            #If location_to_check start if over the end of the reference location, we stop
            if($location2->[1] > $location->[2]) {last;}
            #If location_to_check end if inferior to the start of the reference location, we continue next
            if($location2->[2] < $location->[1]) {next;}

            # Let's check at Gene LEVEL
            if( location_overlap($location, $location2) ){ #location overlap at gene level check now level3

              ###################################
              # let's check at the deeper level #
              my $prot_tag = $prot_omniscient->{'level1'}{$tag2_l1}{$id2_l1}->_tag_value('Name');
              my $align = check_gene_overlap_gffAlign($hash_omniscient, $prot_omniscient, $id1_l1, $id2_l1, $prot_tag); #If contains CDS it has to overlap at CDS level to be merged, otherwise any type of feature level3 overlaping is sufficient to decide to merge the level1 together
              #print Dumper( $align);
              if(@{$align} ){ #check is not empty
                # @list_res = ($gene_id2, $w_overlap12, $w_overlap12_abs, $w_overlap21, $w_overlap21_abs, $overlap12, $overlap12_abs, $overlap21, $overlap21_abs);
                # align contains: level1 id of the protein into the gff file (gene_id2),
                ## overlap percent whole gene model against protein  (w_overlap12),
                # absolute overlap percent whole gene model against the protein  ($w_overlap12_abs). Absolute means we check the cigar annotation of the protein to check the shift in the reading frame and take in account the few nucleotide in more or in less
                ## overlap percent of the protein  against the whole gene model ($w_overlap21),
                # absolute overlap percent of the protein  against the whole gene model ($w_overlap21_abs),
                ## overlap percent cds part of the gene model against protein  (overlap12),
                # absolute overlap percent cds part of the gene model against the protein  ($overlap12_abs). Absolute means we check the cigar annotation of the protein to check the shift in the reading frame and take in account the few nucleotide in more or in less
                ## overlap percent of the protein against the cds part of the gene model ($overlap21),
                # absolute overlap percent of the protein against the cds part of the gene model ($overlap21_abs),
                # absolute overlap percent of the protein  against the whole gene model but rationalized by the total length (prot+gene-overlap) $w_overlap_JD_abs
                # absolute overlap percent of the protein  against the coding part of the gene model but rationalized by the total length (prot+gene-overlap) $overlap_JD_abs
                push @aligns, [@{$align}];
              }
            }
          }
        }
      }

#               +------------------------------------------------------+
#               |+----------------------------------------------------+|
#               ||       NOW WE HAVE ALL THE ALIGNEMENT VALUE         ||
#               |+----------------------------------------------------+|
#               +------------------------------------------------------+
                                # Let's filter them
      #We have at least one prot aligned to this gene model
      my @aligns_filtered;
      if(@aligns){

        # Manage option coding sequence or whole sequence
        my $col = 6;
        if($whole_sequence_opt){$col = 5;}

        #Sort results and filter by the overlap value threshold
        foreach my $result ( sort {$b->[$col] <=> $a->[$col] }  @aligns ){

          # filter by value Threshold
          if($result->[$col] > $valueK){
            push @aligns_filtered, $result;
          }
        }

        #print @aligns_filtered results
        _print2( "\n\nprotein aligned to the gene $id1_l1 over the threshold $valueK:\n", 0);
        foreach my $result ( @aligns_filtered){
          if($result->[$col] > $valueK){
            _print2( "@$result\n", 0 );
          }
        }
        _print2("\n", 0);

        if(@aligns_filtered){

          #########################################
          # 1) filter by pe and specific species
          if($sort_method_by_pe and $sort_method_by_species){
            _print2( "get_result_sort_method_by_pe_and_species case 1 !\n",0);
            $selected = get_result_sort_method_by_pe_and_species(\@aligns_filtered, $col, $sort_method_by_pe, $opt_test, $sort_method_by_species);
            if($selected){$cases{1}++;}
          }

          #####################################################
          # 2.1) filter by pe or if need be by specific species
          if(!$selected  and $priority_opt eq "pe"){

            # filter by protein existence value
            if(! $selected and $sort_method_by_pe){
              _print2( "pe case 2.1.1!\n", 0);
              $selected = get_result_sort_method_by_pe(\@aligns_filtered, $col, $sort_method_by_pe, $opt_test);
              if($selected){$cases{211}++;$cases{21}++;}
            }

            # filter by specific species
            if(! $selected and $sort_method_by_species){
              _print2( "sort_method_by_species case 2.1.2!\n", 0);
              $selected = get_result_sort_method_by_species($sort_method_by_species, \@aligns_filtered, $col);
              if($selected){$cases{212}++;$cases{22}++;}
            }
          }
          #####################################################
          # 2.2) filter by specific species or if need be by pe
          elsif(! $selected  and $priority_opt eq "sp"){
            #########################################
            # filter by specific species
            if(! $selected and $sort_method_by_species){
              _print2( "sort_method_by_species case 2.2.1!\n", 0);
              $selected = get_result_sort_method_by_species($sort_method_by_species, \@aligns_filtered, $col);
               if($selected){$cases{221}++;$cases{22}++;}
            }

            #########################################
            # filter by protein existence value
            if(! $selected and $sort_method_by_pe){
              _print2( "pe case 2.2.2!\n", 0);
              $selected = get_result_sort_method_by_pe(\@aligns_filtered, $col, $sort_method_by_pe, $opt_test);
              if($selected){$cases{222}++;$cases{21}++;}
            }
          }

          #########################################
          # 3) Take the first = the best overlap value
          if(! $selected){
            _print2( "Normal case 3!\n", 0 );
            # read from best value to the lowest one
            $selected = $aligns_filtered[0];
            $cases{3}++;
          }

          _print2( "Protein selected =  $selected\n",0);


#               +------------------------------------------------------+
#               |+----------------------------------------------------+|
#               ||       NOW WE ATTACH THE INFORMATION FOUND          ||
#               |+----------------------------------------------------+|
#               +------------------------------------------------------+

          #Modify l1
          my $name_added=1;
          my $geneName = _get_gn($selected);
          my $feature_l1 = $hash_omniscient->{'level1'}{$tag_l1}{$id1_l1};
          if( $method_opt ne "add"){
            if( ! ($method_opt eq "complete" and $feature_l1->has_tag('Name') ) ){
              create_or_replace_tag($feature_l1, 'Name', $geneName);
            }else{$name_added=0;}
          }
          else{
            create_or_replace_tag($feature_l1, 'lfp_name', $geneName);
          }

          #Now Modify l2
          my $product_added=1;
          my $product = _get_function($selected);
          foreach my $tag_l2 (keys %{$hash_omniscient->{'level2'}}){
            if( exists_keys($hash_omniscient,('level2', $tag_l2, $id1_l1))){
              foreach my $feature_l2 (@{$hash_omniscient->{'level2'}{$tag_l2}{$id1_l1}}){
                #Modify each l2
                if( $method_opt ne "add"){
                  if( ! ($method_opt eq "complete" and $feature_l1->has_tag('product') ) ){
                    create_or_replace_tag($feature_l2, 'product', $product);
                  }else{$product_added=0;}
                }
                else{
                  create_or_replace_tag($feature_l2, 'lfp_product', $product);
                }

                if($name_added or $product_added){
                  if($name_added and $product_added){
                    _print($id1_l1."\t".$feature_l2->_tag_value('ID')."\t".$geneName."\t".$product."\n", 1);
                  }
                  elsif($name_added){
                    _print($id1_l1."\t".$feature_l2->_tag_value('ID')."\t".$geneName."\t-\n", 1);
                  }
                  else{
                    _print($id1_l1."\t".$feature_l2->_tag_value('ID')."\t-\t".$product."\n", 1);
                  }
                }
              }
            }
          }



        }
        else{
          $cases{0}++;
          _print2( "No protein overlap over the threshold $valueK for this gene model: $id1_l1\n", 0);
        }
      }
      else{
        $cases{-1}++;
        _print2( "No protein aligned to this gene model: $id1_l1\n", 0);
      }
    }
  }

}

######################
##### Print the result
print_omniscient( {omniscient => $hash_omniscient, output => $outputTab[2]} );

##########################################
#### SUMMERIZING##########################
if(keys %cases){_print( "\nThe liftover of function from proteins to the gene models has been done as following:\n",0);}
foreach my $key (keys %cases){
  if( $key == 211 or $key == 212 or $key == 221 or $key == 222 ) {
     _print1( "we have $cases{$key} cases $key\n",0);
  }
  else{
    _print1( "we have $cases{$key} cases $key\n",0);
  }
}
_print1( "\nLet's remind the different cases:\n",0);
foreach my $key (sort{$a <=> $b} keys %cases_explanation){
  _print1( "$key => $cases_explanation{$key}\n",0);
}

# Graceful end of script before subroutines
end_script();

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

############################
## PROTEIN PRIORETIZATION

#1) keep score as calculated and select the best

#2) by ProteinExistence value (take the best score in better or equal value of the PE expected)
#   1. Experimental evidence at protein level
#   2. Experimental evidence at transcript level
#   3. Protein inferred from homology
#   4. Protein predicted
#   5. Protein uncertain

#3) Taking into account a specific species following an order of choices
# It's an arbitrary choice
sub _taxid_ref_sorted{

  my %sorted = ( 1 => 'Homo sapiens',
    2 => 'Mus musculus',
    3 => 'Drosophila melanogaster',
    4 => 'Caenorhabditis elegans',
    5 => 'Arabidopsis thaliana'
  );

  return \%sorted;
}

sub get_result_sort_method_by_pe_and_species{
  my ($aligns, $col, $sort_method_by_pe, $opt_test, $sort_method_by_species)=@_;

  my $selected = undef;

  my $counter=1;
  while(! $selected){

    if(! exists ($sort_method_by_species->{$counter}) ){last;}

    # read from best value to the lowest one
    foreach my $result ( @{$aligns} ){

      my $species = _get_species($result);

      if(lc($species) eq lc($sort_method_by_species->{$counter}) ){

        my $pe = _get_pe($result);

        if ($opt_test eq ">"){
          if ($pe > $sort_method_by_pe){
            $selected=$result;last;
          }
        }
        if ($opt_test eq "<"){
          if ($pe < $sort_method_by_pe){
            $selected=$result;last;
          }
        }
        if ($opt_test eq "="){
          if ($pe == $sort_method_by_pe){
            $selected=$result;last;
          }
        }
        if ($opt_test eq "<="){
          if ($pe <= $sort_method_by_pe){
            $selected=$result;last;
          }
        }
        if ($opt_test eq ">="){
          if ($pe >= $sort_method_by_pe){
            $selected=$result;last;
          }
        }
      }
    }
    $counter++;
  }

  return $selected;

}

sub get_result_sort_method_by_species{
  my ($sort_method_by_species, $aligns, $col) = @_;

  my $selected = undef;

  my $found_sp=undef;
  my $counter=1;
  while(! $selected){

    if(! exists ($sort_method_by_species->{$counter}) ){last;}

    # read from best value to the lowest one
    foreach my $result ( @{$aligns} ){

      my $species = _get_species($result);
      #print "<$species> vs >$sort_method_by_species->{$counter}<\n";

      if(lc($species) eq lc($sort_method_by_species->{$counter}) ){
        $selected=$result;last;
      }
    }
    $counter++;
  }

  return $selected;
}

sub get_result_sort_method_by_pe{
  my ($aligns, $col, $sort_method_by_pe, $opt_tes ) = @_;

  my $selected = undef;

  # read from best value to the lowest one
  foreach my $result ( @{$aligns} ){

    my $pe = _get_pe($result);

    if ($opt_test eq ">"){
      if ($pe > $sort_method_by_pe){
        $selected=$result;last;
      }
    }
    if ($opt_test eq "<"){
      if ($pe < $sort_method_by_pe){
        $selected=$result;last;
      }
    }
    if ($opt_test eq "="){
      if ($pe == $sort_method_by_pe){
        $selected=$result;last;
      }
    }
    if ($opt_test eq "<="){
      if ($pe <= $sort_method_by_pe){
        $selected=$result;last;
      }
    }
    if ($opt_test eq ">="){
      if ($pe >= $sort_method_by_pe){
        $selected=$result;last;
      }
    }
  }

  return $selected;
}

#get function
sub _get_function{
  my ($self)=@_;
  my $header = $self->[$#{$self}];
  my $egal = index($header, '=');
  my $function  = substr $header, 0, $egal-2;

  $function =~ s/^\s+|\s+$//g; #removing leading and trailing white spaces

  return $function ;
}

#get species
sub _get_species{
  my ($self)=@_;

  my $species = undef;

  my $header = $self->[$#{$self}];

  my $egal = index($header, '=');
  my $abb = substr $header, $egal-2, 2;
  my $clipped = substr $header, $egal+1;

  while (lc($abb) ne "os") {
    $egal = index($clipped, '=');
    if($egal == -1){last;}
    $abb = substr $clipped, $egal-2, 2;
    $clipped = substr $clipped, $egal+1;
  }
  if($egal == -1){ warn("No species name found in this fasta header: $self");return $species;}
  $egal = index($clipped, '=');
  if($egal != -1){
    $species  =  substr $clipped, 0, $egal-2;
  }
  else{
    $species  =  $clipped;
  }
  $species =~ s/^\s+|\s+$//g; #removing leading and trailing white spaces

  return $species;
}

#get gene name
sub _get_gn{
  my ($self)=@_;

  my $geneName = undef;

  my $header = $self->[$#{$self}];

  my $egal = index($header, '=');
  my $abb = substr $header, $egal-2, 2;
  my $clipped = substr $header, $egal+1;

  while (lc($abb) ne "gn") {
    $egal = index($clipped, '=');
    if($egal == -1){last;} #no = character found

    $abb = substr $clipped, $egal-2, 2;
    $clipped = substr $clipped, $egal+1;
  }
  if($egal == -1){ warn("No gene name found in this fasta header: $self");return $geneName;}
  $egal = index($clipped, '=');
  if($egal != -1){
    $geneName  =  substr $clipped, 0, $egal-2;
  }
  else{
    $geneName  =  $clipped;
  }

  $geneName =~ s/^\s+|\s+$//g; #removing leading and trailing white spaces

  return $geneName ;
}

#get ProteinExistence
sub _get_pe{
  my ($self)=@_;

  my $pe = undef;

  my $header = $self->[$#{$self}];

  my $egal = index($header, '=');
  my $abb = substr $header, $egal-2, 2;
  my $clipped = substr $header, $egal+1;

  while (lc($abb) ne "pe") {
    $egal = index($clipped, '=');
    if($egal == -1){last;} #no = character found

    $abb = substr $clipped, $egal-2, 2;
    $clipped = substr $clipped, $egal+1;
  }
  if($egal == -1){ warn("No pe found in this fasta header: $self"); return $pe; }
  $egal = index($clipped, '=');
  if($egal != -1){
    $pe  =  substr $clipped, 0, $egal-2;
  }
  else{
    $pe  =  $clipped;
  }

  $pe =~ s/^\s+|\s+$//g; #removing leading and trailing white spaces

  return $pe ;
}

#get sequence version
sub _get_sv{
  my ($self)=@_;

  my $sv = undef;

  my $header = $self->[$#{$self}];

  my $egal = index($header, '=');
  my $abb = substr $header, $egal-2, 2;
  my $clipped = substr $header, $egal+1;

  while (lc($abb) ne "sv") {
    $egal = index($clipped, '=');
    if($egal == -1){last;} #no = character found

    $abb = substr $clipped, $egal-2, 2;
    $clipped = substr $clipped, $egal+1;
  }
  if($egal == -1){ warn("No sv found in this fasta header: $self"); return $sv; }
  $egal = index($clipped, '=');
  if($egal != -1){
    $sv  =  substr $clipped, 0, $egal-2;
  }
  else{
    $sv  =  $clipped;
  }

  $sv =~ s/^\s+|\s+$//g; #removing leading and trailing white spaces

  return $sv ;
}

sub  _get_sequence{
  my  ($db, $seq_id) = @_;

  my $sequence = "";
  my $descritpion = "";
  my $seq_id_correct = _take_clean_id($seq_id);
  if( exists $allIDs_prot{lc($seq_id_correct)}){

    my $seq_id_original= $allIDs_prot{lc($seq_id_correct)};

    $sequence = $db->subseq($seq_id_original);
    $descritpion = (split(/\s+/, $db->header($seq_id_original), 2))[1]; #take header and remove the first element wihch is the seq_id_original

    if($sequence eq ""){
      warn "Problem ! no sequence extracted for - $seq_id_correct !\n";  die;
    }
  }
  else{
    warn "Problem ! protein ID $seq_id_correct not found into the protein fasta file!\n";
  }

  return length($sequence), $seq_id_correct, $descritpion;
}

sub _take_clean_id {
  my ($id) = @_ ;

  my $correct_id=$id;

  if($correct_id =~ m/\|/){
    my @tmp = split /\|/, $correct_id;
    $correct_id = $tmp[1];
  }
  if($correct_id =~ m/\./){
    my @tmp = split /\./, $correct_id;
    $correct_id = $tmp[0];
  }
  if($correct_id =~ m/\-/){
    my @tmp = split /\-/, $correct_id;
    $correct_id = $tmp[0];
  }
  return $correct_id;
}

#Check if two kind of L2 overlap at l3
sub check_gene_overlap_gffAlign{
  my  ($hash_omniscient, $prot_omniscient, $gene_id, $gene_id2, $prot_tag)=@_;

#  my $overlap12=undef;
  my $overlap12_abs=undef;

#  my $w_overlap12=undef;
  my $w_overlap12_abs=undef;

#  my $overlap21=undef;
  my $overlap21_abs=undef;

#  my $w_overlap21=undef;
  my $w_overlap21_abs=undef;

  my $w_overlap_JD_abs = -1;
  my $overlap_JD_abs = -1;
  my @list_res;


  #print "IN $gene_id, $gene_id2\n";
  foreach my $l2_type (keys %{$hash_omniscient->{'level2'}} ){

    if(exists_keys($hash_omniscient,('level2', $l2_type, lc($gene_id)))){

      foreach my $mrna_feature ( sort {$a->start <=> $b->start } @{$hash_omniscient->{'level2'}{$l2_type}{lc($gene_id)}}){ # from here bothe feature level2 are the same type
        my $mrna_id = $mrna_feature->_tag_value('ID');


        #hash_omniscient contains match
        foreach my $l1_type (keys %{$prot_omniscient->{'level1'}} ){

          #check full CDS for each mRNA
          if(exists_keys($prot_omniscient,('level1', $l1_type, lc($gene_id2)))){

            #calcul lenght2
            my $lenght2=0;
            foreach my $tag_l2 (keys %{$prot_omniscient->{'level2'}}){

              if(exists_keys($prot_omniscient,('level2', $tag_l2, lc($gene_id2)))){
                foreach my $feature2 (@{$prot_omniscient->{'level2'}{$tag_l2}{lc($gene_id2)}}){
                  #print$feature2->end." -  ".$feature2->start."\n";
                  #print $feature2->end - $feature2->start." + 1\n";
                  $lenght2 += ($feature2->end - $feature2->start + 1);
                }
              }
            }
            #print "lenght protein: $lenght2\n";


            my $w_overlap=0;
            my $w_abs_overlap=0;
            my $w_lenght1=0;
            ####################################
            # CALCUL ONTO THE WHOLE GENE MODEL #
            my @list_tag_l3=('exon');
            if(! exists_keys( $hash_omniscient, ('level3','exon'))){
              warn "No exon found into the annoation file for feature $gene_id, we will use all the other l3 features\n";
              foreach my $tag_l3 (keys %{$hash_omniscient->{'level3'}}){
                push @list_tag_l3,$tag_l3;
              }
            }

            foreach my $tag_l3 ( @list_tag_l3 ){

              if(exists_keys($hash_omniscient,('level3', $tag_l3, lc($mrna_id)))){
                foreach my $feature1 ( sort {$a->start <=> $b->start}  @{$hash_omniscient->{'level3'}{$tag_l3}{lc($mrna_id)}}){

                  #print "annot location: ".$feature1->start." ".$feature1->end."\n";
                  $w_lenght1 += ($feature1->end - $feature1->start + 1);
                  #print $feature1->gff_string."\n";

                  foreach my $tag_l2 (keys %{$prot_omniscient->{'level2'}}){

                    if(exists_keys($prot_omniscient,('level2', $tag_l2, lc($gene_id2)))){
                      foreach my $feature2 (sort {$a->start <=> $b->start }  @{$prot_omniscient->{'level2'}{$tag_l2}{lc($gene_id2)}}){
                        #print "    ".$feature2->gff_string."\n";
                        if(($feature2->start <= $feature1->end) and ($feature2->end >= $feature1->start )){ # they overlap

                          if($feature2->start > $feature1->end) {last;}
                          if($feature2->end < $feature1->start) {next;}
                          #print "prot location: ".$feature2->start." ".$feature2->end."\n";
                          my $start = $feature2->start > $feature1->start ? $feature2->start : $feature1->start;

                          my $end = $feature2->end < $feature1->end ? $feature2->end :  $feature1->end;

                             my $chunck_abs_overlap += get_absolute_match($feature2, $start, $end);
                             $w_abs_overlap += $chunck_abs_overlap;
                             #print "chunck_abs_overlap = $chunck_abs_overlap\n";
                             #print "chunck_overlap= ".($end - $start + 1)."\n";
                             $w_overlap += ($end - $start + 1);

                        }
                      }
                    }
                  }
                }
                #print "gene_id $gene_id, gene_id2 $gene_id2 <=> w_overlap = $w_overlap - w_abs_overlap = $w_abs_overlap\n";
              }
            }

            #1 -> 2
            #$w_overlap12 = sprintf "%.1f", ($w_overlap*100/$w_lenght1);
            $w_overlap12_abs = sprintf "%.1f", ($w_abs_overlap*100/$w_lenght1);


            my $overlap=0;
            my $abs_overlap=0;
            my $lenght1=0;
            # ##########################################################################
            # # CALCUL ONTO THE CODING SEQUENCE PART OF THE GENE MODEL ONLY (SKIP UTR) #
            if(exists_keys( $hash_omniscient, ('level3','cds',lc($mrna_id) ))){

              foreach my $feature1 ( sort {$a->start <=> $b->start}  @{$hash_omniscient->{'level3'}{'cds'}{lc($mrna_id)}}){

            #     #print "annot location: ".$feature1->start." ".$feature1->end."\n";
                $lenght1 += ($feature1->end - $feature1->start + 1);
            #     print $feature1->gff_string."\n";

                foreach my $tag_l2 (keys %{$prot_omniscient->{'level2'}}){

                  if(exists_keys($prot_omniscient,('level2', $tag_l2, lc($gene_id2)))){
                    foreach my $feature2 (sort {$a->start <=> $b->start }  @{$prot_omniscient->{'level2'}{$tag_l2}{lc($gene_id2)}}){

                      if(($feature2->start <= $feature1->end) and ($feature2->end >= $feature1->start )){ # they overlap

                        if($feature2->start > $feature1->end) {last;}
                        if($feature2->end < $feature1->start) {next;}

            #             print "We overlap !\n";
                        #print "prot location: ".$feature2->start." ".$feature2->end."\n";
                        my $start = $feature2->start > $feature1->start ? $feature2->start : $feature1->start;

                        my $end = $feature2->end < $feature1->end ? $feature2->end :  $feature1->end;

                           my $chunck_abs_overlap += get_absolute_match($feature2, $start, $end);
                           $abs_overlap += $chunck_abs_overlap;
                           #print "chunck_abs_overlap = $chunck_abs_overlap\n";
                           #print "chunck_overlap= ".($end - $start + 1)."\n";
                           $overlap += ($end - $start + 1);

                      }
                      #print "overlap = $overlap - abs_overlap = $abs_overlap\n";
                    }
                  }
                }
              }
            #1 -> 2
            #$overlap12 = sprintf "%.1f", ($overlap*100/$lenght1);
            $overlap12_abs = sprintf "%.1f", ($abs_overlap*100/$lenght1);
            }

            ###########################
            # CORRECT THE 21 value by the real length of the protein (Not the whole original protein is aligned into the protein gff file)
            try{
              my ($protLenSeqOriginal, $proteinName, $descritpion) = _get_sequence($db_prot, $prot_tag);
              my $lenght_corrected = $protLenSeqOriginal*3;
              #2 -> 1 whole sequence
              #$w_overlap21 = sprintf "%.1f", ($w_overlap*100/$lenght_corrected);
              $w_overlap21_abs = sprintf "%.1f", ($w_abs_overlap*100/$lenght_corrected);
              #2 -> 1 coding sequence
              #$overlap21 = sprintf "%.1f", ($overlap*100/$lenght_corrected);
              $overlap21_abs = sprintf "%.1f", ($abs_overlap*100/$lenght_corrected);

              $overlap_JD_abs = sprintf "%.1f", ( $abs_overlap*100/($lenght_corrected+$lenght1-$abs_overlap));
              $w_overlap_JD_abs = sprintf "%.1f", ( $abs_overlap*100/($lenght_corrected+$w_lenght1-$abs_overlap));
              #@list_res = ($gene_id2, $w_overlap12, $w_overlap12_abs, $w_overlap21, $w_overlap21_abs, $overlap12, $overlap12_abs, $overlap21, $overlap21_abs, $w_overlap_JD_abs, $overlap_JD_abs , $proteinName, $descritpion);
              @list_res = ($gene_id2, $w_overlap12_abs, $w_overlap21_abs, $overlap12_abs, $overlap21_abs, $w_overlap_JD_abs, $overlap_JD_abs , $proteinName, $descritpion);
              }
            catch{
              _print1 ("We cannot check the real protein length, let's continue without this one: $prot_tag\n");
              #2 -> 1 whole sequence
              #$w_overlap21 = sprintf "%.1f", ($w_overlap*100/$lenght2);
              #$w_overlap21_abs = sprintf "%.1f", ($w_abs_overlap*100/$lenght2);
              #2 -> 1 coding sequence
              #$overlap21 = sprintf "%.1f", ($overlap*100/$lenght2);
              #$overlap21_abs = sprintf "%.1f", ($abs_overlap*100/$lenght2);
            };
          }
        }
      }
    }
  }
  #print "overlap12 $overlap12 and overlap21 $overlap21\n";
  return \@list_res;
}

# Protein could well align with the genome sequene, but less well with the gene model
# So we re-compute the protein overlap according only to the gene model
sub get_absolute_match{
  my  ($feature, $start, $end)=@_;
  my  $absMatch=0;
  my  $nuc_polish=0;

  # We first need to check that the GAP feature is present among the protein attributes
  if(! $feature->has_tag('Gap')){
    warn "I cannot calculate the absolute match because the tag Gap is absent !\n";
  }
  else{

    # Parse the GAP attribute
    my $gap = $feature->_tag_value('Gap');
    my @gap = split/ /, $gap; # Split gap by space


    my $nuc_left=abs($feature->start - $start+1);
    ############################################
    #Have to re-compute the GAP tag from left  # Case where protein match is longer than the overlap with the gene model (left side)
    if($nuc_left > 0){
      my @newGap;

      foreach my $gap (@gap){

        #compute how long is the piece in nucleotide
        my $nuc = nuc_gap_val($gap);

        #if nucleotide to shrink is over the size of the piece we skip the piece, and compute the size to shrink left
        if ($nuc_left >= $nuc){
          $nuc_left -= $nuc;
          next;
        }
        #if nucleotide to shrink is under the size of the piece we recalculate the piece
        elsif ($nuc_left != 0 and $nuc_left < $nuc){
          my $newGapNucPiece = $nuc - $nuc_left;
          my $letter = substr $gap, 0, 1;

          #has to be modulo3
          my $AAval;
          if($gap =~ /^M/ or $gap =~ /^D/){
            $nuc_polish += $newGapNucPiece % 3;
            $AAval = int($newGapNucPiece/3);
          }
          #avoid case wehre we remove a lot and the piece doesnt exist anymore
          if($AAval != 0){
            my $gap_val_ok = $letter.$AAval;
            push @newGap, $gap_val_ok
          }
        }
        else{
         push @newGap, $gap
        }
      }
      @gap = @newGap;
    }

    my $nuc_right=abs($feature->end - $end+1);
    ############################################
    #Have to re-compute the GAP tag from right  #  Case where protein match is longer than the overlap with the gene model (right side)
    if($nuc_right > 0){
      my @newGap;

      foreach my $gap (@gap){

        #compute how long is the piece in nucleotide
        my $nuc = nuc_gap_val($gap);

        #if nucleotide to shrink is over the size of the piece we skip the piece, and compute the size to shrink left
        if ($nuc_right >= $nuc){
          $nuc_right -= $nuc;
          next;
        }
        #if nucleotide to shrink is under the size of the piece we recalculate the piece
        elsif ($nuc_right != 0 and $nuc_right < $nuc){
          my $newGapNucPiece = $nuc - $nuc_right;
          my $letter = substr $gap, 0, 1;

          #has to be modulo3
          my $AAval;
          if($gap =~ /^M/ or $gap =~ /^D/){
            $nuc_polish += $newGapNucPiece % 3;
            $AAval = int($newGapNucPiece/3);
          }
          #avoid case wehre we remove a lot and the piece doesnt exist anymore
          if($AAval != 0){
            my $gap_val_ok = $letter.$AAval;
            push @newGap, $gap_val_ok
          }
        }
        else{
         push @newGap, $gap
        }
      }
      @gap = @newGap;
    }

    my ($match_size, $nuc_polish) = calcul_match_gap(\@gap, $nuc_polish);

    $absMatch += $match_size;
    #print "match_size = $match_size\n";
    #my $plus= int(abs($feature->end - $end)/3);
    #my $modPlus = int(abs($feature->end - $end) % 3);
    #print "modPlus=".$modMinus." plus".$minus."\n";

  }

  return $absMatch*3+$nuc_polish;
}

#@Output : 2 , Match in AA, nuc_polish in nucleotide
sub  calcul_match_gap{
  my ($gap, $nuc_polish)=@_;

  my $match_size=0;

  foreach my $gap (@{$gap}){

        #get value every time it was a match
        if($gap =~ /^M/){ #MATCH - M1 in a protein space is actually an amino acid match (matches 3 bp in nucleotide space)
          substr($gap, 0, 1) = "";
          $match_size += $gap;
        }
        elsif($gap =~ /^D/){ # deletion = insert a gap into the target (delete from reference) - D1 is an amino acid deletion (3bp in nucleotide space)
          substr($gap, 0, 1) = "";
          $nuc_polish -= $gap;
        }
        elsif($gap =~ /^I/){ # insert a gap into the reference sequence - I1 is an amino acid insertion (3bp in nucleotide space)
          substr($gap, 0, 1) = "";
          $nuc_polish -= ($gap*3);
        }
        elsif($gap =~ /^R/){# frameshift reverse in the reference sequence - F and R therefore allow for single bp movement either to the left or right within amino acid space. Sometime this happens in Exonerate where it appears as a slightly shifted codon (codons look stacked ), but it also happens when an amino acid is split across a splice site (1st part of a codon is on one exon and second part on the next exon).
          substr($gap, 0, 1) = "";
          $nuc_polish -= $gap;
        }
        elsif($gap =~ /^F/){# frameshift forward in the reference sequence
          substr($gap, 0, 1) = "";
          $nuc_polish -= $gap;
        }
        else{
          warn "Cannot interpret this CIGAR substring: $gap !\n";
        }
  }

  return $match_size, $nuc_polish;

}

sub nuc_gap_val{
  my ($gap) = @_;

  my $nuc=0;
  if($gap =~ /^M/){ #MATCH - M1 in a protein space is actually an amino acid match (matches 3 bp in nucleotide space)
    $nuc = substr $gap, 1;
    $nuc *=3;
  }
  elsif($gap =~ /^D/){  # deletion = insert a gap into the target (delete from reference) - D1 is an amino acid deletion (3bp in nucleotide space)
    $nuc = substr $gap, 1;
  }
  elsif($gap =~ /^I/){ # insert a gap into the reference sequence - I1 is an amino acid insertion (3bp in nucleotide space)
    $nuc = substr $gap, 1;
    $nuc *=3;
  }
  elsif($gap =~ /^F/){ # frameshift forward in the reference sequence - F and R therefore allow for single bp movement either to the left or right within amino acid space. Sometime this happens in Exonerate where it appears as a slightly shifted codon (codons look stacked ), but it also happens when an amino acid is split across a splice site (1st part of a codon is on one exon and second part on the next exon).
    $nuc = substr $gap, 1;
  }
  elsif($gap =~ /^R/){ # frameshift reverse in the reference sequence -
    $nuc = substr $gap, 1;
  }
  else{
    warn "Cannot interpret this CIGAR substring: $gap !\n";
  }

  return $nuc;
}

sub _print1{
  my ($message, $optionType) = @_;

  $outputTab[$optionType]->print($message) if defined ($optionType);
  dual_print1 $message;

}

sub _print2{
  my ($message, $optionType) = @_;

  $outputTab[$optionType]->print($message) if defined ($optionType && $CONFIG->{verbose} >= 2);
  dual_print2 $message;

}
__END__

=head1 NAME

agat_sp_load_function_from_protein_align.pl

=head1 DESCRIPTION

The script takes an annotation in gff format, a protein alignment in gff format and a protein fasta file as input. It checks if protein alignement overlap gene models, and will load the gene name and/or the function to the gene model according to the user requirements.
The script applies the following steps:
For each gene model structure it take the proteins aligned against, and sort them by an overlaping score. The best coming first.
Then it filters them by applying the overlaping score threshold.
1) If you activated the PE and the species filtering, we will try to find the best protein that follows the defined requirement.
2.1) If you activated the PE filtering or the precedent filtering (1) didn't succeed, we take the best protein according to the PE requirement.
2.2) If you activated the species filtering or the precedent filtering (1) didn't succeed, we take the best protein according to the list of prioritized species defined.
3) If no option or the precedent filtering (1,2.1,2.2)didn't succeed, the best protein will be selected.
You can flip the 2.1 and 2.2 test using the priority option.


=head1 SYNOPSIS

    agat_sp_load_function_from_protein_align.pl -a annotation.gff --pgff protein.gff --pfasta protein.fasta [ -o outfile ]
    agat_sp_load_function_from_protein_align.pl --help

=head1 OPTIONS

=over 8

=item B<-a> or B<--annotation>

Input gtf/gff file of an annotation.

=item B<--pgff>

Input gff file of aligned proteins.

=item B<--pfasta>

Input protein fasta file where the extra information will be retrieved for each aligned protein.

=item B<-m> or B<--method>

Rule to apply to lift function when a protein map properly.
1) replace  => replace or add the product and Name attribute's values.
2) complete => add the product and Name attribute's values only if doesn't exist.
3) add      => add the lfp_product and lfp_name attributes with the corresponding values

=item B<--value>, B<--threshold> or B<-t>

Gene mapping percentage over which a gene must be reported. By default the value is 50.

=item B<-w>

Compute the overlap score based on the whole annotation sequence. By default we use only the coding sequence part.

=item B<--pe>

Protein existence value. We will take the best overlap score protein according to the PE expected
1. Experimental evidence at protein level
2. Experimental evidence at transcript level
3. Protein inferred from homology
4. Protein predicted
5. Protein uncertain

=item B<--test>

Test to apply (> < = >= <=). If you us one of these two character >, <, please don't forget to quote you parameter liket that "<=". Else your terminal will complain.

=item B<--sp>

Species, between the set of the best protein aligned we try first to take the one that follow the species prioritization defined. There is a default one, but you can define you own (quoted and coma separated value)like that: "mus Musculus, Homo Sapiens" from the most important to the less important. In that case Mus will be taken first even if a better overlaping one exist for human.
If none of them is found we take anyway the best overlapping one.

=item B<-p> or B<--priority>

By default the priority is PE test before species test when both are applied. You can flip these two test by activating this option like this: -p species

=item B<-v>

Be verbose.

=item B<-o> , B<--output> or B<--out>

Output GFF file.  If no output file is specified, the output will be
written to STDOUT.

=item B<-thread>, B<threads>, B<cpu>, B<cpus>, B<core>, B<cores>, B<job> or B<jobs>

Integer - Number of parallel processes to use for file input parsing (via forking).

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
