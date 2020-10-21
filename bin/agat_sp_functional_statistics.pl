#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Pod::Usage;
use Statistics::R;
use IO::File;
use List::MoreUtils qw(uniq);
use Bio::Tools::GFF;
use AGAT::Omniscient;

my $header = get_agat_header();
my $gff = undef;
my $opt_output = undef;
my $opt_genomeSize = undef;
my $opt_help= 0;

if ( !GetOptions(
    "help|h" => \$opt_help,
    'g|gs=s' => \$opt_genomeSize,
    'o|output=s'      => \$opt_output,
    "gff|f=s" => \$gff))

{
    pod2usage( { -message => "Failed to parse command line",
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

#### IN / OUT
my $out;
if ($opt_output) {

  if (-f $opt_output){
      print "Cannot create a directory with the name $opt_output because a file with this name already exists.\n";exit();
  }
  if (-d $opt_output){
      print "The output directory choosen already exists. Please give me another Name.\n";exit();
  }
  mkdir $opt_output;

  $out=IO::File->new(">".$opt_output."/report.txt" ) or croak( sprintf( "Can not open '%s' for writing %s", $opt_output."/report.txt", $! ));
  }
else{
  $out = IO::File->new();
  $out->fdopen( fileno(STDOUT), 'w' );
}


                #####################
                #     MAIN          #
                #####################

######################
### Parse GFF input #
print "Reading file $gff\n";
my ($hash_omniscient, $hash_mRNAGeneLink) = slurp_gff3_file_JD({ input => $gff
                                                              });
print "Parsing Finished\n";
### END Parse GFF input #
#########################

###############################################################
### Print Statistics structural first
###############################################################

print "Compute statistics\n";
print_omniscient_statistics ({ input => $hash_omniscient,
															 genome => $opt_genomeSize,
															 output => $out
														 });

###############################################################
### Print Statistics function
###############################################################:
my %names_l1;
my $name_l1_nb=undef;
my %names_l2;
my $name_l2_nb=undef;
my %products;
my $product_l2_nb=undef;
my %descriptions;
my $description_l2_nb=undef;
my %ontology_terms;
my $ontology_term_l2_nb=undef;

my %DB_omni_mrna;
my %DB_omni_gene;

my $nbmRNAwithFunction = 0;
my $nbGeneWithFunction = 0;
my $nbGeneWithProduct = 0;
my $nbGeneWithDescription = 0;
my $total_nb_l1 = 0;
my $total_nb_l2 = 0;

  #################
  # == LEVEL 1 == #
  #################
foreach my $primary_tag_key_level1 (keys %{$hash_omniscient->{'level1'}}){ # primary_tag_key_level1 = gene or repeat etc...
  foreach my $gene_id_tag_key (keys %{$hash_omniscient->{'level1'}{$primary_tag_key_level1}}){
    $total_nb_l1++;
    my $l1_has_function=undef;
    my $l1_has_product=undef;
    my $l1_has_description=undef;
    my $gene_feature=$hash_omniscient->{'level1'}{$primary_tag_key_level1}{$gene_id_tag_key};
    my $id_gene=$gene_feature->_tag_value('ID');

    #Check For NAME
    if($gene_feature->has_tag('Name') ){
      my $value = $gene_feature->_tag_value('Name');
      $names_l1{$value}++;
      $name_l1_nb++;
      #print "l1 has tag name with value:".$value."\n";
    }



    foreach my $primary_tag_key_level2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
      if ( exists_keys( $hash_omniscient, ('level2', $primary_tag_key_level2, $gene_id_tag_key) ) ){
        foreach my $level2_feature ( @{$hash_omniscient->{'level2'}{$primary_tag_key_level2}{$gene_id_tag_key}}) {
          $total_nb_l2++;
          my $l2_has_function=undef;
          my $l2_has_product=undef;
          my $l2_has_description=undef;
          my $id_mrna=$level2_feature->_tag_value('ID');

          #Check For NAME
          if($level2_feature->has_tag('Name') ){
            my $value = $level2_feature->_tag_value('Name');
            $names_l2{$value}++;
            $name_l2_nb++;
            #print "l2 has tag name with value:".$value."\n";
          }

          #Check For product
          if($level2_feature->has_tag('product') ){
            my $value = $level2_feature->_tag_value('product');
            if ($value ne "hypothetical protein"){
                $products{$value}++;
                $product_l2_nb++;
                $l2_has_product=1;
                #print "l2 has tag product with value:".$value."\n";
            }
          }


          #Check For description
          if($level2_feature->has_tag('description') ){
            my $value = $level2_feature->_tag_value('description');
            if ($value ne "hypothetical protein"){
                $descriptions{$value}++;
                $description_l2_nb++;
                $l2_has_description=1;
                #print "l2 has tag descritpion with value:".$value."\n";
            }
          }

          #Check For Ontology_term
          if($level2_feature->has_tag('Ontology_term') ){
            my @values = $level2_feature->get_tag_values('Ontology_term');
            foreach my $tuple (@values){
              my ($type,$value) = split /:/,$tuple;
              $ontology_terms{$value}++;
              $ontology_term_l2_nb++;
              $l2_has_function=1;
              push @{$DB_omni_mrna{'Ontology_term'}{$id_mrna}}, $value;
              $DB_omni_gene{'Ontology_term'}{$id_gene}++;
              #print "l2 has tag ontology_term with value:".$value."\n";
            }
          }

          #Check For Dbxref
          if($level2_feature->has_tag('Dbxref') ){
            my @values = $level2_feature->get_tag_values('Dbxref');
            foreach my $tuple (@values){
              my ($type,$value) = split /:/,$tuple;
              push @{$DB_omni_mrna{$type}{$id_mrna}}, $value;
              $DB_omni_gene{$type}{$id_gene}++;
              $l2_has_function=1;
            }
          }
          elsif($level2_feature->has_tag('db_xref') ){
            my @values = $level2_feature->get_tag_values('db_xref');
            foreach my $tuple (@values){
              my ($type,$value) = split /:/,$tuple;
              push @{$DB_omni_mrna{$type}{$id_mrna}}, $value;
              $DB_omni_gene{$type}{$id_gene}++;
              $l2_has_function=1;
            }
          }

          if($l2_has_function){
            $nbmRNAwithFunction++;
            $l1_has_function=1;
          }
          if($l2_has_product){
            $l1_has_product=1;
          }
          if($l2_has_description){
            $l1_has_description=1;
          }
        }
      }
    }
    if($l1_has_function){
      $nbGeneWithFunction++;
    }
    if($l1_has_product){
      $nbGeneWithProduct++;
    }
    if($l1_has_description){
      $nbGeneWithDescription++;
    }
  }
}

#print result per type within dedicated file when output provided
# create streamOutput
if($opt_output){
  foreach my $type (keys %DB_omni_mrna){
    my $ostreamFunct = IO::File->new();
    $ostreamFunct->open( $opt_output."/$type.txt", 'w' ) or
        croak(
            sprintf( "Can not open '%s' for writing %s", $opt_output."/$type.txt", $! )
        );
    foreach my $seq_id (keys %{$DB_omni_mrna{$type}}){
      print $ostreamFunct $seq_id."\t".join( ',', @{$DB_omni_mrna{$type}{$seq_id}} )."\n";
    }
  }
}

my $nbmRNAwithoutFunction= $total_nb_l2 - $nbmRNAwithFunction;
my $nbGeneWithoutFunction= $total_nb_l1 - $nbGeneWithFunction;
my $nbGeneWithoutProduct= $total_nb_l1 - $nbGeneWithProduct;

my $listOfFunction="";
foreach my $funct (sort keys %DB_omni_mrna){
  $listOfFunction.="$funct,";
}
chop $listOfFunction;

# NOW summerize
my $stringPrint=undef;
my $lineB=       "_______________________________________________________________________________________________________";
$stringPrint .= " ".$lineB."\n";
$stringPrint .= "|".sizedPrint(" ",25)."|".sizedPrint("Nb term linked to mRNA",25)."|".sizedPrint("Nb mRNA with term",25)."|".sizedPrint("Nb gene with term",25)."|\n";
$stringPrint .= "|".$lineB."|\n";

foreach my $type (sort keys %DB_omni_mrna){
    my $total_term_mRNA=0;
    foreach my $id_l2 (keys %{$DB_omni_mrna{$type}} ){
      $total_term_mRNA+=scalar @{$DB_omni_mrna{$type}{$id_l2}};
    }
    my $nbmRNA_with_term = keys %{$DB_omni_mrna{$type}};
    my $nbGenewith_term = keys %{$DB_omni_gene{$type}};

    my $mRNA_type =0; #keys %{$mRNAAssociatedToTerm{$type}};
    my $gene_type =0; #keys %{$GeneAssociatedToTerm{$type}};
    $stringPrint .= "|".sizedPrint(" $type",25)."|".sizedPrint($total_term_mRNA,25)."|".sizedPrint($nbmRNA_with_term,25)."|".sizedPrint($nbGenewith_term,25)."|\n|".$lineB."|\n";
  }



$stringPrint .= "\nnb mRNA without Functional annotation ($listOfFunction) = $nbmRNAwithoutFunction (remind: total mRNA = $total_nb_l2)\n".
                  "nb mRNA with Functional annotation ($listOfFunction) = $nbmRNAwithFunction (remind: total mRNA = $total_nb_l2)\n".
                  "nb gene without Functional annotation ($listOfFunction) = $nbGeneWithoutFunction (remind: total gene = $total_nb_l1)\n".
                  "nb gene with Functional annotation ($listOfFunction) = $nbGeneWithFunction (remind: total gene = $total_nb_l1)\n\n";


#-----name------
if ($name_l1_nb){
  $stringPrint .= "We found $name_l1_nb genes with <Name> attribute. (remind: total gene = $total_nb_l1)\n";
}
else{$stringPrint .= "No gene with <Name> attribute found.\n";}
if ($name_l2_nb){
  $stringPrint .= "We found $name_l2_nb mRNAs with <Name> attribute. They probably have the same names as their parent genes. (remind: total mRNA = $total_nb_l2)\n";
}
else{$stringPrint .= "No mRNA with <Name> attribute found.\n";}

#-----description------
if ($nbGeneWithDescription){
   $stringPrint .= "We found $nbGeneWithDescription genes with <description> attribute.\n";
}
else{$stringPrint .= "No gene with <description> attribute found.\n";}
if ($description_l2_nb){
  $stringPrint .= "We have $description_l2_nb mRNAs with <description> attribute.\n";
}
else{$stringPrint .= "No mRNA with <description> attribute found.\n";}

#-----product------
if($nbGeneWithProduct){
  $stringPrint .= "We found $nbGeneWithProduct genes with <product> attribute.\n";
}
else{$stringPrint .= "No gene with <product> attribute found.\n";}
if ($product_l2_nb){
  $stringPrint .= "We have $product_l2_nb mRNAs with <product> attribute.\n";
}
else{$stringPrint .= "No mRNA with <product> attribute found.\n";}


print $out $stringPrint;
# END STATISTICS #
##################
print "Bye Bye.\n";
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

__END__

=head1 NAME

agat_sp_functional_statistics.pl

=head1 DESCRIPTION

The script aims to summerize functional information stored in the file.

=head1 SYNOPSIS

    agat_sp_functional_statistics.pl --gff file.gff  [ -o outfile ]
    agat_sp_functional_statistics.pl --help

=head1 OPTIONS

=over 8

=item B<--gff> or B<-f>

Input GTF/GFF file.

=item B<--gs> or B<-g>

This option inform about the genome size in oder to compute more statistics.
You can give the size in Nucleotide or directly the fasta file.


=item B<--output> or B<-o>

File where will be written the result. If no output file is specified,
the output will be written to STDOUT.

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
