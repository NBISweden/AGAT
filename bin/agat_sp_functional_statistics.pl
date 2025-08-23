#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use IO::File;
use List::MoreUtils qw(uniq);
use AGAT::AGAT;

my $header = get_agat_header();
my @copyARGV = @ARGV;
my ( $opt, $usage, $config ) = AGAT::AGAT::describe_script_options(
    $header,
    [ 'gff|f|ref|reffile=s', 'Input reference gff file', { required => 1 } ],
    [ 'g|gs=i', 'Genome size', { callbacks => { positive => sub { shift() > 0 or die 'Genome size must be positive' } } } ],
);

my $gff           = $opt->gff;
my $opt_genomeSize = $opt->g;
my $opt_output    = $config->{output} // 'output_functional_statistics';
$config->{output} = $opt_output;
my $opt_verbose   = $config->{verbose};

my $log;
if ( my $log_name = $config->{log_path} ) {
    open( $log, '>', $log_name )
      or die "Can not open $log_name for printing: $!";
    dual_print( $log, $header, 0 );
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    PARAMS    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

if ( -f $opt_output ) {
  dual_print( $log, "Cannot create a directory with the name $opt_output because a file with this name already exists.\n", 0 );
  warn "Cannot create a directory with the name $opt_output because a file with this name already exists.\n"
    if $opt_verbose;
  exit();
}
if ( -d $opt_output ) {
  dual_print( $log, "The output directory choosen already exists. Please give me another Name.\n", 0 );
  warn "The output directory choosen already exists. Please give me another Name.\n"
    if $opt_verbose;
  exit();
}
mkdir $opt_output;

my $stat_file = $opt_output."/stat_features.txt";
my $stat_out = prepare_fileout($stat_file);

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>     MAIN     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

######################
### Parse GFF input #
dual_print( $log, "Reading file $gff\n", $opt_verbose );
my ($hash_omniscient, $hash_mRNAGeneLink) = slurp_gff3_file_JD({ input => $gff,
                                                                 config => $config
                                                              });
dual_print( $log, "Parsing Finished\n", $opt_verbose );
### END Parse GFF input #
#########################

###############################################################
### Print Statistics structural first
###############################################################

dual_print( $log, "Compute statistics\n", $opt_verbose );
print_omniscient_statistics ({ input  => $hash_omniscient,
															 genome => $opt_genomeSize,
															 output => $stat_out
														 });

###############################################################
### Print Statistics function
###############################################################:
my %link_tags; # list children tag to parent tags
my %hash_info;
my %dbxref_db;
## Fill the hash
foreach my $tag_l1 (sort keys %{$hash_omniscient->{'level1'}}){
  foreach my $id_l1 (sort keys %{$hash_omniscient->{'level1'}{$tag_l1}}){
    my $chimerel1l2 = $tag_l1;

    foreach my $tag_l2 (sort keys %{$hash_omniscient->{'level2'}}){
      if ( exists_keys( $hash_omniscient,( 'level2',$tag_l2, $id_l1 ) ) ){
        $chimerel1l2=$tag_l1 ."@".$tag_l2;
        $link_tags{$chimerel1l2}{"l2"}{$tag_l2}++;

        foreach my $feature_l2 (@{$hash_omniscient->{'level2'}{$tag_l2}{$id_l1}}){
          my $level2_ID = lc($feature_l2->_tag_value('ID'));
          analyse_feature(\%hash_info, $chimerel1l2, $feature_l2);

          foreach my $tag_l3 (sort keys %{$hash_omniscient->{'level3'}} ){
            if ( exists_keys( $hash_omniscient, ( 'level3',$tag_l3, $level2_ID ) ) ){
              $link_tags{$chimerel1l2}{"l3"}{$tag_l3}++;
              foreach my $feature_l3 ( @{$hash_omniscient->{'level3'}{$tag_l3}{$level2_ID}} ) {
                analyse_feature(\%hash_info, $chimerel1l2, $feature_l3);
              }
            }
          }
        }
        # Analyse level1 with chimere name used with l2
        $link_tags{$chimerel1l2}{"l1"}=$tag_l1;
        my $feature_l1 = $hash_omniscient->{'level1'}{$tag_l1}{$id_l1};
        analyse_feature(\%hash_info, $chimerel1l2, $feature_l1);
      }
    }

    if ($chimerel1l2 eq $tag_l1 ){ # This is not a chimere, means no l2 linked to this l1, so l1 has not been analysed yet.
      $link_tags{$chimerel1l2}{"l1"}=$tag_l1;
      my $feature_l1 = $hash_omniscient->{'level1'}{$tag_l1}{$id_l1};
      analyse_feature(\%hash_info, $tag_l1, $feature_l1);
    }
  }
}

## Resume the info
foreach my $chimerel1l2 (sort keys %link_tags){
  # get level1 name
  my $tag_l1 = $link_tags{$chimerel1l2}{"l1"};
  #print result per type within dedicated file when output provided
  my $type_output = $opt_output."/".$chimerel1l2;
  mkdir $type_output;
  
  my $table_file = $type_output."/table_per_feature_type.txt";
  my $table_file_out = prepare_fileout($table_file);

  print $table_file_out  "\nFunctional info $chimerel1l2 records:\n";
  table_info_l1l3(\%hash_info, $chimerel1l2, $tag_l1, $table_file_out, $type_output);

  foreach my $tag_l2 ( sort keys %{$link_tags{$chimerel1l2}{"l2"}}){
    table_info_l2(\%hash_info, $chimerel1l2, $tag_l1, $tag_l2, $table_file_out, $type_output);

    # Get tag l3
    foreach my $tag_l3 ( sort keys %{$link_tags{$chimerel1l2}{"l3"}} ){
      table_info_l1l3(\%hash_info, $chimerel1l2, $tag_l3, $table_file_out, $type_output);     
    }
  }
}

# END STATISTICS #
##################
dual_print( $log, "Result available in <$opt_output>. Bye Bye.\n", $opt_verbose );
close $log if $log;
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

# Print table for l1
sub table_info_l1l3{
  my ($hash_info, $chimerel1l2, $tag, $table_file_out, $type_output)=@_;

  my $lineB = "_____________________________________________________________________________";
  my $stringPrint=undef;

  my $stringPrintFooter=undef;
  my $nb_up_ft = keys %{$hash_info{$chimerel1l2}{'type_ft'}{$tag}{"id"}};
  $stringPrintFooter .= "Nb $tag = $nb_up_ft\n";

  # --- HEADER ---
  $stringPrint .= " ".$lineB."\n";
  $stringPrint .= "|".sizedPrint(" ",25)."|".sizedPrint("Nb holded by",25)."|".sizedPrint("Nb $tag",25)."|\n";
  $stringPrint .= "|".sizedPrint(" ",25)."|".sizedPrint("$tag",25)."|".sizedPrint("holding it",25)."|\n";
  $stringPrint .= "|".$lineB."|\n";
  
  # --- general details ---
  foreach my $type ("name", "product", "description", "ontology_term", "dbxref"){

    my $term_holded_by_l1 = 0;
    foreach my $id (keys %{$hash_info->{$chimerel1l2}{$type}{$tag}{"id"}} ) {
      $term_holded_by_l1 +=  $hash_info->{$chimerel1l2}{$type}{$tag}{"id"}{$id};
    }
    my $l1_holding_term = keys %{$hash_info->{$chimerel1l2}{$type}{$tag}{"id"}};
    $stringPrint .= "|".sizedPrint("$type",25)."|".sizedPrint($term_holded_by_l1,25)."|".sizedPrint($l1_holding_term,25)."|\n";
    $stringPrint .= "|".$lineB."|\n";

    # Add information about feature without of ...
    my $nb_up_ft_withoutFunction = $nb_up_ft - $l1_holding_term;
    $stringPrintFooter .= "Nb $tag with <$type> attribute = $l1_holding_term\n";
    $stringPrintFooter .= "Nb $tag without <$type> attribute = $nb_up_ft_withoutFunction\n";

    # --- Print ID \t function file ---
    if($type ne "dbxref" ){
      my $type_output_ok = $type_output."/$type";

      if (exists_keys( $hash_info, ($chimerel1l2, $type, $tag, 'value') ) ){
        mkdir $type_output_ok if (! -d $type_output_ok);
        my $type_file = $type_output_ok."/$type"."_"."$tag.txt";
        my $type_file_out = prepare_fileout($type_file);

        foreach my $id (sort keys %{$hash_info->{$chimerel1l2}{$type}{$tag}{"id"}} ) {
          foreach my $value (sort keys %{$hash_info->{$chimerel1l2}{$type}{$tag}{"value"}{$id}} ) {
            print $type_file_out "$id\t$value\n";
          }
        }
      }
    }
  }

  # --- DBXREF details ---
  # loop over all dbxref db seen
  foreach my $db (sort keys %dbxref_db){
    my $total_db_l1 = 0;
    my $db_l1 = 0;
    if ( exists_keys( $hash_info, ( $chimerel1l2, 'dbxref', $tag, "value", $db, "id" ) ) ){

      
      my $type_output_ok = $type_output."/dbxref";
      mkdir $type_output_ok;
      my $db_file = $type_output_ok."/$db"."_"."$tag.txt";
      my $db_out = prepare_fileout($db_file);

      $db_l1 = keys %{$hash_info->{$chimerel1l2}{'dbxref'}{$tag}{"value"}{$db}{"id"}} ;
      foreach my $id (sort keys %{$hash_info->{$chimerel1l2}{'dbxref'}{$tag}{"value"}{$db}{"id"}}){
        $total_db_l1 += $hash_info->{$chimerel1l2}{'dbxref'}{$tag}{"value"}{$db}{"id"}{$id};
        foreach my $value (sort keys %{$hash_info->{$chimerel1l2}{'dbxref'}{$tag}{"value"}{$db}{"value"}{$id}}){
          print $db_out "$id\t$value\n";
        }
      }
    }
    $stringPrint .= "|".sizedPrint("dbxref:$db",25)."|".sizedPrint($total_db_l1,25)."|".sizedPrint($db_l1,25)."|\n";
    $stringPrint .= "|".$lineB."|\n";

    # Add information about feature without of ...
    my $nb_up_ft_withoutFunction = $nb_up_ft - $db_l1;
    $stringPrintFooter .= "Nb $tag with <$db> dbxref = $db_l1\n";
    $stringPrintFooter .= "Nb $tag without <$db> dbxref = $nb_up_ft_withoutFunction\n";
  }

  # Add info about feature without
  $stringPrint .= $stringPrintFooter;
  #Print results
  print $table_file_out $stringPrint;
}


# Print table for l2
sub table_info_l2{
  my ($hash_info, $chimerel1l2, $tag_l1, $tag_l2, $table_file_out, $type_output)=@_;
  
  my $lineB=       "_______________________________________________________________________________________________________";
  my $stringPrint=undef;
  
  my $stringPrintFooter=undef;
  my $nb_l2_ft = keys %{$hash_info{$chimerel1l2}{'type_ft'}{$tag_l2}{"id"}};
  my $nb_l1_ft = keys %{$hash_info->{$chimerel1l2}{'type_ft'}{$tag_l1}{"id"}};
  $stringPrintFooter .= "Nb $tag_l1 = $nb_l1_ft\n";
  $stringPrintFooter .= "Nb $tag_l2 = $nb_l2_ft\n";

  # --- HEADER ---
  $stringPrint .= " ".$lineB."\n";
  $stringPrint .= "|".sizedPrint(" ",25)."|".sizedPrint("Nb holded by",25)."|".sizedPrint("Nb $tag_l2",25)."|".sizedPrint("Nb $tag_l1 with",25)."|\n";
  $stringPrint .= "|".sizedPrint(" ",25)."|".sizedPrint("$tag_l2",25)."|".sizedPrint("holding it",25)."|".sizedPrint("$tag_l2 holding it",25)."|\n";
  $stringPrint .= "|".$lineB."|\n";
  
  # --- general details ---
  foreach my $type ("name", "product", "description", "ontology_term", "dbxref"){

    my $term_holded_by_l2 = 0;
    foreach my $id (keys %{$hash_info->{$chimerel1l2}{$type}{$tag_l2}{"id"}} ) {
      $term_holded_by_l2 +=  $hash_info->{$chimerel1l2}{$type}{$tag_l2}{"id"}{$id};
    }
    my $l2_holding_term = keys %{$hash_info->{$chimerel1l2}{$type}{$tag_l2}{"id"}};
    my $l1_from_l2_holding_term = keys %{$hash_info->{$chimerel1l2}{$type}{$tag_l2}{"parent"}};
    $stringPrint .= "|".sizedPrint("$type",25)."|".sizedPrint($term_holded_by_l2,25)."|".sizedPrint($l2_holding_term,25)."|".sizedPrint($l1_from_l2_holding_term,25)."|\n";
    $stringPrint .= "|".$lineB."|\n";

    # Add information about feature without of ...
    my $nb_up_ft_withoutFunction = $nb_l1_ft - $l1_from_l2_holding_term;
    $stringPrintFooter .= "Nb $tag_l1 with <$type> attribute = $l1_from_l2_holding_term\n";
    $stringPrintFooter .= "Nb $tag_l1 without <$type> attribute = $nb_up_ft_withoutFunction\n";
    $nb_up_ft_withoutFunction = $nb_l2_ft - $l2_holding_term;
    $stringPrintFooter .= "Nb $tag_l2 with <$type> attribute = $l2_holding_term\n";
    $stringPrintFooter .= "Nb $tag_l2 without <$type> attribute = $nb_up_ft_withoutFunction\n";

    # --- Print ID \t function file ---
    if($type ne "dbxref" ){
      if (exists_keys( $hash_info, ($chimerel1l2, $type, $tag_l2, 'value') ) ){
        my $type_output_ok = $type_output."/$type";
        mkdir $type_output_ok if (! -d $type_output_ok);
        my $type_file = $type_output_ok."/$type"."_"."$tag_l2.txt";
        my $type_file_out = prepare_fileout($type_file);
        foreach my $id (sort keys %{$hash_info->{$chimerel1l2}{$type}{$tag_l2}{"id"}} ) {
          foreach my $value (sort keys %{$hash_info->{$chimerel1l2}{$type}{$tag_l2}{"value"}{$id}} ) {
            print $type_file_out "$id\t$value\n";
          }
        }
      }
    }
  }

  # --- DBXREF details ---
  # loop over all dbxref db seen
  foreach my $db (sort keys %dbxref_db){
    my $db_l2 = 0;
    my $total_db_l2 = 0;
    my $db_l1_from_l2 = 0;
    my $db_l1 = 0;

    if ( exists_keys( $hash_info, ( $chimerel1l2, 'dbxref', $tag_l2, "value", $db, "id" ) ) ){
      $db_l2 = keys %{$hash_info->{$chimerel1l2}{'dbxref'}{$tag_l2}{"value"}{$db}{"id"}} ;

      my $type_output_ok = $type_output."/dbxref";
      mkdir $type_output_ok if (! -d $type_output_ok);
      my $db_file = $type_output_ok."/$db"."_"."$tag_l2.txt";
      my $db_out = prepare_fileout($db_file);

      foreach my $id (sort keys %{$hash_info->{$chimerel1l2}{'dbxref'}{$tag_l2}{"value"}{$db}{"id"}}){
        $total_db_l2 += $hash_info->{$chimerel1l2}{'dbxref'}{$tag_l2}{"value"}{$db}{"id"}{$id};
        foreach my $value (sort keys %{$hash_info->{$chimerel1l2}{'dbxref'}{$tag_l2}{"value"}{$db}{"value"}{$id}}){
          print $db_out "$id\t$value\n";
        }
      }
    }
    if ( exists_keys( $hash_info, ( $chimerel1l2, 'dbxref', $tag_l2, "value", $db, "parent" ) ) ){
      $db_l1_from_l2 = keys %{$hash_info->{$chimerel1l2}{'dbxref'}{$tag_l2}{"value"}{$db}{"parent"}} ;
    }
    $stringPrint .= "|".sizedPrint("dbxref:$db",25)."|".sizedPrint($total_db_l2,25)."|".sizedPrint($db_l2,25)."|".sizedPrint($db_l1_from_l2,25)."|\n";
    $stringPrint .= "|".$lineB."|\n";

    # Add information about feature without of ...
    my $nb_up_ft_withoutFunction = $nb_l1_ft - $db_l1_from_l2;
    $stringPrintFooter .= "Nb $tag_l1 with $db dbxref = $db_l1_from_l2\n";
    $stringPrintFooter .= "Nb $tag_l1 without $db dbxref = $nb_up_ft_withoutFunction\n";
    $nb_up_ft_withoutFunction = $nb_l2_ft - $db_l2;
    $stringPrintFooter .= "Nb $tag_l2 with $db dbxref = $db_l2\n";
    $stringPrintFooter .= "Nb $tag_l2 without $db dbxref = $nb_up_ft_withoutFunction\n";
  }

  # Add info about feature without
  $stringPrint .= $stringPrintFooter;
  #Print results
  print $table_file_out $stringPrint;
}

sub analyse_feature{
  my ($hash_info, $tag_l1, $feature)=@_;

  my $id = $feature->_tag_value('ID');
  my $type_ft = lc($feature->primary_tag);
  my $parent = undef;
  if($feature->has_tag('Parent') ){
    $parent = $feature->_tag_value('Parent');
  }

  # General info
  $hash_info{$tag_l1}{'type_ft'}{$type_ft}{"id"}{$id}++ ;
  $hash_info{$tag_l1}{'type_ft'}{$type_ft}{"parent"}{$parent}++ if ($parent);

  #Check For NAME
  if($feature->has_tag('Name') ){
    my @values = $feature->get_tag_values('Name');
    foreach my $value (@values){
      $hash_info{$tag_l1}{'name'}{$type_ft}{"value"}{$id}{$value}++ ;
      $hash_info{$tag_l1}{'name'}{$type_ft}{"id"}{$id}++ ;
      $hash_info{$tag_l1}{'name'}{$type_ft}{"parent"}{$parent}++ if ($parent);
    }
  }

  #Check For product /!\ can have "hypothetical protein" values
  if($feature->has_tag('product') ){
    my @values = $feature->get_tag_values('product');
    foreach my $value (@values){
      $hash_info{$tag_l1}{'product'}{$type_ft}{"value"}{$id}{$value}++ ;
      $hash_info{$tag_l1}{'product'}{$type_ft}{"id"}{$id}++ ;
      $hash_info{$tag_l1}{'product'}{$type_ft}{"parent"}{$parent}++ if ($parent);
    }
  }

  #Check For description /!\ can have "hypothetical protein" values
  if($feature->has_tag('description') ){
    my @values = $feature->get_tag_values('description');
    foreach my $value (@values){
      $hash_info{$tag_l1}{'description'}{$type_ft}{"value"}{$id}{$value}++ ;
      $hash_info{$tag_l1}{'description'}{$type_ft}{"id"}{$id}++ ;
      $hash_info{$tag_l1}{'description'}{$type_ft}{"parent"}{$parent}++ if ($parent);
    }
  }

  #Check For Ontology_term
  if($feature->has_tag('Ontology_term') ){
    my @values = $feature->get_tag_values('Ontology_term');
    foreach my $tuple (@values){
      my ($type,$value) = split /:/,$tuple;
      $hash_info{$tag_l1}{'ontology_term'}{$type_ft}{"value"}{$id}{$value}++ ;
      $hash_info{$tag_l1}{'ontology_term'}{$type_ft}{"id"}{$id}++ ;
      $hash_info{$tag_l1}{'ontology_term'}{$type_ft}{"parent"}{$parent}++ if ($parent);
      #print "l2 has tag ontology_term with value:".$value."\n";
    }
  }

  #Check For Dbxref
  if($feature->has_tag('Dbxref') ){
    my @values = $feature->get_tag_values('Dbxref');
    foreach my $tuple (@values){
      my ($type,$value) = split /:/,$tuple;
      $dbxref_db{$type}++; 
      $hash_info{$tag_l1}{'dbxref'}{$type_ft}{"value"}{$type}{"value"}{$id}{$value}++ ;
      $hash_info{$tag_l1}{'dbxref'}{$type_ft}{"value"}{$type}{"id"}{$id}++;
      $hash_info{$tag_l1}{'dbxref'}{$type_ft}{"value"}{$type}{"parent"}{$parent}++ if ($parent);
      $hash_info{$tag_l1}{'dbxref'}{$type_ft}{"id"}{$id}++ ;
      $hash_info{$tag_l1}{'dbxref'}{$type_ft}{"parent"}{$parent}++ if ($parent);
    }
  }
  elsif($feature->has_tag('db_xref') ){
    my @values = $feature->get_tag_values('db_xref');
    foreach my $tuple (@values){
      my ($type,$value) = split /:/,$tuple;
      $dbxref_db{$type}++; 
      $hash_info{$tag_l1}{'dbxref'}{$type_ft}{"value"}{$type}{"value"}{$id}{$value}++ ;
      $hash_info{$tag_l1}{'dbxref'}{$type_ft}{"value"}{$type}{"id"}{$id}++;
      $hash_info{$tag_l1}{'dbxref'}{$type_ft}{"value"}{$type}{"parent"}{$parent}++ if ($parent);
      $hash_info{$tag_l1}{'dbxref'}{$type_ft}{"id"}{$id}++ ;
      $hash_info{$tag_l1}{'dbxref'}{$type_ft}{"parent"}{$parent}++ if ($parent);
    }
  }
}

__END__

=head1 NAME

agat_sp_functional_statistics.pl

=head1 DESCRIPTION

The script aims to summerize functional information stored in a GXF/GTF file.

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

Folder where will be written the results. [Default output_functional_statistics]

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
