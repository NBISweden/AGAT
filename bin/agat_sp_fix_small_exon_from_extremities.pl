#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Clone 'clone';
use Getopt::Long;
use Pod::Usage;
use List::MoreUtils qw(uniq);
use Bio::DB::Fasta;
use Bio::SeqIO;
use AGAT::AGAT;

my $header = get_agat_header();
my $config;
my $outfile = undef;
my $gff = undef;
my $file_fasta=undef;
my $codonTableId=1;
my $SIZE_OPT=15;
my $verbose = undef;
my $opt_help= 0;

my @copyARGV=@ARGV;
if ( !GetOptions(
    'c|config=s'               => \$config,
    "h|help" => \$opt_help,
    "gff=s" => \$gff,
    "fasta|fa|f=s" => \$file_fasta,
    "table|codon|ct=i" => \$codonTableId,
    "size|s=i" => \$SIZE_OPT,
    "v!" => \$verbose,
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

if ( ! (defined($gff)) or !(defined($file_fasta)) ){
    pod2usage( {
           -message => "$header\nAt least 2 parameter is mandatory:\nInput reference gff file (--gff) and Input fasta file (--fasta)\n\n",
           -verbose => 0,
           -exitval => 1 } );
}

# --- Manage config ---
$config = get_agat_config({config_file_in => $config});

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    PARAMS    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

my $gffout = prepare_gffout($config, $outfile);

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    EXTRA     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

$codonTableId = get_proper_codon_table($codonTableId);
print "Codon table ".$codonTableId." in use. You can change it using --table option.\n";

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>     MAIN     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

######################
### Parse GFF input #
my ($hash_omniscient, $hash_mRNAGeneLink) = slurp_gff3_file_JD({ input => $gff,
                                                                 config => $config
                                                              });
print ("GFF3 file parsed\n");


####################
# index the genome #
my $db = Bio::DB::Fasta->new($file_fasta);
print ("Fasta file parsed\n");

####################

#counters
my $exonCounter=0;
my $mrnaCounter=0;
my $geneCounter=0;


foreach my $primary_tag_key_level1 (keys %{$hash_omniscient->{'level1'}}){ # primary_tag_key_level1 = gene or repeat etc...
  foreach my $gene_id (keys %{$hash_omniscient->{'level1'}{$primary_tag_key_level1}}){

    my $gene_feature = $hash_omniscient->{'level1'}{$primary_tag_key_level1}{$gene_id};
    my $strand = $gene_feature->strand();
    print "gene_id = $gene_id\n" if $verbose;

    foreach my $primary_tag_key_level2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
      if ( exists_keys( $hash_omniscient, ('level2', $primary_tag_key_level2, $gene_id) ) ){
        my $rnaFix=undef;
        foreach my $level2_feature ( @{$hash_omniscient->{'level2'}{$primary_tag_key_level2}{$gene_id}}) {

          # get level2 id
          my $level2_ID = lc($level2_feature->_tag_value('ID'));

          my $exonFix=undef;
          if ( exists_keys( $hash_omniscient, ('level3', 'exon', $level2_ID) ) ){
            my @exon_sorted = sort {$a->start <=> $b->start} @{$hash_omniscient->{'level3'}{'exon'}{$level2_ID}};

            my $number_exon=$#{$hash_omniscient->{'level3'}{'exon'}{$level2_ID}}+1;

            #####################
            #start with left exon
            my $left_exon = $exon_sorted[0];
            my $exon_size = ($left_exon->end - $left_exon->start +1);

            if($exon_size < $SIZE_OPT){

              my $original_exon_start = $left_exon->start;
              my $new_exon_start = $left_exon->start-($SIZE_OPT - $exon_size );

              #modify the exon start
              $left_exon->start($new_exon_start);
              $exonCounter++;
              $exonFix=1;

              print "left_exon start fixed\n" if $verbose;

              #take care of CDS if needed
              if ( exists_keys( $hash_omniscient, ('level3', 'cds', $level2_ID) ) ){
                my @cds_sorted = sort {$a->start <=> $b->start} @{$hash_omniscient->{'level3'}{'cds'}{$level2_ID}};

                #Check if the exon modification could affect the CDS
                if($original_exon_start == $cds_sorted[0]->start()){

                  my $original_cds_start = $original_exon_start;

                  #get the sequence
                  my $sequence = $db->seq( $gene_feature->seq_id() );
                  #get codon table
                  my $codonTable = Bio::Tools::CodonTable->new( -id => $codonTableId);

                  #extract the codon
                  my $this_codon = substr( $sequence, $original_cds_start-1, 3);

                  if($strand eq "+" or $strand == "1"){
                     #Check if it is not terminal codon, otherwise we have to extend the CDS.

                    if(! $codonTable->is_start_codon( $this_codon )){
                      print "first exon plus strand : this is not a start codon\n";exit;
                    }

                  }
                  if($strand eq "-" or $strand == "-1"){
                    #reverse complement
                    my $seqobj = Bio::Seq->new(-seq => $this_codon);
                    $this_codon = $seqobj->revcom()->seq;

                    #Check if it is not terminal codon, otherwise we have to extend the CDS.
                    if(! $codonTable->is_ter_codon( $this_codon )){
                      print "first exon minus strand : this is not a terminal codon\n";exit;
                    }
                  }
                }
              }

            }
            ################
            #then right exon
            if($number_exon > 1){

              my $right_exon =  $exon_sorted[$#exon_sorted];
              my $exon_size = ($right_exon->end - $right_exon->start +1);

              if($exon_size < $SIZE_OPT){

                my $original_exon_end = $right_exon->end;
                my $new_exon_end = $right_exon->end+($SIZE_OPT - $exon_size );

                #modify the exon end
                $right_exon->end($new_exon_end);
                $exonCounter++;
                $exonFix=1;

                print "right_exon end fixed\n" if $verbose;

                #take care of CDS if needed
                if ( exists_keys( $hash_omniscient, ('level3', 'cds', $level2_ID) ) ){
                  my @cds_sorted = sort {$a->start <=> $b->start} @{$hash_omniscient->{'level3'}{'cds'}{$level2_ID}};

                  #Check if the exon modification could affect the CDS
                  if($original_exon_end == $cds_sorted[$#cds_sorted]->end()){

                    my $original_cds_end = $original_exon_end;

                    #get the sequence
                    my $sequence = $db->seq( $gene_feature->seq_id() );
                    #get codon table
                    my $codonTable = Bio::Tools::CodonTable->new( -id => $codonTableId);

                    #extract the codon
                    my $this_codon = substr( $sequence, $original_cds_end-3, 3);

                    if($strand eq "+" or $strand == "1"){
                      print "last plus strand\n" if $verbose;
                       #Check if it is not terminal codon, otherwise we have to extend the CDS.

                      if(! $codonTable->is_ter_codon( $this_codon )){

                        print "last exon plus strand : $this_codon is not a stop codon\n";exit;
                      }

                    }
                    if($strand eq "-" or $strand == "-1"){
                      print "last minus strand\n" if $verbose;

                      #reverse complement
                      my $seqobj = Bio::Seq->new(-seq => $this_codon);
                      $this_codon = $seqobj->revcom()->seq;

                      #Check if it is not terminal codon, otherwise we have to extend the CDS.
                      if(! $codonTable->is_start_codon( $this_codon )){
                        print "last exon minus strand : $this_codon is not a start codon\n";exit;
                      }
                    }
                  }
                }
              }
            }
          }
          if($exonFix){
            $mrnaCounter++;
          }
        }
        if($rnaFix){
          $geneCounter++;
        }
      }
    }
  }
}

check_all_level2_locations( { omiscient => $hash_omniscient } ); # review all the feature L2 to adjust their start and stop according to the extrem start and stop from L3 sub features.
check_all_level1_locations( { omiscient => $hash_omniscient } ); # Check the start and end of level1 feature based on all features level2.

#END
my $string_to_print="usage: $0 @copyARGV\n";
$string_to_print .="Results:\n";
$string_to_print .="nb gene affected: $geneCounter\n";
$string_to_print .="nb rna affected: $mrnaCounter\n";
$string_to_print .="nb exon affected: $exonCounter\n";
print $string_to_print;

print_omniscient( {omniscient => $hash_omniscient, output => $gffout} );

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

agat_fix_small_exon_from_extremities.pl

=head1 DESCRIPTION

The script aims to extend the small exons to make them longer.
When submitting annotation to ENA they expect exon size of 15 nt minimum.
Currently we extend only the exon from extremities, otherwise we risk to break the predicted ORF.
/!\ When we extend an exon and the CDS has to be extended too (because is was a partial CDS), we exit;

=head1 SYNOPSIS

    agat_fix_small_exon_from_extremities.pl -gff infile.gff --fasta genome.fa [ -o outfile ]
    agat_fix_small_exon_from_extremities.pl --help

=head1 OPTIONS

=over 8

=item B<-gff>

Input GTF/GFF file.

=item B<-fa> or B<--fasta>

Genome fasta file
The name of the fasta file containing the genome to work with.

=item B<--ct> or B<--table> or B<--codon>

This option allows specifying the codon table to use - It expects an integer (1 by default = standard)

=item B<--size> or B<-s>

Minimum exon size accepted in nucleotide. All exon below this size will be extended to this size. Default value = 15.

=item B<-o> , B<--output> , B<--out> or B<--outfile>

Output GFF file.  If no output file is specified, the output will be
written to STDOUT.

=item B<-v>

Verbose option, make it easier to follow what is going on for debugging purpose.

=item B<-c> or B<--config>

String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any, 
otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
The --config option gives yo the possibility to use your own AGAT config file (located elsewhere or named differently).

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
