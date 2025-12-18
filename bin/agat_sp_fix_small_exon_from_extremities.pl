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

start_script();
my $header = get_agat_header();
# ---------------------------- OPTIONS ----------------------------
my $outfile = undef;
my $gff = undef;
my $file_fasta=undef;
my $codonTableId=1;
my $SIZE_OPT=15;
my $opt_help= 0;

# OPTION MANAGEMENT: partition @ARGV into shared vs script options via library
my ($shared_argv, $script_argv) = split_argv_shared_vs_script(\@ARGV);

# Parse script-specific options from its own list
my $script_parser = Getopt::Long::Parser->new;
$script_parser->configure('bundling','no_auto_abbrev');
if ( ! $script_parser->getoptionsfromarray(
  $script_argv,
  'h|help!'                 => \$opt_help,
  'gff=s'                   => \$gff,
  'fasta|fa|f=s'            => \$file_fasta,
  'table|codon|ct=i'        => \$codonTableId,
  'size|s=i'                => \$SIZE_OPT,
  'output|out|o=s'  => \$outfile))

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
my $shared_opts = parse_shared_options($shared_argv);
initialize_agat({ config_file_in => ( $shared_opts->{config} ), input => $gff, shared_opts => $shared_opts });

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    PARAMS    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    PARAMS    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

my $gffout = prepare_gffout( $outfile );

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    EXTRA     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# --- Check codon table
$codonTableId = get_proper_codon_table($codonTableId);

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>     MAIN     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

######################
### Parse GFF input #
my ($hash_omniscient) = slurp_gff3_file_JD({ input => $gff });

####################
# index the genome #
my $db = Bio::DB::Fasta->new($file_fasta);
dual_print1 "Fasta file parsed\n";

####################

#counters
my $exonCounter=0;
my $mrnaCounter=0;
my $geneCounter=0;


foreach my $primary_tag_key_level1 (keys %{$hash_omniscient->{'level1'}}){ # primary_tag_key_level1 = gene or repeat etc...
  foreach my $gene_id (keys %{$hash_omniscient->{'level1'}{$primary_tag_key_level1}}){

    my $gene_feature = $hash_omniscient->{'level1'}{$primary_tag_key_level1}{$gene_id};
    my $strand = $gene_feature->strand();
    dual_print2 "gene_id = $gene_id\n";

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

              dual_print2 "left_exon start fixed\n";

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
                      die "first exon plus strand : this is not a start codon\n";
                    }

                  }
                  if($strand eq "-" or $strand == "-1"){
                    #reverse complement
                    my $seqobj = Bio::Seq->new(-seq => $this_codon);
                    $this_codon = $seqobj->revcom()->seq;

                    #Check if it is not terminal codon, otherwise we have to extend the CDS.
                    if(! $codonTable->is_ter_codon( $this_codon )){
                      die "first exon minus strand : this is not a terminal codon\n";
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

                dual_print2 "right_exon end fixed\n";

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
                      dual_print2 "last plus strand\n";
                       #Check if it is not terminal codon, otherwise we have to extend the CDS.

                      if(! $codonTable->is_ter_codon( $this_codon )){

                        die "last exon plus strand : $this_codon is not a stop codon\n";
                      }

                    }
                    if($strand eq "-" or $strand == "-1"){
                      dual_print2 "last minus strand\n";

                      #reverse complement
                      my $seqobj = Bio::Seq->new(-seq => $this_codon);
                      $this_codon = $seqobj->revcom()->seq;

                      #Check if it is not terminal codon, otherwise we have to extend the CDS.
                      if(! $codonTable->is_start_codon( $this_codon )){
                        die "last exon minus strand : $this_codon is not a start codon\n";
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

check_all_level2_locations( { omniscient => $hash_omniscient } ); # review all the feature L2 to adjust their start and stop according to the extrem start and stop from L3 sub features.
check_all_level1_locations( { omniscient => $hash_omniscient } ); # Check the start and end of level1 feature based on all features level2.

#END
my $string_to_print .="Results:\n";
$string_to_print .="nb gene affected: $geneCounter\n";
$string_to_print .="nb rna affected: $mrnaCounter\n";
$string_to_print .="nb exon affected: $exonCounter\n";
dual_print1 "$string_to_print";

print_omniscient( {omniscient => $hash_omniscient, output => $gffout} );

# --- final messages ---
end_script();

# ----------------------------------------------------------------------------
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

=item B<-gff> <file>
Input GTF/GFF file.

=item B<-fa> or B<--fasta> <file>

Genome fasta file
The name of the fasta file containing the genome to work with.

=item B<--ct> or B<--table> or B<--codon> <int>

This option allows specifying the codon table to use - It expects an integer (1 by default = standard)

=item B<--size> or B<-s> <int>

Minimum exon size accepted in nucleotide. All exon below this size will be extended to this size. Default value = 15.

=item B<-o>, B<--out> or B<--output> <file>

Output file to create (default GFF3 - see config to modify output format).
If no output file is specified, the output will be written to STDOUT.

=item B<-h> or B<--help>

Display this helpful text.

=back

=head1 SHARED OPTIONS

Shared options are defined in the AGAT configuration file and can be overridden via the command line for this script only.
Common shared options are listed below; for the full list, please refer to the AGAT agat_config.yaml.

=over 8

=item B<--config> <file>

Path to a custom AGAT configuration file.  
By default, AGAT uses `agat_config.yaml` from the working directory if present, otherwise the default file shipped with AGAT
(available locally via `agat config --expose`).

=item B<--cpu>, B<--core>, B<--job> or B<--thread> <int>

Number of parallel processes to use for file input parsing (via forking).

=item B<-v> or B<--verbose> <int>

Verbosity, choice are 0,1,2,3,4. 0 is quiet, 1 is normal, 2,3,4 is more verbose. Default 1.

=back

=head1 FEEDBACK

For questions, suggestions, or general discussions about AGAT, please use the AGAT community forum:
https://github.com/NBISweden/AGAT/discussions

=head1 BUG REPORTING

Bug reports should be submitted through the AGAT GitHub issue tracker:
https://github.com/NBISweden/AGAT/issues

=cut

AUTHOR - Jacques Dainat
