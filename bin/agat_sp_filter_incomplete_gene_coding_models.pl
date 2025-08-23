#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use File::Basename;
use List::MoreUtils qw(uniq);
use Bio::DB::Fasta;
use Bio::SeqIO;
use AGAT::AGAT;

my $header = get_agat_header();
my @copyARGV  = @ARGV;

my ( $opt, $usage, $config ) = AGAT::AGAT::describe_script_options(
    $header,
    [ 'gff=s',              'Input reference gff file',   { required => 1 } ],
    [ 'fasta|fa|f=s',       'Input reference fasta file', { required => 1 } ],
    [ 'table|codon|ct=i',   'Codon translation table',    {
        default   => 1,
        callbacks => { positive => sub { shift() > 0 or die 'Codon translation table must be positive' } },
    } ],
    [ 'add_flag|af!',              'Add incomplete attribute flag' ],
    [ 'skip_start_check|sstartc!', 'Skip start codon check' ],
    [ 'skip_stop_check|sstopc!',   'Skip stop codon check' ],
);

my $gff              = $opt->gff;
my $file_fasta       = $opt->fasta;
my $codonTableId     = $opt->table;
my $add_flag         = $opt->add_flag;
my $skip_start_check = $opt->skip_start_check;
my $skip_stop_check  = $opt->skip_stop_check;
my $opt_output       = $config->{output};
my $opt_verbose      = $config->{verbose};

my $log;
if ( my $log_name = $config->{log_path} ) {
    open( $log, '>', $log_name )
      or die "Can not open $log_name for printing: $!";
    dual_print( $log, $header, 0 );
}

# --- Check codon table ---
$codonTableId = get_proper_codon_table($codonTableId, $log, $opt_verbose);

my $codonTable = Bio::Tools::CodonTable->new( -id => $codonTableId);

######################
# Manage output file #
my $gffout_file;
my $gffout_incomplete_file;
if ($opt_output) {
  my ($filename,$path,$ext) = fileparse($opt_output,qr/\.[^.]*/);

  $gffout_file = $path.$filename.$ext;
  $gffout_incomplete_file = $path.$filename."_incomplete".$ext;
}

my $gffout = prepare_gffout($config, $gffout_file);
my $gffout_incomplete = prepare_gffout($config, $gffout_incomplete_file);

                #####################
                #     MAIN          #
                #####################

######################
### Parse GFF input #
my ($hash_omniscient, $hash_mRNAGeneLink) = slurp_gff3_file_JD({ input => $gff,
                                                                 config => $config
                                                              });
dual_print( $log, "GFF3 file parsed\n", $opt_verbose );


####################
# index the genome #
my $nbFastaSeq=0;
my $db = Bio::DB::Fasta->new($file_fasta);
my @ids      = $db->get_all_primary_ids;
my %allIDs; # save ID in lower case to avoid cast problems
foreach my $id (@ids ){$allIDs{lc($id)}=$id;}
dual_print( $log, "Fasta file parsed\n", $opt_verbose );
####################

#counters
my %mrnaCounter=( 0 => 0, 1 => 0, 2 => 0, 3 => 0);
my $geneCounter=0;
my %omniscient_incomplete; initialize_omni_from(\%omniscient_incomplete, $hash_omniscient);
my @incomplete_mRNA;


foreach my $primary_tag_key_level1 (keys %{$hash_omniscient->{'level1'}}){ # primary_tag_key_level1 = gene or repeat etc...
  foreach my $gene_id (keys %{$hash_omniscient->{'level1'}{$primary_tag_key_level1}}){
    my $gene_feature = $hash_omniscient->{'level1'}{$primary_tag_key_level1}{$gene_id};
    my $strand = $gene_feature->strand();
    dual_print( $log, "gene_id = $gene_id\n", $opt_verbose );

    my @level1_list=();
    my @level2_list=();
    my @level3_list=();

    my $ncGene=1;
    foreach my $primary_tag_key_level2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
      if ( exists_keys( $hash_omniscient, ('level2', $primary_tag_key_level2, $gene_id) ) ){

        my $geneInc=undef;
        foreach my $level2_feature ( @{$hash_omniscient->{'level2'}{$primary_tag_key_level2}{$gene_id}}) {
          my $start_missing=undef;
          my $stop_missing=undef;

          # get level2 id
          my $level2_ID = lc($level2_feature->_tag_value('ID'));

          if ( exists_keys( $hash_omniscient, ('level3', 'cds', $level2_ID) ) ){
            $ncGene=undef;
            $mrnaCounter{'0'}++;

            my $seqobj = extract_cds(\@{$hash_omniscient->{'level3'}{'cds'}{$level2_ID}}, $db);
            my $length_CDS = length($seqobj->seq);
            #------------- check CDS length -------------
            if ($length_CDS >= 3){
              #------------- check start -------------
              if (! $skip_start_check){
                my $start_codon = $seqobj->subseq(1,3);
                if(! $codonTable->is_start_codon( $start_codon )){
                  dual_print( $log, "start= $start_codon  is not a valid start codon\n", $opt_verbose );
                  $start_missing="true";
                  if($add_flag){
                    create_or_replace_tag($level2_feature, 'incomplete', '1');
                  }
                }
              }
              #------------- check stop --------------
              if (! $skip_stop_check){
                my $seqlength  = length($seqobj->seq());
                my $stop_codon = $seqobj->subseq($seqlength - 2, $seqlength) ;

                if(! $codonTable->is_ter_codon( $stop_codon )){
                  dual_print( $log, "stop= $stop_codon is not a valid stop codon\n", $opt_verbose );
                  $stop_missing="true";
                  if($add_flag){
                    if($start_missing){
                      create_or_replace_tag($level2_feature, 'incomplete', '3');
                    }
                    else{
                      create_or_replace_tag($level2_feature, 'incomplete', '2');
                    }
                  }
                }
              }
            }
            else{ #short CDS
              dual_print( $log, "CDS too short ($length_CDS nt) we skip it\n", $opt_verbose );
            }
          }
          else{ #No CDS
            dual_print( $log, "Not a coding rna (no CDS) we skip it\n", $opt_verbose );
          }

          if($start_missing or $stop_missing){
            #Keep track counter
            if ($start_missing and $stop_missing) {
              $mrnaCounter{'3'}++;
            }
            elsif($start_missing){
              $mrnaCounter{'1'}++;
            }
            else{
              $mrnaCounter{'2'}++;
            }
            $geneInc="true";

            if(! $add_flag){
              push(@incomplete_mRNA, $level2_ID); # will be removed at the end
              push(@level2_list, $level2_feature); # will be appended to omniscient_incomplete
              foreach my $primary_tag_l3 (keys %{$hash_omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...
                if ( exists ($hash_omniscient->{'level3'}{$primary_tag_l3}{$level2_ID} ) ){
                  push(@level3_list, @{$hash_omniscient->{'level3'}{$primary_tag_l3}{$level2_ID}})
                }
              }
            }
          }
        }
        if($geneInc){
          $geneCounter++;
          #Save the mRNA and parent and child features
          if(! $add_flag){
            @level1_list=($gene_feature);
            append_omniscient(\%omniscient_incomplete, \@level1_list, \@level2_list, \@level3_list);
          }
        }
      }
    }
    #after checking all mRNA of a gene
    if($ncGene){
      dual_print( $log, "This is a non coding gene (no cds to any of its RNAs)", $opt_verbose );
    }
  }
}


#END
my $string_to_print="usage: $0 @copyARGV\n";
$string_to_print .="Results:\n";

if ($geneCounter) {
  $string_to_print .="We checked ".$mrnaCounter{0}." mRNAs.\n";
  $string_to_print .="There are ".$mrnaCounter{3}." mRNAs without start and stop codons.\n";
  $string_to_print .="There are ".$mrnaCounter{2}." mRNAs without stop codons.\n";
  $string_to_print .="There are ".$mrnaCounter{1}." mRNAs without start codons.\n";
  $string_to_print .="Number of gene affected: $geneCounter\n";
}
else{
  $string_to_print .="No gene with incomplete mRNA!\n";
}
dual_print( $log, $string_to_print, $opt_verbose );

if(! $add_flag){
  #clean for printing
  if (@incomplete_mRNA){
    check_all_level2_locations( { omniscient => \%omniscient_incomplete } ); # review all the feature L2 to adjust their start and stop according to the extrem start and stop from L3 sub features.
    check_all_level1_locations( { omniscient => \%omniscient_incomplete } );

    remove_omniscient_elements_from_level2_ID_list($hash_omniscient, \@incomplete_mRNA);
		check_all_level2_locations( { omniscient => $hash_omniscient } ); # review all the feature L2 to adjust their start and stop according to the extrem start and stop from L3 sub features.
    check_all_level1_locations( { omniscient => $hash_omniscient } ); # Check the start and end of level1 feature based on all features level2.
  }
}

dual_print( $log, "Now printing complete models\n", $opt_verbose );
print_omniscient( {omniscient => $hash_omniscient, output => $gffout} );

if(@incomplete_mRNA){
  dual_print( $log, "Now printing incomplete models\n", $opt_verbose );
  print_omniscient( {omniscient => \%omniscient_incomplete, output => $gffout_incomplete} );
}

dual_print( $log, "Bye Bye.\n", $opt_verbose );
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

sub extract_cds{
  my($feature_list, $db)=@_;

  my @sortedList = sort {$a->start <=> $b->start} @$feature_list;
  my $sequence="";
  foreach my $feature ( @sortedList ){
    $sequence .= get_sequence($db, $feature->seq_id, $feature->start, $feature->end);
  }

  #create sequence object
  my $seq  = Bio::Seq->new( '-format' => 'fasta' , -seq => $sequence);

  #check if need to be reverse complement
  if($sortedList[0]->strand eq "-1" or $sortedList[0]->strand eq "-"){
    $seq=$seq->revcom;
  }
  return $seq;
}

sub  get_sequence{
  my  ($db, $seq_id, $start, $end) = @_;

  my $sequence="";
  my $seq_id_correct = undef;
  if( exists $allIDs{lc($seq_id)}){

    $seq_id_correct = $allIDs{lc($seq_id)};

    $sequence = $db->subseq($seq_id_correct, $start, $end);

    if($sequence eq ""){
      my $msg = "Problem ! no sequence extracted for - $seq_id !\n";
      dual_print( $log, $msg, 0 );
      warn $msg if $opt_verbose;
      exit;
    }
    if( length($sequence) != ($end-$start+1) ){
      my $wholeSeq = $db->subseq($seq_id_correct);
      $wholeSeq = length($wholeSeq);
      my $msg = "Problem ! The size of the sequence extracted " . length($sequence) .
                " is different than the specified span: " . ($end-$start+1) . ".\n" .
                "That often occurs when the fasta file does not correspond to the annotation file." .
                " Or the index file comes from another fasta file which had the same name and haven't been removed.\n" .
                "As last possibility your gff contains location errors (Already encountered for a Maker annotation)\n" .
                "Supplement information: seq_id=$seq_id ; seq_id_correct=$seq_id_correct ; start=$start ; end=$end ; sequence length: $wholeSeq )\n";
      dual_print( $log, $msg, 0 );
      warn $msg if $opt_verbose;
    }
  }
  else{
    my $msg = "Problem ! ID $seq_id not found !\n";
    dual_print( $log, $msg, 0 );
    warn $msg if $opt_verbose;
  }

  return $sequence;
}

__END__

=head1 NAME

agat_sp_filter_incomplete_gene_coding_models.pl

=head1 DESCRIPTION

The script aims to remove incomplete gene models. An incomplete gene coding model
is a gene coding with start and/or stop codon missing in its cds.
You can modify the behavior using the skip_start_check or skip_stop_check options.

=head1 SYNOPSIS

    agat_sp_filter_incomplete_gene_coding_models.pl --gff infile.gff --fasta genome.fa [ -o outfile ]
    agat_sp_filter_incomplete_gene_coding_models.pl --help

=head1 OPTIONS

=over 8

=item B<-gff>

Input GTF/GFF file.

=item B<-fa> or B<--fasta>

Genome fasta file.
The name of the fasta file containing the genome to work with.

=item B<--ct> or B<--table> or B<--codon>

This option allows specifying the codon table to use.
It expects an integer [default 1]

=item B<--af> or B<--add_flag>

Instead of filter the result into two output files, write only one and add the flag <incomplete> in the gff.(tag = inclomplete, value = 1, 2, 3.  1=start missing; 2=stop missing; 3=both)

=item B<--skip_start_check> or B<--sstartc>

Gene model must have a start codon. Activated by default.

=item B<--skip_stop_check> or B<--sstopc>

Gene model must have a stop codon. Activated by default.

=item B<-o> , B<--output> , B<--out> or B<--outfile>

Output GFF file.  If no output file is specified, the output will be
written to STDOUT.

=item B<-v>

Verbose option, make it easier to follow what is going on for debugging purpose.

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
