#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use List::MoreUtils qw(uniq);
use Bio::DB::Fasta;
use Bio::SeqIO;
use AGAT::AGAT;

my $header = get_agat_header();
my $config;
my $cpu;
my $outfile = undef;
my $gff = undef;
my $file_fasta=undef;
my $codonTableId=1;
my $skip_start_check=undef;
my $skip_stop_check=undef;
my $add_flag=undef;
my $verbose = undef;
my $opt_help= 0;

my @copyARGV=@ARGV;
if ( !GetOptions(
    'c|config=s'                => \$config,
                    'thread|threads|cpu|cpus|core|cores|job|jobs=i' => \$cpu,
    "h|help"                    => \$opt_help,
    "gff=s"                     => \$gff,
    "fasta|fa|f=s"              => \$file_fasta,
    "table|codon|ct=i"          => \$codonTableId,
    "add_flag|af!"              => \$add_flag,
    "skip_start_check|sstartc!" => \$skip_start_check,
    "skip_stop_check|sstopc!"   => \$skip_stop_check,
    "v!"                        => \$verbose,
    "output|outfile|out|o=s"    => \$outfile))

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
initialize_agat({ config_file_in => $config, input => $gff });
$CONFIG->{cpu} = $cpu if defined($cpu);

# --- Check codon table ---
$codonTableId = get_proper_codon_table($codonTableId);

my $codonTable = Bio::Tools::CodonTable->new( -id => $codonTableId);

######################
# Manage output file #
my $gffout_file;
my $gffout_incomplete_file;
if ($outfile) {
  my ($filename,$path,$ext) = fileparse($outfile,qr/\.[^.]*/);

	$gffout_file = $path.$filename.$ext;
  $gffout_incomplete_file = $path.$filename."_incomplete".$ext;
}

my $gffout = prepare_gffout( $gffout_file );
my $gffout_incomplete = prepare_gffout( $gffout_incomplete_file );

                #####################
                #     MAIN          #
                #####################

######################
### Parse GFF input #
my ($hash_omniscient) = slurp_gff3_file_JD({ input => $gff });

####################
# index the genome #
my $nbFastaSeq=0;
my $db = Bio::DB::Fasta->new($file_fasta);
my @ids      = $db->get_all_primary_ids;
my %allIDs; # save ID in lower case to avoid cast problems
foreach my $id (@ids ){$allIDs{lc($id)}=$id;}
print ("Fasta file parsed\n");
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
    print "gene_id = $gene_id\n" if $verbose;

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
                  print "start= $start_codon  is not a valid start codon\n" if ($verbose);
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
                  print "stop= $stop_codon is not a valid stop codon\n" if ($verbose);
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
              print "CDS too short ($length_CDS nt) we skip it\n" if ($verbose);
            }
          }
          else{ #No CDS
            print "Not a coding rna (no CDS) we skip it\n" if ($verbose);
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
      print "This is a non coding gene (no cds to any of its RNAs)" if ($verbose);
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
print $string_to_print;

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

print "Now printing complete models\n";
print_omniscient( {omniscient => $hash_omniscient, output => $gffout} );

if(@incomplete_mRNA){
  print "Now printing incomplete models\n";
  print_omniscient( {omniscient => \%omniscient_incomplete, output => $gffout_incomplete} );
}

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
      warn "Problem ! no sequence extracted for - $seq_id !\n";  exit;
    }
    if( length($sequence) != ($end-$start+1) ){
      my $wholeSeq = $db->subseq($seq_id_correct);
      $wholeSeq = length($wholeSeq);
      warn "Problem ! The size of the sequence extracted ".length($sequence).
			" is different than the specified span: ".($end-$start+1).".\n".
			"That often occurs when the fasta file does not correspond to the annotation file.".
			" Or the index file comes from another fasta file which had the same name and haven't been removed.\n".
      "As last possibility your gff contains location errors (Already encountered for a Maker annotation)\n".
			"Supplement information: seq_id=$seq_id ; seq_id_correct=$seq_id_correct ; start=$start ; end=$end ; sequence length: $wholeSeq )\n";
    }
  }
  else{
    warn "Problem ! ID $seq_id not found !\n";
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

=item B<--ad> or B<--add_flag>

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
