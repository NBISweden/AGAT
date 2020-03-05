#!/usr/bin/env perl

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use Bio::Tools::GFF;
use AGAT::Omniscient;

my $header = get_agat_header();
my $outfile = undef;
my $gff = undef;
my $nf =undef;
my $help;

if( !GetOptions(
    "help" => \$help,
    "gff|in=s" => \$gff,
    "nf" => \$nf,
    "outfile|output|o|out|gtf=s" => \$outfile))
{
    pod2usage( { -message => "Failed to parse command line.",
                 -verbose => 1,
                 -exitval => 1 } );
}

# Print Help and exit
if ($help) {
    pod2usage( { -message => "$header\n",
                 -verbose => 99,
                 -exitval => 0 } );
}

if ( ! (defined($gff)) ){
    pod2usage( {
           -message => "$header\nAt least 1 parameters is mandatory:\nInput gff file (--gff).\n\n",
           -verbose => 0,
           -exitval => 1 } );
}

## Manage output file
my $gtf_out;
if ($outfile) {
open(my $fh, '>', $outfile) or die "Could not open file '$outfile' $!";
  $gtf_out = Bio::Tools::GFF->new(-fh => $fh, -gff_version => 2.5);
}
else{
  $gtf_out = Bio::Tools::GFF->new(-fh => \*STDOUT, -gff_version => 2.5);
}

######################
### Parse GFF input #
### Read gff input file.
my ($hash_omniscient, $hash_mRNAGeneLink) = slurp_gff3_file_JD({ input => $gff
                                                              });

# rebuild gene_id and transcript_id feature;

  my $gene_id=undef;
  #################
  # == LEVEL 1 == #
  #################
  foreach my $primary_tag_key_level1 (keys %{$hash_omniscient->{'level1'}}){ # primary_tag_key_level1 = gene or repeat etc...
    foreach my $id_tag_key_level1 (keys %{$hash_omniscient->{'level1'}{$primary_tag_key_level1}}){
      my $feature_level1=$hash_omniscient->{'level1'}{$primary_tag_key_level1}{$id_tag_key_level1};

      # Gene ID level1
      my $gene_id_att=undef;
      if($feature_level1->has_tag('gene_id')){
        $gene_id_att=$feature_level1->_tag_value('gene_id');
      }

      my $transcript_id=undef;
      my $level3_gene_id=undef;
      #################
      # == LEVEL 2 == #
      #################
      foreach my $primary_tag_key_level2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...

        if ( exists ($hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1} ) ){
          foreach my $feature_level2 ( @{$hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1}}) {



            # Gene ID level2
            my $gene_id_mrna_att=undef;
            if($feature_level2->has_tag('gene_id')){
              $gene_id_mrna_att=$feature_level2->_tag_value('gene_id');
            }

            my $transcript_id_mrna_att=undef;
            if($feature_level2->has_tag('transcript_id')){
              $transcript_id_mrna_att=$feature_level2->_tag_value('transcript_id');
            }

            # get gff3 feature (ID)
            my $level2_ID = lc($feature_level2->_tag_value('ID'));

            my $level3_transcript_id=undef;
            #################
            # == LEVEL 3 == #
            #################

            ############
            # Go through one time to check if  gene_id and transcript_id are present and save them
            foreach my $primary_tag_key_level3 (keys %{$hash_omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...
              if ( exists ($hash_omniscient->{'level3'}{$primary_tag_key_level3}{$level2_ID} ) ){
                foreach my $feature_level3 ( @{$hash_omniscient->{'level3'}{$primary_tag_key_level3}{$level2_ID}}) {

                  #Get level3 gene_id
                  if(! $level3_gene_id){
                    if($feature_level3->has_tag('gene_id')){
                      $level3_gene_id=$feature_level3->_tag_value('gene_id');
                    }
                  }

                  #Get level3 transcript_id
                  if(! $level3_transcript_id){
                    if($feature_level3->has_tag('transcript_id')){
                      $level3_transcript_id=$feature_level3->_tag_value('transcript_id');
                    }
                  }
                  if($level3_gene_id and $level3_transcript_id){last;}
                }
              }
              if($level3_gene_id and $level3_transcript_id){last;}
            }

            #################
            # CHOOSE the gene_id. We take the first from level1 to level3.
            if($gene_id_att){
              $gene_id=$gene_id_att;
            }
            elsif($gene_id_mrna_att){
              $gene_id=$gene_id_mrna_att
            }
            elsif($level3_gene_id){
              $gene_id=$level3_gene_id;
            }
            else{ # We didn't find any gene_id we will the ID of level1 as gene_id.
              $gene_id=$feature_level1->_tag_value('ID');
            }

            #################
            # CHOOSE the transcript_id. We take the first from level2 to level3.
            if($transcript_id_mrna_att){
              $transcript_id=$transcript_id_mrna_att;
            }
            elsif($level3_transcript_id){
              $transcript_id=$level3_transcript_id;
            }
            else{ # We didn't find any gene_id we will the ID of level2 as transcript_id.
              $transcript_id=$feature_level2->_tag_value('ID');
            }

            ##############
            # Second pass of level3 features
            # Add gene_id and transcript_id to level3 feature that don't have this information
            foreach my $primary_tag_key_level3 (keys %{$hash_omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...
              if ( exists ($hash_omniscient->{'level3'}{$primary_tag_key_level3}{$level2_ID} ) ){
                foreach my $feature_level3 ( @{$hash_omniscient->{'level3'}{$primary_tag_key_level3}{$level2_ID}}) {

                  #Check add gene_id
                  if(! $feature_level3->has_tag('gene_id')) {
                    $feature_level3->add_tag_value('gene_id', $gene_id);
                  }
                  elsif($feature_level3->_tag_value('gene_id') ne $gene_id) { #gene_id different, we replace it.
                    warn("We replace the transcript_id ".$feature_level3->_tag_value('gene_id')." by ".$gene_id.". Is it normal ?\n");exit;
                    $feature_level3->add_tag_value('gene_id', $gene_id);
                  }
                  #Check add transcript_id
                  if(! $feature_level3->has_tag('transcript_id')){
                    $feature_level3->add_tag_value('transcript_id', $transcript_id);
                  }
                  elsif($feature_level3->_tag_value('transcript_id') ne $transcript_id){ #transcript_id different, we replace it.
                    warn("We replace the transcript_id ".$feature_level3->_tag_value('transcript_id')." by ".$transcript_id.". Is it normal ?\n");exit;
                    $feature_level3->add_tag_value('transcript_id', $transcript_id);
                  }
                }
              }
            }

            ## add level2 missing information gene_id
            if(! $feature_level2->has_tag('gene_id')) {
               $feature_level2->add_tag_value('gene_id', $gene_id);
            }
            elsif($feature_level2->_tag_value('gene_id') ne $gene_id) { #gene_id different, we replace it.
              warn("We replace the transcript_id ".$feature_level2->_tag_value('gene_id')." by ".$gene_id.". Is it normal ?\n");exit;
              $feature_level2->add_tag_value('gene_id', $gene_id);
            }
            # add level2 missing information transcript_id
            if(! $feature_level2->has_tag('transcript_id')){
              $feature_level2->add_tag_value('transcript_id', $transcript_id);
            }
            elsif($feature_level2->_tag_value('transcript_id') ne $transcript_id){ #gene_id transcript_id, we replace it.
              warn("We replace the transcript_id ".$feature_level2->_tag_value('transcript_id')." by ".$transcript_id.". Is it normal ?\n");exit;
              $feature_level2->add_tag_value('transcript_id', $transcript_id);
            }
          }
        }
      }
      ## add level1 missing information gene_id
      if(! $feature_level1->has_tag('gene_id')) {
        $feature_level1->add_tag_value('gene_id', $gene_id);
      }
      elsif($feature_level1->_tag_value('gene_id') ne $gene_id) { #gene_id different, we replace it.
        warn("We replace the transcript_id ".$feature_level1->_tag_value('gene_id')." by ".$gene_id.". Is it normal ?\n");exit;
        $feature_level1->add_tag_value('gene_id', $gene_id);
      }
    }
  }
  # print results
  print_omniscient($hash_omniscient, $gtf_out);


if($nf){
  if($outfile){
    `cp $outfile tmp`;
    `sed 's/  / /g' tmp > tmp2`;
    `sed 's/ ;/;/g' tmp2 > $outfile`;
    `rm tmp tmp2`;
  }
  else{print "!! option nf can be used only to write result in a file. Doesn't work with STDOUT\n";}
}
else{
  print "\nKeep in mind that the current format of attibutes/values of the 9th column  is like that:\n".
  "attribute1 value1 ; attribute2 value2\nSome tools (i.e Kraken) struggle with the space between <value> and <;>. If you want to remove it relaunch the script using the <nf> option.\n";
}

print "Bye Bye\n";

__END__

=head1 NAME

agat_convert_sp_gff2gtf.pl

=head1 DESCRIPTION

The script aims to convert any GTF/GFF file into a proper GTF file.
Full information about the format can be found here: https://github.com/NBISweden/GAAS/blob/master/annotation/knowledge/gxf.md
The last descrption of the fomat specify only 9 acctepeted feature  type (3rd colum):
gene, transcript, exon, CDS, Selenocysteine, start_codon, stop_codon, three_prime_utr and five_prime_utr
Nevertheless if your file contains other type of features they will not be removed,
as long as the parser can deal with them.

To be fully GTF compliant all feature need to have a gene_id and a transcript_id attribute.
The gene_id	is unique identifier for the genomic source of the transcript, which is
used to group transcripts into genes.
The transcript_id	is a unique identifier for the predicted transcript,
which is used to group features into transcripts.


Keep in mind that some bioperl versions forget to add the header (##gff-version 2) in the output.
Check the output to add it if missing, it will avoid you troubles during your downstream analyses.

=head1 SYNOPSIS

    agat_convert_sp_gff2gtf.pl --gff infile.gtf [ -o outfile ]
		agat_convert_sp_gff2gtf -h

=head1 OPTIONS

=over 8

=item B<--gff> or B<--in>

Input GFF file that will be read

=item B<--att> or B<-a>

With this option, attributes "gene_id" and "transcript_id" will be created when they are missing.

=item B<--nf>

With this option, attibutes/values of the 9th column are written
"attribute1 value1; attribute2 value2" instead of "attribute1 value1 ; attribute2 value2".
(The difference is the space before the semilicon)

=item B<-o> , B<--output> , B<--out> , B<--outfile> or B<--gtf>

Output GTF file. If no output file is specified, the output will be
written to STDOUT.

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
