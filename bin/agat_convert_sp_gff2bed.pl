#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Bio::Tools::GFF;
use Pod::Usage;
use AGAT::Omniscient;

my $header = get_agat_header();
my $outfile = undef;
my $gff = undef;
my $sub = "exon";
my $help;

if( !GetOptions(
    "help|h" => \$help,
    "gff=s" => \$gff,
		"sub=s" => \$sub,
    "outfile|output|out|o=s" => \$outfile))
{
    pod2usage( { -message => "Failed to parse command line.",
                 -verbose => 1,
                 -exitval => 1 } );
}
# Print Help and exit
if ($help) {
    pod2usage( { -message => "$header",
                 -verbose => 99,
                 -exitval => 0 } );
}

if ( ! (defined($gff)) ){
    pod2usage( {
           -message => "$header\nAt least 1 parameters is mandatory. Input gff file (--gff)\n",
           -verbose => 0,
           -exitval => 1 } );
}

## Manage output file
my $bedout;
if ($outfile) {
open(my $fh, '>', $outfile) or die "Could not open file '$outfile' $!";
  #$bedout= Bio::FeatureIO->new(-fh => $fh, -format => 'bed' );
  $bedout=$fh;
}
else{
  #$bedout = Bio::FeatureIO->new(-fh => \*STDOUT,  -format => 'bed');
  $bedout=\*STDOUT ;
}

### Parse GTF input file
my ($hash_omniscient, $hash_mRNAGeneLink) = slurp_gff3_file_JD({ input => $gff
                                                              });
# END parsing


# NOT USED BECAUSE in Bio::FeatureIO::bed:
#  my $block_count = '';  #not implemented, used for sub features
#  my $block_sizes = '';  #not implemented, used for sub features
#  my $block_starts = ''; #not implemented, used for sub feature
#   #################
#   # == LEVEL 1 == #
#   #################
#   foreach my $tag_l1 (keys %{$hash_omniscient->{'level1'}}){ # tag_l1 = gene or repeat etc...
#     foreach my $id_l1 (keys %{$hash_omniscient->{'level1'}{$tag_l1}}){

#       #my feature is a Bio::SeqFeature::Generic
#       my $feature_l1=$hash_omniscient->{'level1'}{$tag_l1}{$id_l1};

#       #create a new  Bio::SeqFeature::Annotated object;
#       my $newObj = Bio::SeqFeature::Annotated->new();
#       #initialize this object with the contents of another feature
#       $newObj->from_feature($feature_l1);
#       #print the new object
#       print $bedout->write_feature($newObj);
#   }
# }

# sort by seq id
my $hash_sortBySeq = gather_and_sort_l1_by_seq_id($hash_omniscient);

#################
# == LEVEL 1 == #
#################
foreach my $seqid (sort { (($a =~ /(\d+)$/)[0] || 0) <=> (($b =~ /(\d+)$/)[0] || 0) } keys %{$hash_sortBySeq}){ # loop over all the feature level1


	foreach my $tag_l1 (sort {$a cmp $b} keys %{$hash_omniscient->{'level1'}}){
		foreach my $feature_l1 ( @{$hash_sortBySeq->{$seqid}{$tag_l1}} ){
			my $id_l1 = lc($feature_l1->_tag_value('ID'));


	    #################
	    # == LEVEL 2 == #
	    #################
			my $field1_chrom = undef;
			my $field2_chromStart = undef;
			my $field3_chromEnd = undef;
			my $field4_name = undef;
			my $field5_score = undef;
			my $field6_strand = undef;
			my $field7_thickStart = undef;
			my $field8_thickEnd = undef;
			my $field9_itemRgb = "255,0,0";
			my $field10_blockCount = undef;
			my $field11_blockSizes = undef;
			my $field12_blockStarts = undef;

	    foreach my $tag_l2 ( sort {$a cmp $b} keys %{$hash_omniscient->{'level2'}}){ # tag_l2 = mrna or mirna or ncrna or trna etc...
				if( exists_keys( $hash_omniscient, ('level2', $tag_l2, $id_l1 ) ) ) {
	        foreach my $feature_l2 ( sort {$a->start <=> $b->start} @{$hash_omniscient->{'level2'}{$tag_l2}{$id_l1}}) {

						$field1_chrom = $feature_l2->seq_id;
						$field2_chromStart = $feature_l2->start;
						$field3_chromEnd = $feature_l2->end;
						$field4_name = $feature_l2->_tag_value('ID');
						$field5_score = $feature_l2->score ;
						if(!$field5_score or $field5_score < 0){$field5_score = "0";}
						$field6_strand = $feature_l2->strand;
							$field6_strand = "+" if ( $field6_strand eq "1");
							$field6_strand = "-" if ( $field6_strand eq "-1");
						$field7_thickStart = $feature_l2->start;
						$field8_thickEnd = $feature_l2->end;

	          #################
	          # == LEVEL 3 == #
	          #################
	          my $level2_ID = lc( $feature_l2->_tag_value('ID') );

						if( exists_keys( $hash_omniscient, ('level3', lc($sub), $level2_ID ) ) ) {

							foreach my $feature_l3 ( sort { $a->start <=> $b->start } @{$hash_omniscient->{'level3'}{lc($sub)}{$level2_ID}}) {
									$field10_blockCount++;

									my $size_l3 = $feature_l3->end - $feature_l3->start; #No +1 as in GFF feature size because 0-based format in bed
									$field11_blockSizes .= $size_l3.",";

									my $start_l3 = $feature_l3->start - 1; #No +1 as in GFF feature size because 0-based format in bed
									$field12_blockStarts .= $start_l3.",";
	            }
	          }
	        }
	      }
	    }

			if ($field11_blockSizes){
				$field11_blockSizes=~ s/,+$//; #removing trailing coma
			}
			if ($field12_blockStarts){
				$field12_blockStarts=~ s/,+$//; #removing trailing coma
			}

			# skip topfeatures e.g chromosome, location, etc. because do not have level2 features.
			if($field1_chrom){
				my $line = $field1_chrom."\t".$field2_chromStart."\t".$field3_chromEnd."\t".$field4_name."\t".$field5_score."\t".$field6_strand."\t".$field7_thickStart."\t".$field8_thickEnd."\t".$field9_itemRgb;
				if($field10_blockCount){
					$line .= "\t".$field10_blockCount."\t".$field11_blockSizes."\t".$field12_blockStarts;
				}
				$line .= "\n";

				print $bedout $line;
			}
	  }
	}
}

__END__


=head1 NAME

agat_convert_sp_gff2bed.pl

=head1 DESCRIPTION

The script aims to convert GTF/GXF file into bed file.
It will convert level2 feature from gff (mRNA, transcripts) into bed feature.
If level2 subfeatures selected (defaut: exon) exist, they will be reported in the
block fields (9-12th colum in bed).

Definintion of the bed format:
# 1 chrom - The name of the chromosome (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10671).
# 2 chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
# 3 chromEnd - The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99.
########## OPTIONAL fields ##########
# 4 name - Defines the name of the BED line. This label is displayed to the left of the BED line in the Genome Browser window when the track is open to full display mode or directly to the left of the item in pack mode.
# 5 score - A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray).
# 6 strand - Defines the strand - either '+' or '-'.
# 7 thickStart - The starting position at which the feature is drawn thickly
# 8 thickEnd - The ending position at which the feature is drawn thickly
# 9 itemRgb - An RGB value of the form R,G,B (e.g. 255,0,0). If the track line itemRgb attribute is set to "On", this RBG value will determine the display color of the data contained in this BED line. NOTE: It is recommended that a simple color scheme (eight colors or less) be used with this attribute to avoid overwhelming the color resources of the Genome Browser and your Internet browser.
# 10 blockCount - The number of blocks (exons) in the BED line.
# 11 blockSizes - A comma-separated list of the block sizes. The number of items in this list should correspond to blockCount.
# 12 blockStarts

=head1 SYNOPSIS

    agat_convert_sp_gff2bed.pl --gff file.gff  [ -o outfile ]
    agat_convert_sp_gff2bed.pl --help

=head1 OPTIONS

=over 8

=item B<--gff>

Input GFF3 file that will be read

=item B<--sub>

Define the subfeature (level3, e.g exon,cds,utr,etc...) to report as blocks in the bed output.
Defaut: exon.

=item B<--outfile>, B<--out>, B<--output>, or B<-o>

File where will be written the result. If no output file is specified, the output will be written to STDOUT.

=item B<-h> or B<--help>

Display this helpful text.

=back

=cut
