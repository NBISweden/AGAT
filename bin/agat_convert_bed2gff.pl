#!/usr/bin/env perl

use strict;
use warnings;
use Clone;
use Pod::Usage;
use Getopt::Long;
use AGAT::AGAT;

my $header = get_agat_header();
my $config;
my $cpu;
my $outfile = undef;
my $bed = undef;
my $source_tag = "data";
my $primary_tag = "gene";
my $inflating_off = undef;
my $inflate_type = "exon";
my $verbose = undef;
my $help;


if( !GetOptions(  	'c|config=s'     => \$config,
                    'thread|threads|cpu|cpus|core|cores|job|jobs=i' => \$cpu,
					"h|help"         => \$help,
					"bed=s"          => \$bed,
					"source=s"       => \$source_tag,
					"verbose|v!"     => \$verbose,
					"primary_tag=s"  => \$primary_tag,
					"inflate_off!"   => \$inflating_off,
					"inflate_type=s" => \$inflate_type,
					"outfile|output|o|out|gff=s" => \$outfile ) )
{
    pod2usage( { -message => "Failed to parse command line.\n",
                 -verbose => 1,
                 -exitval => 1 } );
}

# Print Help and exit
if ($help) {
	pod2usage( {-message => "$header\n",
	            -verbose => 99,
	            -exitval => 0 } );
}

if ( ! (defined($bed)) ){
    pod2usage( {
           -message => "$header\nAt least 1 parameter is mandatory:\nInput bed file (--bed).\n\n",
           -verbose => 0,
           -exitval => 1 } );
}

# --- Manage config ---
initialize_agat({ config_file_in => $config, input => $bed });
$CONFIG->{cpu} = $cpu if defined($cpu);

## Manage output file
my $gffout = prepare_gffout( $outfile );

# Ask for specific GFF information
if (!$source_tag or !$primary_tag){
	print "Some information needed in a GFF3 file dont exist in a BED file...",
	" as we cannot guess them, please fill the needed information:\n\n";

}
if(! $source_tag ){
	print "What is the data source? The source corresponds to the tool that produced the data. Example: Stringtie,Maker,Augustus,etc. [default: data]\n";
	my $source_tag = <STDIN>;
	chomp $source_tag;
	if ($source_tag eq '') {$source_tag = 'data';}
	if ($source_tag =~ /\s/) {die("Whitespace not allowed in $source_tag\n") }
}

if(! $primary_tag ){
	print "What is the data type? Example: gene,mRNA,CDS,etc.  [default: gene]\n";
	my $primary_tag = <STDIN>;
	chomp $primary_tag;
	if ($primary_tag eq '') {$primary_tag = 'gene';}
	if ($primary_tag =~ /\s/) {die("Whitespace not allowed in $primary_tag\n") }
}

### Read bed input file.
open my $fh, $bed or die "Could not open $bed: $!";

                                #######################
                                #        MAIN         #
                                #######################

my %bedOmniscent;
my $UniqID=0;
my $inflate_cpt=0;
my $inflate_cds_cpt=0;
my $inflate_left_cpt=0;
my $inflate_right_cpt=0;
while( my $line = <$fh>)  {
  chomp $line;

	if ($line =~ /#/){ print "skip commented line: $line" if ($verbose); next; } #skip commented lines

  my @fields = split /\t/, $line;
	if (! skip_line($fields[0])){

		# Check if space delimited format. BEED accept both tabulated and space delimited format
		if ($#fields == 0){
			@fields = split /\s/, $line;
		}

    my $fieldNumber=$#fields+1;
    if($fieldNumber < 3 or $fieldNumber >12){
      print "Problem with that line:\n$line\nA bed file has at least three required fields ! 9 others fields are optional. So, a maximum of 12 fields is allowed !",
      "\n Your line contains $fieldNumber fields. Check the sanity of your file. Bye Bye.\n";exit;
    }

    my $cptField=0;
    $UniqID++;
    foreach my $field (@fields){
      $cptField++;
			$field=~ s/,+$//; #removing trailing coma

      ##########
      #MANDATORY fields
      ###########

      # 1 chrom - The name of the chromosome (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10671).
      if($cptField == 1){
        $bedOmniscent{$UniqID}{'chrom'}=$field;
      }

      # 2 chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
      if($cptField == 2){
        $bedOmniscent{$UniqID}{'chromStart'}=$field;
      }

      # 3 chromEnd - The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99.
      if($cptField == 3){
        $bedOmniscent{$UniqID}{'chromEnd'}=$field;
      }

      ##########
      # OPTIONAL fields
      ##########

      # 4 name - Defines the name of the BED line. This label is displayed to the left of the BED line in the Genome Browser window when the track is open to full display mode or directly to the left of the item in pack mode.
      if($cptField == 4){
        $bedOmniscent{$UniqID}{'name'}=$field;
      }

      # 5 score - A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray).
      if($cptField == 5){
        $bedOmniscent{$UniqID}{'score'}=$field;
      }

      # 6 strand - Defines the strand - either '+' or '-'.
      if($cptField == 6){
        $bedOmniscent{$UniqID}{'strand'}=$field;
      }

      # 7 thickStart - The starting position at which the feature is drawn thickly
			# (for example, the start codon in gene displays). When there is no thick part,
			# thickStart and thickEnd are usually set to the chromStart position.
      if($cptField == 7){
        $bedOmniscent{$UniqID}{'thickStart'}=$field;
      }

      # 8 thickEnd - The ending position at which the feature is drawn thickly
			# (for example, the stop codon in gene displays).
      if($cptField == 8){
        $bedOmniscent{$UniqID}{'thickEnd'}=$field;
      }

      # 9 itemRgb - An RGB value of the form R,G,B (e.g. 255,0,0). If the track line itemRgb attribute is set to "On", this RBG value will determine the display color of the data contained in this BED line. NOTE: It is recommended that a simple color scheme (eight colors or less) be used with this attribute to avoid overwhelming the color resources of the Genome Browser and your Internet browser.
      if($cptField == 9){
        $bedOmniscent{$UniqID}{'itemRgb'}=$field;
      }

      # 10 blockCount - The number of blocks (exons) in the BED line.
      if($cptField == 10){
        $bedOmniscent{$UniqID}{'blockCount'}=$field;
      }

      # 11 blockSizes - A comma-separated list of the block sizes. The number of items in this list should correspond to blockCount.
      if($cptField == 11){
        $bedOmniscent{$UniqID}{'blockSizes'}=$field;
      }

      # 12 blockStarts - A comma-separated list of block starts. All of the blockStart positions should be calculated relative to chromStart. The number of items in this list should correspond to blockCount.
      if($cptField == 12){
        $bedOmniscent{$UniqID}{'blockStarts'}=$field;
      }
		}
  }
}

###
# MANAGE BED OMNISCIEN FOR OUTPUT
foreach my $id ( sort {$a <=> $b} keys %bedOmniscent){
#  foreach my $key (keys %{$bedOmniscent{$id}}){

    my $seq_id=$bedOmniscent{$id}{'chrom'};

    #my $source_tag; #fill at the beginning

    #my $primary_tag; #fill at the beginning

    my $start=$bedOmniscent{$id}{'chromStart'}+1; # shift from 0-base to 1-based

    my $end=$bedOmniscent{$id}{'chromEnd'};

    my $frame=".";

    my $score;
    if( exists_keys (\%bedOmniscent, ($id, 'score') ) ){
      $score=$bedOmniscent{$id}{'score'};
    }

    my $strand;
    if( exists_keys (\%bedOmniscent, ($id, 'strand') ) ){
      $strand=$bedOmniscent{$id}{'strand'};
    }

    my $feature = Bio::SeqFeature::Generic->new(-seq_id => $seq_id,
                                                -source_tag => $source_tag,
                                                -primary_tag => $primary_tag,
                                                -start => $start,
																								-end => $end ,
                                                -frame => $frame ,
                                                -strand =>$strand,
                                                -tag => {'ID' => $id}
                                                ) ;

    if( exists_keys ( \%bedOmniscent, ($id, 'name') ) ){
      $feature->add_tag_value('Name',$bedOmniscent{$id}{'name'});
    }
		if( exists_keys ( \%bedOmniscent, ($id, 'thickStart') ) ){
      $feature->add_tag_value('thickStart',$bedOmniscent{$id}{'thickStart'});
    }
		if( exists_keys ( \%bedOmniscent, ($id ,'thickEnd') ) ){
			$feature->add_tag_value('thickEnd',$bedOmniscent{$id}{'thickEnd'});
		}
		if( exists_keys ( \%bedOmniscent, ($id, 'itemRgb') ) ){
			$feature->add_tag_value('itemRgb',$bedOmniscent{$id}{'itemRgb'});
		}
		if( exists_keys ( \%bedOmniscent, ($id, 'blockCount') ) ){
			$feature->add_tag_value('blockCount',$bedOmniscent{$id}{'blockCount'});
		}
		if( exists_keys ( \%bedOmniscent, ($id, 'blockSizes') ) ){
			$feature->add_tag_value('blockSizes',$bedOmniscent{$id}{'blockSizes'});
		}
		if( exists_keys ( \%bedOmniscent, ($id, 'blockStarts') ) ){
			$feature->add_tag_value('blockStarts',$bedOmniscent{$id}{'blockStarts'});
		}

    $gffout->write_feature($feature);

		if ( exists_keys ( \%bedOmniscent, ($id, 'blockCount') ) and ! $inflating_off){
			print "inflating $inflating_off\n" if ($verbose);
			my $l3_start_line = $bedOmniscent{$id}{'blockStarts'};
			$l3_start_line =~ s/^\s+//; # remove spaces
			my @l3_start_list = split /,/, $l3_start_line;

			my $l3_size_line = $bedOmniscent{$id}{'blockSizes'};
			$l3_size_line =~ s/^\s+//; # remove spaces
			my @l3_size_list = split /,/, $l3_size_line;

			if ($#l3_size_list != $#l3_start_list){warn "Error: Number of elements in blockSizes (11th column) blockStarts (12th column) is different!\n";}

			my $l3_indice=-1;
			my $phase = "." ;
			my @list_l3;
			foreach my $l3_start (sort {$a <=> $b} @l3_start_list){
				$l3_indice++;
				$inflate_cpt++;

				$l3_start = $start+$l3_start; # positions should be calculated relative to chromStart
				my $l3_end = $l3_start+$l3_size_list[$l3_indice]-1;
				my $Inflate_ID=$inflate_type.$inflate_cpt;
				my $feature = Bio::SeqFeature::Generic->new(-seq_id => $seq_id,
																										-source_tag => $source_tag,
																										-primary_tag => $inflate_type,
																										-start => $l3_start,
																										-end => $l3_end ,
																										-frame => "." ,
																										-strand =>$strand,
																										tag => {'ID' => $Inflate_ID, 'Parent' => $id}
																										) ;
				push @list_l3, $feature;
			}
			# set outside in case minus strand because could not have been computed on the fly
			if (lc($inflate_type) eq "cds"){
				$phase = 0;
				@list_l3 = $strand eq "+" ? sort {$a->start <=> $b->start} @list_l3 : sort {$b->start <=> $a->start} @list_l3	;
				# set phases
				# compute the phase. Assuming it always start at 0. No fragmented prediction. Has to be done last
				foreach my $feature (@list_l3){
					$feature->frame($phase);
					# compute the phase. Assuming it always start at 0. No fragmented prediction. Has to be done last
					my $cds_length = $feature->end - $feature->start + 1;
					$phase = (3-(($cds_length-$phase)%3))%3; #second modulo allows to avoid the frame with 3. Instead we have 0.
				}
			}
			# print l3 features
			foreach my $feature (sort {$a->start <=> $b->start}  @list_l3){
					$gffout->write_feature($feature);
			}

			# add CDS if inflating is $inflate_type=exon
			if ( ($inflate_type eq "exon") and  exists_keys( \%bedOmniscent, ($id, 'thickStart') ) and exists_keys( \%bedOmniscent, ($id, 'thickEnd') ) ){
				my $phase = ".";
				my $l3_indice=-1;
				my $write_cds = undef;
				my $left_utr_start = undef;
				my $left_utr_end = undef;
				my $right_utr_start = undef;
				my $right_utr_end = undef;
				my $left_utr_type = ( $strand eq "+") ? "five_prime_UTR" : "three_prime_UTR" ;
				my $right_utr_type = ( $strand eq "+") ? "three_prime_UTR" : "five_prime_UTR" ;
				my @cds;
				my @utrs;
				foreach my $l3_start (sort {$a <=> $b} @l3_start_list){

					$l3_indice++;

					my $thickStart = $bedOmniscent{$id}{'thickStart'};
					my $thickEnd = $bedOmniscent{$id}{'thickEnd'};

					# exon position
					my $l3_end = $l3_start+$l3_size_list[$l3_indice]-1;
					#left utr
					my $left_utr_start = ( $l3_start < $thickStart ) ? $l3_start : undef ;
					$left_utr_end = $l3_end;
					#right utr
					my $right_utr_start = ( $l3_end > $thickEnd ) ? $l3_start : undef ;
					$right_utr_end = $l3_end;
					# set default cds
					my $cds_start = $l3_start;
					my $cds_end = $l3_end;
					# set first CDS
					if ( ( $thickStart <= $l3_end ) and ( $thickStart >= $l3_start ) ){ # CDS overlaps exon
						$write_cds = 1;
						$cds_start = $thickStart+1;
						$left_utr_end = $thickStart;
					}
					elsif( $thickStart < $l3_start ){
						$write_cds = 1;
					}
					# set last CDS
					#print "$thickEnd <= $l3_end and $thickEnd >= $l3_start\n";
					if ( ( $thickEnd <= $l3_end ) and ( $thickEnd >= $l3_start ) ){ # CDS overlaps exon
						$write_cds = 2;
						$cds_end = $thickEnd;
						$right_utr_start = $thickEnd + 1;
					}
					elsif( $thickEnd < $l3_start ){
						$write_cds = undef;
					}

					if ($left_utr_start){
						$inflate_left_cpt++;
						my $Inflate_ID = $left_utr_type.$inflate_left_cpt;
						my $feature = Bio::SeqFeature::Generic->new(-seq_id => $seq_id,
																												-source_tag => $source_tag,
																												-primary_tag => $left_utr_type,
																												-start => $left_utr_start,
																												-end => $left_utr_end ,
																												-frame => $phase ,
																												-strand =>$strand,
																												tag => {'ID' => $Inflate_ID, 'Parent' => $id}
																												) ;
						push @utrs, $feature;
					}
					if ($write_cds){
						$inflate_cds_cpt++;
						my $Inflate_ID = "CDS".$inflate_cds_cpt;
						my $feature = Bio::SeqFeature::Generic->new(-seq_id => $seq_id,
																												-source_tag => $source_tag,
																												-primary_tag => "CDS",
																												-start => $cds_start,
																												-end => $cds_end ,
																												-frame => $phase ,
																												-strand =>$strand,
																												tag => {'ID' => $Inflate_ID, 'Parent' => $id}
																												) ;
						push @cds, $feature;
						$write_cds = undef if ($write_cds == 2); #deactivate to stop before UTR
					}
					if ($right_utr_start){
						$inflate_right_cpt++;
						my $Inflate_ID = $right_utr_type.$inflate_right_cpt;
						my $feature = Bio::SeqFeature::Generic->new(-seq_id => $seq_id,
																												-source_tag => $source_tag,
																												-primary_tag => $right_utr_type,
																												-start => $right_utr_start,
																												-end => $right_utr_end ,
																												-frame => $phase ,
																												-strand =>$strand,
																												tag => {'ID' => $Inflate_ID, 'Parent' => $id}
																												) ;
						push @utrs, $feature;
					}
				}
				$phase = 0;
				@cds = $strand eq "+" ? sort {$a->start <=> $b->start} @cds : sort {$a->start <=> $b->start} @cds;
				# set phases
				foreach my $feature (@cds){
					$feature->frame($phase);
					# compute the phase. Assuming it always start at 0. No fragmented prediction. Has to be done last
					my $cds_length = $feature->end - $feature->start + 1;
					$phase = (3-(($cds_length-$phase)%3))%3; #second modulo allows to avoid the frame with 3. Instead we have 0.
				}
				foreach my $feature (sort {$a->start <=> $b->start}  @cds){
						$gffout->write_feature($feature);
				}
				foreach my $feature (sort {$a->start <=> $b->start} @utrs){
						$gffout->write_feature($feature);
				}
			}
		}
}

close $fh;

# check if the line has to be skipped or not
sub skip_line{
	my ($field0)=@_;

	my $skip = undef;

	if($field0 eq ""){
		$skip=1;
	}
	if($field0 =~ /^track/){
		print "Skip track line, we skip it because we cannot render it properly in a gff file.\n" if ($verbose);
		$skip=1;
	}
	if($field0 =~ /^browser/){
		print "Skip browser line, we skip it because we cannot render it properly in a gff file.\n" if ($verbose);
		$skip=1;
	}
	return $skip;
}

__END__

=head1 NAME

agat_convert_bed2gff.pl

=head1 DESCRIPTION

The script takes a bed file as input, and will translate it in gff format.
The BED format is described here: https://genome.ucsc.edu/FAQ/FAQformat.html#format1
The script converts 0-based, half-open [start-1, end) bed file to
1-based, closed [start, end] General Feature Format v3 (GFF3).

=head1 SYNOPSIS

    agat_convert_bed2gff.pl --bed infile.bed [ -o outfile ]
    agat_convert_bed2gff.pl -h

=head1 OPTIONS

=over 8

=item B<--bed>

Input bed file that will be converted.

=item B<--source>

The source informs about the tool used to produce the data and is stored in 2nd field of a gff file.
Example: Stringtie,Maker,Augustus,etc. [default: data]

=item B<--primary_tag>

The primary_tag corresponds to the data type and is stored in 3rd field of a gff file.
Example: gene,mRNA,CDS,etc.  [default: gene]

=item B<--inflate_off>

By default we inflate the block fields (blockCount, blockSizes, blockStarts) to create subfeatures
of the main feature (primary_tag). The type of subfeature created is based on the
inflate_type parameter. If you do not want this inflating behaviour you can deactivate it
by using the --inflate_off option.

=item B<--inflate_type>

Feature type (3rd column in gff) created when inflate parameter activated [default: exon].

=item B<--verbose>

add verbosity

=item B<-o> , B<--output> , B<--out> , B<--outfile> or B<--gff>

Output GFF file. If no output file is specified, the output will be
written to STDOUT.

=item B<-thread>, B<threads>, B<cpu>, B<cpus>, B<core>, B<cores>, B<job> or B<jobs>

Integer — Number of parallel processes to use for file input parsing (via forking).

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
