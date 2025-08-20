#!/usr/bin/env perl

use strict;
use warnings;
use Clone 'clone';
use Pod::Usage;
use Getopt::Long;
use Bio::Tools::Genscan;
use AGAT::AGAT;

my $header = get_agat_header();
my $config;
my $outfile = undef;
my $genscan = undef;
my $seq_id = "unknown";
my @list_exon_types = ("Init", "Intr", "Term", "Sngl");
my %hash_raw;
my %hash_uniqID;
my $verbose = undef;
my $help;

my $common = parse_common_options() || {};
$config   = $common->{config};
$outfile  = $common->{output};
$verbose  = $common->{verbose};
$help     = $common->{help};


if( !GetOptions(    'c|config=s'                => \$config,
					"h|help"                    => \$help,
					"g|genscan=s"               => \$genscan,
					"seqid=s"                   => \$seq_id,
					"verbose|v!"                => \$verbose,
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

if ( ! (defined($genscan)) ){
    pod2usage( {
           -message => "$header\nAt least 1 parameter is mandatory:\nInput genscan file (--genscan).\n\n",
           -verbose => 0,
           -exitval => 1 } );
}

# --- Manage config ---
$config = get_agat_config({config_file_in => $config});

my $log;
my $log_name = get_log_path($common, $config);
open($log, '>', $log_name) or die "Can not open $log_name for printing: $!";
dual_print($log, $header, 0);

## Manage output file
my $gffout = prepare_gffout($config, $outfile);

## MAIN ##############################################################

read_genscan($genscan);
#use Data::Dumper; print Dumper (\%hash_raw);
convert_genscan();

## SUBROUTINES #######################################################

sub read_genscan {
	open(INPUT, "<", $genscan) or die ("$!\n");
	# Open Mfannot file for reading
	while (<INPUT>) {
			chomp;

			$_=~ s/^\s+//; #removing leading spaces
			$_=~ s/\s+$//; #removing leading spaces

			if ($_ =~ /^\d/) {

				my @splitline = split /\s+/, $_;
				my @field1 = split /\./, $splitline[0];

				$hash_raw{$field1[0]}{$splitline[1]}{$field1[1]} = $_;
			}
	}
	close(INPUT);
}

sub convert_genscan{

	my @list_exons;

	foreach my $genscan_id ( sort { $a <=> $b } keys %hash_raw){

		my $gene_id = create_uniq_id("gene");
		my $mrna_id = create_uniq_id("mrna");
		@list_exons=();

		foreach my $tag ( keys %{$hash_raw{$genscan_id} } ){
			# first create EXONS
			#print $tag."\n";
			if ( 	grep( /^$tag/, @list_exon_types ) ){
				foreach my $nb (keys %{$hash_raw{$genscan_id}{$tag} } ){

					my @splitline = split /\s+/, $hash_raw{$genscan_id}{$tag}{$nb};

					my $start=$splitline[3];
					my $end=$splitline[4];
					if($start > $end){
						my $tmp_start = $start;
						$start = $end;
						$end = $tmp_start;
					}
					# phase and frame from genscan cannot be used
					my $strand=$splitline[2];
					my $id = "fake";
					my $score = $splitline[$#splitline];

					my $feature = Bio::SeqFeature::Generic->new(-seq_id => $seq_id,
																											-source_tag => "genscan",
																											-primary_tag => "exon",
																											-start => $start,
																											-end => $end ,
																											-frame => ".",
																											-strand =>$strand,
																											-score =>$score,
																											- tag => {'ID' => $id, 'Parent' => $mrna_id}
																											) ;
					push @list_exons, $feature;
				}
			}
			else{
				#Prom = Promoter (TATA box / initation site)
				#PlyA = poly-A signal (consensus: AATAAA)

			}
		}
		@list_exons = sort {$a->start <=> $b->start} @list_exons;
		if (@list_exons) {

			# CREATE GENE
			my $gene = clone($list_exons[0]);
			create_or_replace_tag($gene, 'ID', $gene_id);
			$gene->remove_tag('Parent');
			$gene->primary_tag('gene');
			$gene->end($list_exons[$#list_exons]->end());

			# CREATE MRNA
			my $mrna = clone($gene);
			create_or_replace_tag($mrna, 'Parent', $gene_id);
			create_or_replace_tag($mrna, 'ID', $mrna_id);
			$mrna->primary_tag('mRNA');


			#print gene
			$gffout->write_feature($gene);
			#print mrna
			$gffout->write_feature($mrna);
			#print exons
			my $strand = $list_exons[0]->strand;

			foreach my $exon (@list_exons){
				my $ordered_exon_id = create_uniq_id("exon");
				create_or_replace_tag($exon, 'ID', $ordered_exon_id);
				$gffout->write_feature($exon);
			}
			#-------------- print cds --------------
			my @cds;
			foreach my $exon (@list_exons){
				my $cds_id = create_uniq_id("cds");
				my $cds = clone($exon);
				create_or_replace_tag($cds, 'ID', $cds_id);
				$cds->primary_tag('CDS');
				push @cds, $cds;
			}

			if(($strand eq "+") or ($strand eq "1")){
				@cds=sort {$a->start <=> $b->start} @cds;
			}else{
				@cds=sort {$b->start <=> $a->start} @cds;
			}

			# compute the phase. Assuming it always start at 0. No fragmented prediction
			my $phase = 0;
			foreach my $cds_feature ( @cds) {
				$cds_feature->frame($phase);
				my $cds_length=$cds_feature->end-$cds_feature->start +1;
				$phase=(3-(($cds_length-$phase)%3))%3; #second modulo allows to avoid the frame with 3. Instead we have 0.
			}

			foreach my $cds ( sort {$a->start <=> $b->start} @cds){
				$gffout->write_feature($cds);
			}
			#-------------- END print cds --------------
		}
	}
}

sub create_uniq_id{
	my ( $tag ) = @_;

	my $uniqID;

	if(! exists_keys(\%hash_uniqID,($tag) ) ){
		$uniqID=$tag."_1";
		$hash_uniqID{$tag}=1;
	}
	else{
		$hash_uniqID{$tag}++;
		$uniqID=$tag."_".$hash_uniqID{$tag};
	}

	return $uniqID;
}

__END__

genscan format Explanation

Gn.Ex : gene number, exon number (for reference)
Type  : Init = Initial exon (ATG to 5' splice site)
        Intr = Internal exon (3' splice site to 5' splice site)
        Term = Terminal exon (3' splice site to stop codon)
        Sngl = Single-exon gene (ATG to stop)
        Prom = Promoter (TATA box / initation site)
        PlyA = poly-A signal (consensus: AATAAA)
S     : DNA strand (+ = input strand; - = opposite strand)
Begin : beginning of exon or signal (numbered on input strand)
End   : end point of exon or signal (numbered on input strand)
Len   : length of exon or signal (bp)
Fr    : reading frame (a forward strand codon ending at x has frame x mod 3)
Ph    : net phase of exon (exon length modulo 3)
I/Ac  : initiation signal or 3' splice site score (tenth bit units)
Do/T  : 5' splice site or termination signal score (tenth bit units)
CodRg : coding region score (tenth bit units)
P     : probability of exon (sum over all parses containing exon)
Tscr  : exon score (depends on length, I/Ac, Do/T and CodRg scores)

Comments

The SCORE of a predicted feature (e.g., exon or splice site) is a
log-odds measure of the quality of the feature based on local sequence
properties. For example, a predicted 5' splice site with
score > 100 is strong; 50-100 is moderate; 0-50 is weak; and
below 0 is poor (more than likely not a real donor site).

The PROBABILITY of a predicted exon is the estimated probability under
GENSCAN's model of genomic sequence structure that the exon is correct.
This probability depends in general on global as well as local sequence
properties, e.g., it depends on how well the exon fits with neighboring
exons.  It has been shown that predicted exons with higher probabilities
are more likely to be correct than those with lower probabilities.

=head1 NAME

agat_convert_genscan2gff.pl

=head1 DESCRIPTION

The script takes a genscan file as input, and will translate it in gff format.
The genscan format is described here: http://genome.crg.es/courses/Bioinformatics2003_genefinding/results/genscan.html
/!\ vvv Known problem vvv /!\
You must have submited only DNA sequence, wihtout any header!!
Indeed the tool expects only DNA sequences and does not crash/warn if an header
is submited along the sequence.
e.g If you have an header ">seq" s-e-q are seen as the 3 first nucleotides of the sequence.
Then all prediction location are shifted accordingly.
(checked only on the online version http://argonaute.mit.edu/GENSCAN.html. I don't
know if there is the same pronlem elsewhere.)
/!\ ^^^ Known problem ^^^^ /!\

=head1 SYNOPSIS

    agat_convert_genscan2gff.pl --genscan infile.bed [ -o outfile ]
    agat_convert_genscan2gff.pl -h

=head1 OPTIONS

=over 8

=item B<--genscan> or B<-g>

Input genscan bed file that will be convert.

=item B<--seqid>

String - Sequence ID. [default: unknown]

=item B<--verbose> or B<-v>

add verbosity

=item B<-o> , B<--output> , B<--out> , B<--outfile> or B<--gff>

Output GFF file. If no output file is specified, the output will be
written to STDOUT.

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
