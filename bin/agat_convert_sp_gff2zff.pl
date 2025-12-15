#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Bio::DB::Fasta;
use File::Basename;
use Bio::SeqIO;
use AGAT::AGAT;

start_script();
my $header = get_agat_header();
# -------------------------------- LOAD OPTIONS --------------------------------
my $outfile = undef;
my $opt_gff = undef;
my $model_id = -1;
my $fasta = undef;
my $help;

# OPTION MANAGEMENT: partition @ARGV into shared vs script options via library
my ($shared_argv, $script_argv) = split_argv_shared_vs_script(\@ARGV);

# Parse script-specific options from its own list
my $script_parser = Getopt::Long::Parser->new;
$script_parser->configure('bundling','no_auto_abbrev');
if ( ! $script_parser->getoptionsfromarray(
		$script_argv,
		'h|help'                  => \$help,
		'gff=s'                   => \$opt_gff,
		'fasta=s'                 => \$fasta,
		'outfile|output|out|o=s'  => \$outfile,
	) ) {
	pod2usage({
		-message => 'Failed to parse command line.',
		-verbose => 1,
		-exitval => 1
	});
}
# Print Help and exit
if ($help) {
	pod2usage({ -message => "$header",
				-verbose => 99,
				-exitval => 0 });
}

if ( ! defined($opt_gff) or ! defined($fasta) ){
	pod2usage({
		   -message => "$header\nAt least 2 parameters are mandatory:\n  Input gff file (--gff)\nInput fasta file (--fasta)\n",
		   -verbose => 0,
		   -exitval => 1 });
}

# Parse shared options
my ($shared_opts) = parse_shared_options($shared_argv);

# --- Manage config ---
initialize_agat({ config_file_in => ($shared_opts->{config}), input => $opt_gff, shared_opts => $shared_opts });
# -----------------------------------------------------------------------------------------------

## Manage output file
my $zffout;
my $fastaout;
if ($outfile) {
	my ($outfile_prefix,$path,$ext) = fileparse($outfile,qr/\.[^.]*/);
	my $outfile_zff = $outfile_prefix.".ann";
	open(my $fh, '>', $outfile_zff) or die "Could not open file '$outfile_zff' $!";
	$zffout=$fh;
	my $outfile_fasta = $outfile_prefix.".dna";
	open(my $fh2, '>', $outfile_fasta) or die "Could not open file '$outfile_fasta' $!";
	$fastaout= Bio::SeqIO->new(-fh => $fh2, -format => 'Fasta' );
}
else{
  $zffout=\*STDOUT ;
	$fastaout = Bio::SeqIO->new(-fh => \*STDOUT, -format => 'Fasta');
}

#### read fasta
my $nbFastaSeq=0;
my $db = Bio::DB::Fasta->new($fasta);
my @ids      = $db->get_all_primary_ids;
my %allIDs; # save ID in lower case to avoid cast problems
foreach my $id (@ids ){$allIDs{lc($id)}=$id;}

### Parse GTF input file
my ($hash_omniscient) = slurp_gff3_file_JD({ input => $opt_gff });

# sort by seq id
my $hash_sortBySeq = gather_and_sort_l1_by_seq_id($hash_omniscient);

#################
# == LEVEL 1 == #
#################
foreach my $seqid (sort { (($a =~ /(\d+)$/)[0] || 0) <=> (($b =~ /(\d+)$/)[0] || 0) } keys %{$hash_sortBySeq}){ # loop over all the feature level1

	#print fasta sequence
	my $seq_id_correct = undef;
	if( exists $allIDs{lc($seqid)}){
		$seq_id_correct = $allIDs{lc($seqid)};
		my $seqObj = Bio::Seq->new( '-format' => 'fasta' , -seq =>  $db->seq($seq_id_correct));
		$seqObj->id($seq_id_correct);
		$fastaout->write_seq($seqObj);
	}
	else{
		warn "$seqid sequence ID not found among the fasta file!\n";
	}

	print $zffout ">".$seqid."\n";

	foreach my $tag_l1 (sort {$a cmp $b} keys %{$hash_omniscient->{'level1'}}){
		foreach my $feature_l1 ( @{$hash_sortBySeq->{$seqid}{$tag_l1}} ){
			my $id_l1 = lc($feature_l1->_tag_value('ID'));


	    #################
	    # == LEVEL 2 == #
	    #################
	    foreach my $tag_l2 ( sort {$a cmp $b} keys %{$hash_omniscient->{'level2'}}){ # tag_l2 = mrna or mirna or ncrna or trna etc...
				if( exists_keys( $hash_omniscient, ('level2', $tag_l2, $id_l1 ) ) ) {
	        foreach my $feature_l2 ( sort {$a->start <=> $b->start} @{$hash_omniscient->{'level2'}{$tag_l2}{$id_l1}}) {

	          #################
	          # == LEVEL 3 == #
	          #################
	          my $level2_ID = lc( $feature_l2->_tag_value('ID') );

						if( exists_keys( $hash_omniscient, ('level3', 'cds', $level2_ID ) ) ) {
							$model_id++;
							my $exon_list_ref = $hash_omniscient->{'level3'}{'cds'}{$level2_ID};

							#deal with single exon
							if ( scalar @{ $exon_list_ref } == 1 ){
								my $start = $exon_list_ref->[0]->start;
								my $end = $exon_list_ref->[0]->end;
								if ( ($exon_list_ref->[0]->strand eq '-') or $exon_list_ref->[0]->strand eq '-1')
								{
									($start, $end) = ($end, $start);
								}
								print $zffout "Esngl\t$start\t$end\tMODEL$model_id\n";
								next;
							}

							#deal multi exons
							# check strand for sorting features
							my $strand = $exon_list_ref->[0]->strand;
							if ( ($strand eq "+" ) or ($strand eq "1" ) ) {
								@{$exon_list_ref} = sort {$a->start <=> $b->start} @{$exon_list_ref};
							}else{
								@{$exon_list_ref} = sort {$b->start <=> $a->start} @{$exon_list_ref};
							}

							my $cpt = 1;
							foreach my $feature_l3 ( @{$exon_list_ref} ) {

								# check strand to deal with location
								my $start = $feature_l3->start;
								my $end = $feature_l3->end;
								if ( ($strand eq '-') or ($strand eq '-1') ){
										($start, $end) = ($end, $start);
								}

								#deal multi exons
								if ( $cpt eq 1 ){ #first exon
									print $zffout "Einit\t$start\t$end\tMODEL$model_id\n";
								}
								elsif ( $cpt eq  scalar @{ $exon_list_ref } ){ #last exon
									print $zffout "Eterm\t$start\t$end\tMODEL$model_id\n";
								}
								else{
									print $zffout "Exon\t$start\t$end\tMODEL$model_id\n";
								}
								$cpt++;
							}
	          }
	        }
	      }
	    }
	  }
	}
}

# --- final messages ---
end_script();

__END__


=head1 NAME

agat_convert_sp_gff2zff.pl

=head1 DESCRIPTION

The script converts GTF/GFF file into zff file a format used by the ab initio
tool SNAP. The script produces a .ann file containing the annotation and .dna
file containing the fasta file. The .ann and .dna are identicaly sorted by
sequence identifier (This is mandatory for usage with fathom).

=head1 SYNOPSIS

    agat_convert_sp_gff2zff.pl --gff file.gff  --fasta file.fasta [ -o outfile ]
    agat_convert_sp_gff2zff.pl --help

=head1 OPTIONS

=over 8

=item B<--gff>

Input GTF/GFF file

=item B<--fasta>

Input fasta file

=item B<--outfile>, B<--out>, B<--output>, or B<-o>

File prefix where will be written the results (e.g. outfile.ann and outfile.dna).
If no output file is specified, the output will be written to STDOUT.

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
