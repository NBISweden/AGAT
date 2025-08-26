#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use warnings;
use Pod::Usage;
use Getopt::Long::Descriptive;
use IO::File;
use Bio::SeqIO;
use AGAT::AGAT;

my $header = get_agat_header();
my $start_run = time();
my @copyARGV = @ARGV;
my ( $opt, $usage, $config ) = AGAT::AGAT::describe_script_options( $header,
    [ 'gff=s', 'Input GFF file', { required => 1 } ],
    [ 'tsv=s', 'Input TSV file', { required => 1 } ],
    [ 'csv!',  'TSV file is comma-separated' ],
);

my $input_gff = $opt->gff;
my $input_tsv = $opt->tsv;
my $csv       = $opt->csv;

my $log;
if ( my $log_name = $config->{log_path} ) {
    open( $log, '>', $log_name ) or die "Can not open $log_name for printing: $!";
    dual_print( $log, $header,  3 );
}

# Manage Output
my $outputFile = $config->{output};
my $gffout = prepare_gffout($config, $outputFile);

# Manage GFF Input
my $format = $config->{force_gff_input_version};
if(! $format ){ $format = select_gff_format($input_gff); }
dual_print($log, "Reading $input_gff using format GFF$format\n");
my $gff_in = AGAT::BioperlGFF->new(-file => $input_gff, -gff_version => $format);

# Manage tsv input
open(INPUT, "<", $input_tsv) or die ("$!\n");

# Open tsv file for reading
my $line = 0;
my %tsv;
my %header;
my $nb_header = undef;
while (<INPUT>) {
	chomp;
	$line++;

	$_=~ s/^\s+//; #removing leading spaces
	$_=~ s/\s+$//; #removing trailing spaces

	# split line
	my @splitline;
	if ($csv){
		@splitline = split /,/, $_;
	}
  else{
	  @splitline = split /\t/, $_; # split at tabulation
  }

	if ($line == 1){
		$nb_header = scalar @splitline;
             dual_print($log, "$nb_header headers\n");
		my $cpt = 0;
		foreach my $header_title (@splitline){
			$header{$cpt++} = $header_title;
		}
	}
	else{
		my $nb_column = scalar @splitline;
             if($nb_column != $nb_header) { dual_print($log, "Number of header ($nb_header) different to number of columm ($nb_column) line $line\n"); }
		for(my $i = 1; $i <= $#splitline; $i++){
			$tsv{lc($splitline[0])}{$header{$i}} = $splitline[$i];
		}
	}
}

my $cpt_line = 0;
while (my $feature = $gff_in->next_feature() ) {
	$cpt_line++;
	if($feature->has_tag($header{"0"})){
		my $id = lc($feature->_tag_value($header{"0"}));
		if (exists_keys(\%tsv,($id))){
			foreach my $att ( sort keys %{$tsv{$id}} ){
				#attribute already exists
				if($feature->has_tag($att)){
					my @originalvalues = $feature->get_tag_values($att);
					if (grep( /$tsv{$id}{$att}/, @originalvalues)){
                                           dual_print($log, "Value $tsv{$id}{$att} already exists for attribute $att in feature with ID $id\n", 2);
					}
					#add attribute
					else{
                                           dual_print($log, "Attribute $att exists for feature with ID $id, we add the new value $tsv{$id}{$att} to it.\n", 2);
						$feature->add_tag_value($att, $tsv{$id}{$att});
					}
				}
				# new attribute
				else{
                                   dual_print($log, "New attribute $att with value $tsv{$id}{$att} added for feature with ID $id.\n", 2);
					$feature->add_tag_value($att, $tsv{$id}{$att});
				}
			}
		}
	}
	$gffout->write_feature($feature);
}

my $end_run = time();
my $run_time = $end_run - $start_run;
dual_print($log, "Job done in $run_time seconds\n");

close $log if $log;

__END__

=head1 NAME

agat_sq_add_attributes_from_tsv.pl

=head1 DESCRIPTION

The purpose of the script is to add information from a tsv/csv file to the attributes of a gff file (9th column).
e.g. an attribute looks like this in a GFF3 file: tag=value1,value2 
The first line of the tsv/csv file must contain the headers (corresponding to an attribute tag in the GFF/GTF file),
while the other lines contain the values (corresponding to an attribute value in the GFF/GTF file).
The first column is used to synchronize information between the tsv file and the GFF/GTF file. In other words, 
it's used to determine which feature we're going to add attributes to.
The other columns will be added as attribute in the GFF/GTF file. The header becomes the tag for the new attribute, 
and the value is that defined for the corresponding feature line. 
(If the tag already exists, we append the value only if the value doesn't already exist).

--- example ---

* input.tsv:
ID	annot_type1
gene1	anot_x
cds1	anot_y

* input gff:
chr1	irgsp	gene	1000	2000	.	+	.	ID=gene1
chr1	irgsp	CDS	2983	3268	.	+	.	ID=cds1

* output.gff:
chr1	irgsp	gene	1000	2000	.	+	.	ID=gene1;annot_type1=anot_x
chr1	irgsp	CDS	2983	3268	.	+	.	ID=cds1;annot_type1=anot_y

--- example2 ---

* input.tsv:
gene_id	annot_type1
gene1	anot_x
cds1	anot_y

* input gff:
chr1	irgsp	gene	1000	2000	.	+	.	gene_id=gene1
chr1	irgsp	CDS	2983	3268	.	+	.	gene_id=cds1

* output.gff:
chr1	irgsp	gene	1000	2000	.	+	.	gene_id=gene1;annot_type1=anot_x
chr1	irgsp	CDS	2983	3268	.	+	.	gene_id=cds1;annot_type1=anot_y

=head1 SYNOPSIS

    agat_sq_add_attributes_from_tsv.pl --gff input.gff --tsv input.tsv [ -o output.gff3 ]
    agat_sq_add_attributes_from_tsv.pl --help

=head1 OPTIONS

=over 8

=item B<--gff>

STRING: Input GTF/GFF file.

=item B<--tsv>

STRING: Input tsv file

=item B<--csv>

BOLEAN: Inform the script that the tsv input file is actually a csv (coma-separated).

=item B<-v> or B<--verbose>

BOLEAN: Add verbosity

=item B<-o> or B<--output>

STRING: Output file. If no output file is specified, the output will be written
to STDOUT. The result is in tabulate format.

=item B<-c> or B<--config>

String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any, 
otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

=item B<--help> or B<-h>

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
