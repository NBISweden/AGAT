#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use warnings;
use Pod::Usage;
use Getopt::Long;
use IO::File ;
use Bio::SeqIO;
use AGAT::AGAT;

start_script();
my $header = get_agat_header();
# -------------------------------- LOAD OPTIONS --------------------------------
my $input_gff;
my $input_tsv;
my $outputFile;
my $csv;
my $opt_help = 0;

# OPTION MANAGEMENT: partition @ARGV into shared vs script options via library
my ($shared_argv, $script_argv) = split_argv_shared_vs_script(\@ARGV);

# Parse script-specific options from its own list
my $script_parser = Getopt::Long::Parser->new;
$script_parser->configure('bundling','no_auto_abbrev');
if ( ! $script_parser->getoptionsfromarray(
		$script_argv,
		'gff=s'        => \$input_gff,
		'o|output=s'   => \$outputFile,
		'tsv=s'        => \$input_tsv,
		'csv!'         => \$csv,
		'h|help!'      => \$opt_help,
	) ) {
	pod2usage( { -message => 'Failed to parse command line',
				 -verbose => 1,
				 -exitval => 1 } );
}

if ($opt_help) {
    pod2usage( { -verbose => 99,
                 -exitval => 0,
                 -message => "$header\n" } );
}

if (! $input_gff or ! $input_tsv){
	pod2usage( {  -message => "$header\nAt least 2 input file are mandatory:\n".
					  "--gff input.gff\n--tsv input.tsv",
					  -verbose => 0,
					  -exitval => 1 } );
}

# Parse shared options (CPU, config, etc.)
my ($shared_opts) = parse_shared_options($shared_argv);

# --- Manage config ---
initialize_agat({ config_file_in => ( $shared_opts->{config} ), input => $input_gff, shared_opts => $shared_opts });

# ------------------------------------------------------------------------------

# Manage Output
my $gffout = prepare_gffout( $outputFile );

# Manage GFF Input
my $format = $CONFIG->{force_gff_input_version};
if(! $format ){ $format = select_gff_format($input_gff); }
dual_print1 "Reading $input_gff using format GFF$format\n";
my $inputfh = open_maybe_gz($input_gff);
my $gff_in = AGAT::BioperlGFF->new(-fh => $inputfh, -gff_version => $format);

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
		dual_print2 "$nb_header headers\n";
		my $cpt = 0;
		foreach my $header_title (@splitline){
			$header{$cpt++} = $header_title;
		}
	}
	else{
		my $nb_column = scalar @splitline;
		if($nb_column != $nb_header) { dual_print1 "Number of header ($nb_header) different to number of columm ($nb_column) line $line\n"; }
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
						dual_print2 "Value $tsv{$id}{$att} already exists for attribute $att in feature with ID $id\n";
					}
					#add attribute
					else{
						dual_print2 "Attribute $att exists for feature with ID $id, we add the new value $tsv{$id}{$att} to it.\n";
						$feature->add_tag_value($att, $tsv{$id}{$att});
					}
				}
				# new attribute
				else{
					dual_print2 "New attribute $att with value $tsv{$id}{$att} added for feature with ID $id.\n";
					$feature->add_tag_value($att, $tsv{$id}{$att});
				}
			}
		}
	}
	$gffout->write_feature($feature);
}

# print fasta in asked and any
write_fasta($gffout, $gff_in);

# --- final messages ---
end_script();

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

=item B<-o> or B<--output>

STRING: Output file. If no output file is specified, the output will be written
to STDOUT. The result is in tabulate format.


=item B<--help> or B<-h>

Display this helpful text.

=back

=head1 SHARED OPTIONS

Shared options are defined in the AGAT configuration file and can be overridden via the command line for this script only.
Common shared options are listed below; for the full list, please refer to the AGAT agat_config.yaml.
Note: For _sq_ scripts, only the following options are supported: verbose, output_format, gff_output_version, gtf_output_version, progress_bar, and tabix.

=over 8

=item B<--config>

String - Path to a custom AGAT configuration file.  
By default, AGAT uses `agat_config.yaml` from the working directory if present, otherwise the default file shipped with AGAT
(available locally via `agat config --expose`).

=item B<-v> or B<--verbose>

Integer - Verbosity, choice are 0,1,2,3,4. 0 is quiet, 1 is normal, 2,3,4 is more verbose. Default 1.

=back

=head1 FEEDBACK

For questions, suggestions, or general discussions about AGAT, please use the AGAT community forum:
https://github.com/NBISweden/AGAT/discussions

=head1 BUG REPORTING

Bug reports should be submitted through the AGAT GitHub issue tracker:
https://github.com/NBISweden/AGAT/issues

=cut

AUTHOR - Jacques Dainat
