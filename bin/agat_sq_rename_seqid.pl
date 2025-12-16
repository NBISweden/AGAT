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
# ---------------------------- OPTIONS ----------------------------
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
		'gff=s'           => \$input_gff,
		'o|output=s'      => \$outputFile,
		'tsv=s'           => \$input_tsv,
		'csv!'            => \$csv,
		'h|help!'         => \$opt_help )  )
{
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
initialize_agat({ config_file_in => $shared_opts->{config}, input => $input_gff, shared_opts => $shared_opts });

# Manage Output
my $gffout = prepare_gffout( $outputFile);

# Manage GFF Input
my $format = $CONFIG->{force_gff_input_version};
if(! $format ){ $format = select_gff_format($input_gff); }
my $inputfh = open_maybe_gz($input_gff);
my $gff_in = AGAT::BioperlGFF->new(-fh => $inputfh, -gff_version => $format);

# Manage tsv input
open(INPUT, "<", $input_tsv) or die ("$!\n");

# Open tsv file for reading
my %tsv;
while (<INPUT>) {
	chomp;

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

	$tsv{$splitline[0]} = $splitline[1];
}

while (my $feature = $gff_in->next_feature() ) {

	if(exists_keys(\%tsv,  ($feature->seq_id() ) ) ){
		$feature->seq_id($tsv{$feature->seq_id()});
	}
	$gffout->write_feature($feature);
}

# write_fasta in asked and any
write_fasta($gffout, $gff_in);

# --- final messages ---
end_script();

# ---------------------------- FUNCTIONS ----------------------------

__END__

=head1 NAME

agat_sq_rename_seqid.pl

=head1 DESCRIPTION

The script aims to modify seqid (1st column) of a GFF/GTF file efficiently. 
Indeed, when the number of chromosomes or scaffolding is large, 
replacement using e.g. sed command can be time-consuming. 
You must provide a file (tsv or csv) without header and with 
one renaming information by line: The first value is the original sequence identifier (1st column of the GFF/GTF file), 
the second value is the new sequence identifier to use.


number of chromosomes or scaffolding is large, sed replacement is time-consuming

=head1 SYNOPSIS

    agat_sq_rename_seqid.pl --gff input.gff --tsv input.tsv [ -o output.gff3 ]
    agat_sq_rename_seqid.pl --help

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
