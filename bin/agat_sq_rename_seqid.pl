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

my $header = get_agat_header();
my $config;
my $start_run = time();
my $input_gff;
my $input_tsv;
my $outputFile;
my $verbose;
my $csv;
my $opt_help = 0;

Getopt::Long::Configure ('bundling');
if ( !GetOptions (  'gff=s'           => \$input_gff,
                    'o|output=s'      => \$outputFile,
	                'tsv=s'           => \$input_tsv,
                    'csv!'            => \$csv,
			        'v|verbose!'      => \$verbose,
                    'c|config=s'      => \$config,
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

# --- Manage config ---
initialize_agat({ config_file_in => $config, input => $input_gff });

# Manage Output
my $gffout = prepare_gffout( $outputFile);

# Manage GFF Input
my $format = $CONFIG->{force_gff_input_version};
if(! $format ){ $format = select_gff_format($input_gff); }
my $gff_in = AGAT::BioperlGFF->new(-file => $input_gff, -gff_version => $format);

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

my $end_run = time();
my $run_time = $end_run - $start_run;
print "Job done in $run_time seconds\n";

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
