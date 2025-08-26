#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use warnings;
use Pod::Usage;
use IO::File ;
use Bio::SeqIO;
use AGAT::AGAT;

my $header = get_agat_header();
my ( $opt, $usage, $config ) = AGAT::AGAT::describe_script_options(
    $header,
    [ 'gff=s', 'Input GTF/GFF file', { required => 1 } ],
    [ 'tsv=s', 'Input tsv file',      { required => 1 } ],
    [ 'delimiter=s',
      'Field delimiter for TSV/CSV file',
      { default   => "\t",
        callbacks => {
            allowed => sub {
                my $d = shift;
                ( $d eq "\t" || $d eq "," )
                  or die 'delimiter must be "\t" or ","';
            }
        } } ],
    [ 'csv!' => 'Input file is comma-separated', { implies => { delimiter => ',' } } ],
);

my $input_gff = $opt->gff;
my $input_tsv = $opt->tsv;
my $outfile   = $config->{output};
my $delimiter = $opt->delimiter;

my $log;
if ( my $log_name = $config->{log_path} ) {
    open( $log, '>', $log_name )
      or die "Can not open $log_name for printing: $!";
    dual_print( $log, $header,  3 );
}

my $start_run = time();

# Manage Output
my $gffout = prepare_gffout( $config, $outfile );

# Manage GFF Input
my $format = $config->{force_gff_input_version};
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
        my @splitline = split /\Q$delimiter\E/, $_;

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
dual_print( $log, "Job done in $run_time seconds\n");

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

    agat_sq_rename_seqid.pl --gff input.gff --tsv mapping.tsv [ -o output.gff3 ]
    agat_sq_rename_seqid.pl --gff input.gff --tsv mapping.csv --csv
    agat_sq_rename_seqid.pl --gff input.gff --tsv mapping.tsv --delimiter ","
    agat_sq_rename_seqid.pl --help

=head1 OPTIONS

=over 8

=item B<--gff>

STRING: Input GTF/GFF file.

=item B<--tsv>

STRING: Input tsv file

=item B<--csv>

Boolean: convenience flag setting the delimiter to a comma.

=item B<--delimiter>

STRING: Column separator in the mapping file. Allowed values are "\t" and ",".
Defaults to a tab character. Using C<--csv> automatically sets this to a comma.

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
