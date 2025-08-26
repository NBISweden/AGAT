#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use AGAT::AGAT;

my $header = get_agat_header();
my ( $opt, $usage, $config ) = AGAT::AGAT::describe_script_options(
    $header,
    [ 'gff|g|file|input=s', 'Input reference gff file', { required => 1 } ],
    [ 'feature_type|ft=s',  'Top feature type of the group', { default => 'gene' } ],
    [ 'interval|i=i',       'Number of top features per file', { default => 1000, callbacks => { positive => sub { shift > 0 or die "Interval must be positive" } } } ],
);

my $opt_gff      = $opt->gff;
my $feature_type = $opt->feature_type;
my $interval     = $opt->interval;
my $opt_output   = $config->{output};

my $log;
if ( my $log_name = $config->{log_path} ) {
    open( $log, '>', $log_name ) or die "Can not open $log_name for printing: $!";
    dual_print( $log, $header, 0 );
}
my $opt_verbose = $config->{verbose};

my $start_run = time();

my $format = $config->{force_gff_input_version};
if ( !$format ) { $format = select_gff_format($opt_gff); }
my $ref_in = AGAT::BioperlGFF->new( -file => $opt_gff, -gff_version => $format );

if ( -d $opt_output ) {
    warn "The output directory <$opt_output> already exists.\n" if $opt_verbose;
    exit;
} else {
    my ( $path, $ext );
    ( $opt_output, $path, $ext ) = fileparse( $opt_output, qr/\.[^.]*$/ );
    dual_print( $log, "Creating the $opt_output folder\n");
    mkdir $opt_output;
}

dual_print( $log,
    "I will split the file into files containing $interval group of feature. The top feature of the group of feature is currently defined by <$feature_type>.\n");

my $startP = time;
my $nbLine = `wc -l < $opt_gff`;
$nbLine =~ s/ //g;
chomp $nbLine;
dual_print( $log, "$nbLine line to process...\n");
my $line_cpt = 0;

my $count_feature = 0;
my $count_file    = 1;
my ( $file_name, $path, $ext ) = fileparse( $opt_gff, qr/\.[^.]*$/ );

my $gffout = prepare_gffout( $config, $opt_output . "/" . $file_name . "_" . $count_file . ".gff" );

while ( my $feature = $ref_in->next_feature() ) {
    $line_cpt++;

    if ( $feature->primary_tag eq $feature_type ) {
        if ( $count_feature == $interval ) {
            close $gffout;
            $count_file++;
            $gffout =
              prepare_gffout( $config, $opt_output . "/" . $file_name . "_" . $count_file . ".gff" );
            $count_feature = 0;
        }
        $count_feature++;
    }
    $gffout->write_feature($feature);

    if ( ( 30 - ( time - $startP ) ) < 0 ) {
        my $done = ( $line_cpt * 100 ) / $nbLine;
        $done = sprintf( '%.0f', $done );
        dual_print( $log, "\rProgression : $done % processed.\n");
        $startP = time;
    }
}
close $gffout;

my $end_run  = time();
my $run_time = $end_run - $start_run;
dual_print( $log, "Job done in $run_time seconds\n");

__END__

=head1 NAME

agat_sq_split.pl

=head1 DESCRIPTION

split gff3 file into several files.
By default we create files containing 1000 genes and all sub-features associated.
GFF3 input file must be sequential.

=head1 SYNOPSIS

    agat_sq_split.pl -i <input file> -o <output file>
    agat_sq_split.pl --help

=head1 OPTIONS

=over 8

=item B<--gff>, B<--file> or B<--input>

STRING: Input GTF/GFF file.

=item B<-i> or B<--interval>
Integer.  Number of group of feature to include in each file. 1000 by default.

=item B<--ft> or B<--feature_type>
The top feature of the feature group. By default "gene".

=item B<-o> or B<--output>

STRING: Output file.  If no output file is specified, the output will be written to STDOUT. The result is in tabulate format.

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

