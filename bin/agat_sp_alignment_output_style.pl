#!/usr/bin/env perl

use strict;
use warnings;
use AGAT::AGAT;

my $header = get_agat_header();
my $start_run = time();
my ( $opt, $usage, $config ) = AGAT::AGAT::describe_script_options(
    $header,
    [ 'gff|g=s', 'Input GTF/GFF file', { required => 1 } ],
);

my $opt_gfffile = $opt->gff;

my $log;
if ( my $log_name = $config->{log_path} ) {
    open( $log, '>', $log_name )
      or die "Can not open $log_name for printing: $!";
    dual_print( $log, $header, 0 );
}

######################
# Manage output file #
my $gffout = prepare_gffout( $config, $config->{output} );

                #####################
                #     MAIN          #
                #####################

######################
### Parse GFF input #
my ($hash_omniscient, $hash_mRNAGeneLink) =  slurp_gff3_file_JD({
                                                               input => $opt_gfffile,
                                                               config => $config
                                                               });
dual_print( $log, "GFF3 file parsed\n", $config->{verbose} );

###
# Print result
print_omniscient_as_match( {omniscient => $hash_omniscient, output => $gffout} );

my $end_run = time();
my $run_time = $end_run - $start_run;
dual_print( $log, "Job done in $run_time seconds\n", $config->{verbose} );
__END__

=head1 NAME

agat_sp_alignment_output_style.pl

=head1 DESCRIPTION

The script takes a normal gtf/gff annotation format file and convert it
to gff3 alignment format. It means it add a structure of match / match_part
as relationship between the different features.

=head1 SYNOPSIS

    agat_sp_alignment_output_style.pl -g infile.gff [ -o outfile ]
    agat_sp_alignment_output_style.pl --help

=head1 OPTIONS

=over 8

=item B<-g>, B<--gff> or B<-ref>

Input GTF/GFF file.

=item B<-v>

Verbose option to see the warning messages when parsing the gff file.

=item B<-o> or B<--output>

Output GFF file.  If no output file is specified, the output will be
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
