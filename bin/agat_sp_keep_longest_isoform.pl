#!/usr/bin/env perl

use strict;
use warnings;
use IO::File;
use List::MoreUtils qw(uniq);
use AGAT::AGAT;

my $header = get_agat_header();
my @copyARGV = @ARGV;
my ( $opt, $usage, $config ) = AGAT::AGAT::describe_script_options(
    $header,
    [ 'gff|f|ref|reffile=s', 'Input GTF/GFF file', { required => 1 } ],
);

my $gff        = $opt->gff;
my $opt_output = $config->{output};
my $opt_verbose = $config->{verbose};

my $log;
if ( my $log_name = $config->{log_path} ) {
    open( $log, '>', $log_name )
      or die "Can not open $log_name for printing: $!";
    dual_print( $log, $header, 0 );
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    PARAMS    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

my $gffout = prepare_gffout($config, $opt_output);

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>     MAIN     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

######################
### Parse GFF input #
dual_print( $log, "Reading file $gff\n", $opt_verbose );
my ($hash_omniscient, $hash_mRNAGeneLink) = slurp_gff3_file_JD({ input => $gff,
                                                                 config => $config
                                                              });
dual_print( $log, "Parsing Finished\n", $opt_verbose );
### END Parse GFF input #
#########################

# clean omniscient to remove isoforms
my ($nb_iso_removed_cds,  $nb_iso_removed_exon) = remove_shortest_isoforms($hash_omniscient);

# print omniscientNew containing only the longest isoform per gene
print_omniscient( {omniscient => $hash_omniscient, output => $gffout} );

dual_print( $log, $nb_iso_removed_cds . " L2 isoforms with CDS removed (shortest CDS)\n", $opt_verbose );
dual_print( $log, $nb_iso_removed_exon . " L2 isoforms wihtout CDS removed (Either no isoform has CDS, we removed those with shortest concatenated exons, or at least one isoform has CDS, we removed those wihtout)\n", $opt_verbose );

# END STATISTICS #
##################
dual_print( $log, "Done\n", $opt_verbose );

close $log if $log;




__END__

=head1 NAME

agat_sp_keep_longest_isoform.pl

=head1 DESCRIPTION

The script aims to filter isoforms when present. For a locus:
- when all isoforms have CDS we keep the one with the longest CDS.
- when some isoforms have CDS some others not, we keep the one with the longest CDS.
- when none of the isoforms have CDS, we keep the one with the longest concatenated exons. 

=head1 SYNOPSIS

    agat_sp_keep_longest_isoform.pl -gff file.gff  [ -o outfile ]
    agat_sp_keep_longest_isoform.pl --help

=head1 OPTIONS

=over 8

=item B<--gff> or B<-f>

GTF/GFF file.

=item B<--output> or B<-o>

File where will be written the result. If no output file is specified, the output will be written to STDOUT.

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
