#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use IO::File;
use List::MoreUtils qw(uniq);
use AGAT::AGAT;

my $header = get_agat_header();
my ($config, $gff, $opt_output, $opt_genomeSize, $opt_plot, $opt_help);

my $common = parse_common_options() || {};
$config     = $common->{config};
$opt_output = $common->{output};
$opt_help   = $common->{help};

if ( !GetOptions(
    "gff|f=s"         => \$gff))

{
    pod2usage( { -message => "Failed to parse command line",
                 -verbose => 1,
                 -exitval => 1 } );
}

# Print Help and exit
if ($opt_help) {
    pod2usage( { -verbose => 99,
                 -exitval => 0,
                 -message => "$header\n" } );
}

if ( ! (defined($gff)) ){
    pod2usage( {
           -message => "$header\nAt least 1 parameter is mandatory:\nInput reference gff file (--gff) \n\n",
           -verbose => 0,
           -exitval => 2 } );
}

# --- Manage config ---
$config = get_agat_config({config_file_in => $config});

my $log;
if ($config->{log}) {
  my ($file) = $0 =~ /([^\/]+)$/;
  my $log_name = $file . ".agat.log";
  open($log, '>', $log_name) or die "Can not open $log_name for printing: $!";
  dual_print($log, $header, 0);
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    PARAMS    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

my $gffout = prepare_gffout($config, $opt_output);

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>     MAIN     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

######################
### Parse GFF input #
dual_print($log, "Reading file $gff\n");
my ($hash_omniscient, $hash_mRNAGeneLink) = slurp_gff3_file_JD({ input => $gff,
                                                                 config => $config
                                                              });
dual_print($log, "Parsing Finished\n");
### END Parse GFF input #
#########################

# clean omniscient to remove isoforms
my ($nb_iso_removed_cds,  $nb_iso_removed_exon) = remove_shortest_isoforms($hash_omniscient);

# print omniscientNew containing only the longest isoform per gene
print_omniscient( {omniscient => $hash_omniscient, output => $gffout} );

dual_print($log, $nb_iso_removed_cds." L2 isoforms with CDS removed (shortest CDS)\n");
dual_print($log, $nb_iso_removed_exon." L2 isoforms wihtout CDS removed (Either no isoform has CDS, we removed those with shortest concatenated exons, or at least one isoform has CDS, we removed those wihtout)\n");

# END STATISTICS #
##################
dual_print($log, "Done\n");

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
