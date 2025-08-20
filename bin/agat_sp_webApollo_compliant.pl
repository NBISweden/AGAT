#!/usr/bin/env perl

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use AGAT::AGAT;

my $header = get_agat_header();
my $config;
my $start_run = time();
my $opt_gfffile;
my $opt_output;
my $opt_help = 0;

my $common = parse_common_options() || {};
$config     = $common->{config};
$opt_output = $common->{output};
my $verbose = $common->{verbose};
$opt_help   = $common->{help};

# OPTION MANAGMENT
if ( !GetOptions( 'g|gff=s' => \$opt_gfffile ) )
{
    pod2usage( { -message => 'Failed to parse command line',
                 -verbose => 1,
                 -exitval => 1 } );
}

# Print Help and exit
if ($opt_help) {
    pod2usage( { -verbose => 99,
                 -exitval => 0,
                 -message => "$header\n" } );
}

if (! defined($opt_gfffile) ){
    pod2usage( {
           -message => "\nAt least 1 parameter is mandatory:\nInput reference gff file (-g).\n\n".
           "Ouptut is optional. Look at the help documentation to know more.\n",
           -verbose => 0,
           -exitval => 1 } );
}

# --- Manage config ---
$config = get_agat_config({config_file_in => $config});

my $log;
my $log_name = get_log_path($common, $config);
open($log, '>', $log_name) or die "Can not open $log_name for printing: $!";
dual_print($log, $header, 0);

######################
# Manage output file #

my $gffout = prepare_gffout($config, $opt_output);

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>     MAIN     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

######################
### Parse GFF input #
my ($hash_omniscient, $hash_mRNAGeneLink) = slurp_gff3_file_JD({ input => $opt_gfffile,
                                                                 config => $config
                                                              });
print ("GFF3 file parsed\n");

########
# Transform thing needed for webapollo.
webapollo_compliant($hash_omniscient);

#############
# Print result
print_omniscient( {omniscient => $hash_omniscient, output => $gffout} );

my $end_run = time();
my $run_time = $end_run - $start_run;
print "Job done in $run_time seconds\n";
__END__

=head1 NAME

agat_sp_webApollo_compliant.pl

=head1 DESCRIPTION

This script aim to remove useless/problematic information for webapollo,
change some featuree type to avoid problem whem loading them into webapollo,
and optimize some attribute for a nice displaying.

=head1 SYNOPSIS

    agat_sp_webApollo_compliant.pl -g infile.gff [ -o outfile ]
    agat_sp_webApollo_compliant.pl --help

=head1 OPTIONS

=over 8

=item B<-g>, B<--gff> or B<-ref>

Input GTF/GFF file.

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
