#!/usr/bin/env perl

# script similar to agat_sp_gxf_to_gff3.pl

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use File::Basename;
use AGAT::AGAT;

start_script();
my $header = get_agat_header();
# -------------------------------- LOAD OPTIONS --------------------------------
my $opt_gfffile;
my $opt_output;
my $opt_help;

# OPTION MANAGEMENT: partition @ARGV into shared vs script options via library
my ($shared_argv, $script_argv) = split_argv_shared_vs_script(\@ARGV);

# Parse script-specific options from its own list
my $script_parser = Getopt::Long::Parser->new;
$script_parser->configure('bundling','no_auto_abbrev');
if ( !$script_parser->getoptionsfromarray(
        $script_argv,
        'g|gxf|gtf|gff=s'          => \$opt_gfffile,
        'o|out|output=s'           => \$opt_output,
        'h|help!'                  => \$opt_help,
    ) ) {
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
           -message => "$header\nAt least 1 parameter is mandatory:\n --gff (Input reference gff file).\n\n".
           "Invoke the help for more information (--help).\n",
           -verbose => 0,
           -exitval => 1 } );
}
# Parse shared options without pass_through for strong type errors. CPU and config are handled there.
my ($shared_opts) = parse_shared_options($shared_argv);

# --- Load config file into global CONFIG ---
initialize_agat({config_file_in => ( $shared_opts->{config} ), input => $opt_gfffile, shared_opts => $shared_opts });

                  

######################
# Manage output file #
my $gffout = prepare_gffout( $opt_output );

                #####################
                #     MAIN          #
                #####################

######################
### Parse GFF input #
my ($hash_omniscient) = slurp_gff3_file_JD({ input  => $opt_gfffile });

###
# Print result
print_omniscient( {omniscient => $hash_omniscient,
                   output => $gffout
                } );

# --- final messages ---
end_script();


__END__
=head1 NAME

agat_convert_sp_gxf2gxf.pl

=head1 DESCRIPTION

This script fixes and/or standardizes any GTF/GFF file into full sorted GTF/GFF file.
It AGAT parser removes duplicate features, fixes duplicated IDs, adds missing ID and/or Parent attributes,
deflates factorized attributes (attributes with several parents are duplicated with uniq ID),
add missing features when possible (e.g. add exon if only CDS described, add UTR if CDS and exon described),
fix feature locations (e.g. check exon is embedded in the parent features mRNA, gene), etc...

All AGAT's scripts with the _sp_ prefix use the AGAT parser, before to perform any supplementary task.
So, it is not necessary to run this script prior the use of any other _sp_ script. 

=head1 SYNOPSIS

    agat_convert_sp_gxf2gxf.pl -g infile.gff [ -o outfile ]
    agat_convert_sp_gxf2gxf.pl --help

=head1 OPTIONS

=over 8

=item B<-g>, B<--gtf>, B<--gff> or B<--gxf> <file>

Input GTF/GFF file. Compressed file with .gz extension is accepted.

=item B<-o>, B<--out> or B<--output> <file>

Output file to create (default GFF3 - see config to modify output format).
If no output file is specified, the output will be written to STDOUT.

=item B<-h> or B<--help>

Display this helpful text.

=back

=head1 SHARED OPTIONS

Shared options are defined in the AGAT configuration file and can be overridden via the command line for this script only.
Common shared options are listed below; for the full list, please refer to the AGAT agat_config.yaml.

=over 8

=item B<--config> <file>

Path to a custom AGAT configuration file.  
By default, AGAT uses `agat_config.yaml` from the working directory if present, otherwise the default file shipped with AGAT
(available locally via `agat config --expose`).

=item B<--cpu>, B<--core>, B<--job> or B<--thread> <int>

Number of parallel processes to use for file input parsing (via forking).

=item B<-v> or B<--verbose> <int>

Verbosity, choice are 0,1,2,3,4. 0 is quiet, 1 is normal, 2,3,4 is more verbose. Default 1.

=back

=head1 FEEDBACK

For questions, suggestions, or general discussions about AGAT, please use the AGAT community forum:
https://github.com/NBISweden/AGAT/discussions

=head1 BUG REPORTING

Bug reports should be submitted through the AGAT GitHub issue tracker:
https://github.com/NBISweden/AGAT/issues

=cut

AUTHOR - Jacques Dainat
