#!/usr/bin/env perl

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use AGAT::AGAT;

start_script();
my $header = get_agat_header();
# -----------------------------------------------------------------------------------------------
my $opt_gfffile;
my $opt_output;
my $opt_help = 0;

# OPTION MANAGMENT
my @copyARGV=@ARGV;
my ($shared_argv, $script_argv) = split_argv_shared_vs_script(\@ARGV);
my $script_parser = Getopt::Long::Parser->new;
$script_parser->configure('bundling','no_auto_abbrev');
if ( !$script_parser->getoptionsfromarray(
    $script_argv,
    'g|gff=s'     => \$opt_gfffile,
    'o|out|output=s'  => \$opt_output,
    'h|help!'     => \$opt_help ) )
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
           -message => "$header\nAt least 1 parameter is mandatory:\nInput reference gff file (-g).\n\n".
           "Ouptut is optional. Look at the help documentation to know more.\n",
           -verbose => 0,
           -exitval => 1 } );
}

my ($shared_opts) = parse_shared_options($shared_argv);
initialize_agat({ config_file_in => $shared_opts->{config}, input => $opt_gfffile, shared_opts => $shared_opts });
# -----------------------------------------------------------------------------------------------

######################
# Manage output file #
my $gffout = prepare_gffout( $opt_output );

                #####################
                #     MAIN          #
                #####################

######################
### Parse GFF input #
my ($hash_omniscient) =  slurp_gff3_file_JD({ input => $opt_gfffile });

#######################
# Convert FULL standard gff3 to ensembl gff type
convert_omniscient_to_ensembl_style($hash_omniscient);

###
# Print result
print_omniscient( {omniscient => $hash_omniscient, output => $gffout} );
# --- final messages ---
end_script();
# -----------------------------------------------------------------------------------------------
__END__

=head1 NAME

agat_sp_ensembl_output_style.pl

=head1 DESCRIPTION

This script takes a normal gff3 annotation format file and convert it to gff3
like ensembl format.

=head1 SYNOPSIS

    agat_sp_ensembl_output_style.pl -g infile.gff [ -o outfile ]
    agat_sp_ensembl_output_style.pl --help

=head1 OPTIONS

=over 8

=item B<-g>, B<--gff> or B<-ref> <file>

Input GTF/GFF file.

=item B<-o> or B<--output> <file>

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
