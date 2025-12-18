#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use IO::File;
use List::MoreUtils qw(uniq);
use AGAT::AGAT;

start_script();
my $header = get_agat_header();
# ---------------------------- OPTIONS ----------------------------
my $gff = undef;
my $opt_output = undef;
my $opt_help= 0;

# ---------------------------- OPTIONS ----------------------------
my ($shared_argv, $script_argv) = split_argv_shared_vs_script(\@ARGV);

my $script_parser = Getopt::Long::Parser->new;
$script_parser->configure('bundling','no_auto_abbrev');
if ( ! $script_parser->getoptionsfromarray(
    $script_argv,
    'h|help!'          => \$opt_help,
    'o|out|output=s'   => \$opt_output,
    'gff|f=s'          => \$gff ))
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

if ( ! (defined($gff)) ){
    pod2usage( {
           -message => "$header\nAt least 1 parameter is mandatory:\nInput reference gff file (--gff) \n\n",
           -verbose => 0,
           -exitval => 2 } );
}

# --- Manage config ---
my $shared_opts = parse_shared_options($shared_argv);
initialize_agat({ config_file_in => ( $shared_opts->{config} ), input => $gff, shared_opts => $shared_opts });

# ----------------------------------------------------------------------------

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    PARAMS    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

my $gffout = prepare_gffout( $opt_output );

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>     MAIN     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

######################
### Parse GFF input #
my ($hash_omniscient) = slurp_gff3_file_JD({ input => $gff });
### END Parse GFF input #
#########################

# clean omniscient to remove isoforms
my ($nb_iso_removed_cds,  $nb_iso_removed_exon) = remove_shortest_isoforms($hash_omniscient);

# print omniscientNew containing only the longest isoform per gene
print_omniscient( {omniscient => $hash_omniscient, output => $gffout} );

dual_print1 $nb_iso_removed_cds." L2 isoforms with CDS removed (shortest CDS)\n";
dual_print1 $nb_iso_removed_exon." L2 isoforms wihtout CDS removed (Either no isoform has CDS, we removed those with shortest concatenated exons, or at least one isoform has CDS, we removed those wihtout)\n";

# --- final messages ---
end_script();

# ----------------------------------------------------------------------------


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

=item B<--gff> or B<-f> <file>

GTF/GFF file.

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
