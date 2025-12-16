#!/usr/bin/env perl

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use AGAT::AGAT;

start_script();
my $header = get_agat_header();
# -------------------------------- LOAD OPTIONS --------------------------------
my $opt_output;
my $opt_gff;
my $gtf_version;
my $help;

# OPTION MANAGEMENT: partition @ARGV into shared vs script options via library
my ($shared_argv, $script_argv) = split_argv_shared_vs_script(\@ARGV);

# Parse script-specific options from its own list
my $script_parser = Getopt::Long::Parser->new;
$script_parser->configure('bundling','no_auto_abbrev');
if ( !$script_parser->getoptionsfromarray(
        $script_argv,
    'h|help!'                => \$help,
    'gff|gtf|i=s'            => \$opt_gff,
    'gtf_version=s'          => \$gtf_version,
    'outfile|output|o|out=s' => \$opt_output,
    ) ) {
    pod2usage( { -message => 'Failed to parse command line',
                 -verbose => 1,
                 -exitval => 1 } );
}

# Print Help and exit
if ($help) {
    pod2usage( { -message => "$header\n",
                 -verbose => 99,
                 -exitval => 0 } );
}

if ( ! (defined($opt_gff)) ){
    pod2usage( {
           -message => "$header\nAt least 1 parameters is mandatory:\nInput gff/gtf file (--gff or --gtf).\n\n",
           -verbose => 0,
           -exitval => 1 } );
}

# Parse shared options without pass_through for strong type errors. CPU and config are handled there.
my ($shared_opts) = parse_shared_options($shared_argv);

# --- Load config file into global CONFIG ---
initialize_agat({config_file_in => ( $shared_opts->{config} ), input => $opt_gff, shared_opts => $shared_opts });

# ------------------------------------------------------------------------------

# check GTF versions
if ($gtf_version){
    my @gtf_version_list = (1, 2, 2.1, 2.2, 2.5, 3, "relax");
    my %gtf_version_hash = map { $_ => 1 } @gtf_version_list;
    if(! exists_keys (\%gtf_version_hash, ("$gtf_version") ) ) {
        die "$gtf_version is not a valid GTF version. Please choose one among this list: @gtf_version_list\n";
    }
    dual_print1 "GTF version $gtf_version selected by command line interface.\n";
} else {
    $gtf_version = $CONFIG->{gtf_output_version};
    dual_print1 "GTF version $gtf_version selected from the agat config file.\n";
}

# Update config
$CONFIG->{"gtf_output_version"}=$gtf_version;
$CONFIG->{"output_format"}="gtf";

# Manage output file 
my $gffout = prepare_gffout( $opt_output );

######################
### Read gff input file.
my ($hash_omniscient) = slurp_gff3_file_JD({ input => $opt_gff });
dual_print1 "converting to GTF$gtf_version\n";

# Now print  omniscient
print_omniscient( {omniscient => $hash_omniscient, output => $gffout} );

# --- final messages ---
end_script();

#################################### methods ####################################

__END__

=head1 NAME

agat_convert_sp_gff2gtf.pl

=head1 DESCRIPTION

The script aims to convert any GTF/GFF file into a proper GTF file.
Full information about the format can be found here: https://agat.readthedocs.io/en/latest/gxf.html
You can choose among 7 different GTF types (1, 2, 2.1, 2.2, 2.5, 3 or relax).
Depending the version selected the script will filter out the features that are not accepted.
For GTF2.5 and 3, every level1 feature (e.g nc_gene pseudogene) will be converted into
gene feature and every level2 feature (e.g mRNA ncRNA) will be converted into
transcript feature.
Using the "relax" option you will produce a GTF-like output keeping all
original feature types (3rd column). No modification will occur e.g. mRNA to transcript. 

To be fully GTF compliant all feature have a gene_id and a transcript_id attribute.
The gene_id	is unique identifier for the genomic source of the transcript, which is
used to group transcripts into genes.
The transcript_id is a unique identifier for the predicted transcript,
which is used to group features into transcripts.

=head1 SYNOPSIS

    agat_convert_sp_gff2gtf.pl --gff infile.gff [ -o outfile ]
    agat_convert_sp_gff2gtf -h

=head1 OPTIONS

=over 8

=item B<--gff>, B<--gtf> or B<-i>

Input GFF/GTF file that will be read

=item B<--gtf_version>
version of the GTF output (1,2,2.1,2.2,2.5,3 or relax). Default value from AGAT config file (relax for the default config). The script option has the higher priority.

relax: all feature types are accepted.

GTF3 (9 feature types accepted): gene, transcript, exon, CDS, Selenocysteine, start_codon, stop_codon, three_prime_utr and five_prime_utr

GTF2.5 (8 feature types accepted): gene, transcript, exon, CDS, UTR, start_codon, stop_codon, Selenocysteine

GTF2.2 (9 feature types accepted): CDS, start_codon, stop_codon, 5UTR, 3UTR, inter, inter_CNS, intron_CNS and exon

GTF2.1 (6 feature types accepted): CDS, start_codon, stop_codon, exon, 5UTR, 3UTR

GTF2 (4 feature types accepted): CDS, start_codon, stop_codon, exon

GTF1 (5 feature types accepted): CDS, start_codon, stop_codon, exon, intron

=item B<-o>, B<--output>, B<--out> or B<--outfile>

Output GTF file. If no output file is specified, the output will be
written to STDOUT.

=item B<-h> or B<--help>

Display this helpful text.

=back

=head1 SHARED OPTIONS

Shared options are defined in the AGAT configuration file and can be overridden via the command line for this script only.
Common shared options are listed below; for the full list, please refer to the AGAT agat_config.yaml.

=over 8

=item B<--config>

String - Path to a custom AGAT configuration file.  
By default, AGAT uses `agat_config.yaml` from the working directory if present, otherwise the default file shipped with AGAT
(available locally via `agat config --expose`).

=item B<--cpu>, B<--core>, B<--job> or B<--thread>

Integer - Number of parallel processes to use for file input parsing (via forking).

=item B<-v> or B<--verbose>

Integer - Verbosity, choice are 0,1,2,3,4. 0 is quiet, 1 is normal, 2,3,4 is more verbose. Default 1.

=back

=head1 FEEDBACK

For questions, suggestions, or general discussions about AGAT, please use the AGAT community forum:
https://github.com/NBISweden/AGAT/discussions

=head1 BUG REPORTING

Bug reports should be submitted through the AGAT GitHub issue tracker:
https://github.com/NBISweden/AGAT/issues

=cut

AUTHOR - Jacques Dainat
