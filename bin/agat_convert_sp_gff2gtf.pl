#!/usr/bin/env perl

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use AGAT::AGAT;

my $header = get_agat_header();
my $config = get_agat_config();
my $outfile = undef;
my $gff = undef;
my $relax  = undef;
my $gtf_version = 3;
my $verbose = undef;
my $help;


if( !GetOptions(
    "help" => \$help,
    "gff|gtf|i=s" => \$gff,
	"gtf_version=s" => \$gtf_version,
	"verbose|v!" => \$verbose,
    "outfile|output|o|out=s" => \$outfile))
{
    pod2usage( { -message => "Failed to parse command line.",
                 -verbose => 1,
                 -exitval => 1 } );
}

# Print Help and exit
if ($help) {
    pod2usage( { -message => "$header\n",
                 -verbose => 99,
                 -exitval => 0 } );
}

if ( ! (defined($gff)) ){
    pod2usage( {
           -message => "$header\nAt least 1 parameters is mandatory:\nInput gff/gtf file (--gff or --gtf).\n\n",
           -verbose => 0,
           -exitval => 1 } );
}

## Manage output file
my $gtf_out;
if ($outfile) {
open(my $fh, '>', $outfile) or die "Could not open file '$outfile' $!";
  $gtf_out = Bio::Tools::GFF->new(-fh => $fh, -gff_version => 2.5);
}
else{
  $gtf_out = Bio::Tools::GFF->new(-fh => \*STDOUT, -gff_version => 2.5);
}


# check GTF versions
my @gtf_version_list = (1, 2, 2.1, 2.2, 2.5, 3, "relax");
my %gtf_version_hash = map { $_ => 1 } @gtf_version_list;
if(! exists_keys (\%gtf_version_hash, ("$gtf_version") ) ) {
	print "$gtf_version is not a valid GTF version. Please choose one among this list: @gtf_version_list\n"; exit;
}

print "Update config\n";
$config->{"gtf_output_version"}=$gtf_version;
$config->{"output_format"}="gtf";

print "Reading input file\n";
######################
### Parse GFF input #
### Read gff input file.
my ($hash_omniscient, $hash_mRNAGeneLink) = slurp_gff3_file_JD({ input => $gff,
                                                                 config => $config });
print "converting to GTF$gtf_version\n";
# Now print  omniscient
print_omniscient( {omniscient => $hash_omniscient, output => $gtf_out} );


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
version of the GTF output (1,2,2.1,2.2,2.5,3 or relax). Default 3.

relax: all feature types are accepted.

GTF3 (9 feature types accepted): gene, transcript, exon, CDS, Selenocysteine, start_codon, stop_codon, three_prime_utr and five_prime_utr

GTF2.5 (8 feature types accepted): gene, transcript, exon, CDS, UTR, start_codon, stop_codon, Selenocysteine

GTF2.2 (9 feature types accepted): CDS, start_codon, stop_codon, 5UTR, 3UTR, inter, inter_CNS, intron_CNS and exon

GTF2.1 (6 feature types accepted): CDS, start_codon, stop_codon, exon, 5UTR, 3UTR

GTF2 (4 feature types accepted): CDS, start_codon, stop_codon, exon

GTF1 (5 feature types accepted): CDS, start_codon, stop_codon, exon, intron

=item B<-o> , B<--output> , B<--out> , B<--outfile> or B<--gtf>

Output GTF file. If no output file is specified, the output will be
written to STDOUT.

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
