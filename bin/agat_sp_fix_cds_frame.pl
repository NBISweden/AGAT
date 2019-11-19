#!/usr/bin/env perl

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use Bio::Tools::GFF;
use Bio::DB::Fasta;
use Agat::Omniscient;


my $header = get_agat_header();
my $start_run = time();
my $opt_fasta = undef;
my $opt_gfffile;
my $opt_verbose;
my $codonTable=0;
my $opt_output;
my $opt_help = 0;

# OPTION MANAGMENT
if ( !GetOptions( 'g|gff=s'     => \$opt_gfffile,
                  'o|output=s'  => \$opt_output,
                  "fasta|fa=s" => \$opt_fasta,
                  "v|vebose!" => \$opt_verbose,
                  "table|codon|ct=i" => \$codonTable,
                  'h|help!'         => \$opt_help ) )
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

if (! defined($opt_gfffile) or ! defined($opt_fasta)){
    pod2usage( {
           -message => "$header\nAt least 2 parameters are mandatory:\nInput reference gff file (-g) and Input fasta file (--fasta).\n\n".
           "Ouptut is optional. Look at the help documentation to know more.\n",
           -verbose => 0,
           -exitval => 1 } );
}

######################
# Manage output file #

my $gffout;
if ($opt_output) {
  $opt_output=~ s/.gff//g;
  open(my $fh, '>', $opt_output.".gff") or die "Could not open file '$opt_output' $!";
  $gffout= Bio::Tools::GFF->new(-fh => $fh, -gff_version => 3 );
  }
else{
  $gffout = Bio::Tools::GFF->new(-fh => \*STDOUT, -gff_version => 3);
}

if($codonTable<0 and $codonTable>25){
  print "$codonTable codon table is not a correct value. It should be between 0 and 25 (0,23 and 25 can be problematic !)\n";
}
else{
  print "We will use the codon table ".$codonTable.". If it is not what you want please stop the tool and use the --table option. \n";
}

                #####################
                #     MAIN          #
                #####################

####################
# index the genome #
my $db = Bio::DB::Fasta->new($opt_fasta);
print ("Genome fasta parsed\n");

######################
### Parse GFF input #
my ($hash_omniscient, $hash_mRNAGeneLink) = slurp_gff3_file_JD({ input => $opt_gfffile
                                                            });
print ("GFF3 file parsed\n");

###
# Fix frame
fil_cds_frame($hash_omniscient, $db, $opt_verbose);

###
# Print result
print_omniscient($hash_omniscient, $gffout); #print gene modified

my $end_run = time();
my $run_time = $end_run - $start_run;
print "Job done in $run_time seconds\n";
__END__

=head1 NAME

agat_sp_fix_cds_frame.pl

=head1 DESCRIPTION

This script aims to fix the cds phases.

=head1 SYNOPSIS

    agat_sp_fix_cds_frame.pl -g infile.gff -f fasta [ -o outfile ]
    agat_sp_fix_cds_frame.pl --help

=head1 OPTIONS

=over 8

=item B<-g>, B<--gff> or B<-ref>

Input GTF/GFF file.

=item B<-fa> or B<--fasta>

Genome fasta file.

=item B<--ct>, B<--codon> or B<--table>

Codon table to use. 0 By default.

=item B<-v> or B<--verbose>

Add verbosity.

=item B<-o> or B<--output>

Output GFF file.  If no output file is specified, the output will be
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
