#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long;
use Pod::Usage;
use List::MoreUtils qw(uniq);
use Bio::Tools::GFF;
use Agat::Omniscient;

my $header = get_agat_header();
my $outfile = undef;
my @opt_files;
my $file2 = undef;
my $opt_help= 0;

if ( !GetOptions(
    "help|h" => \$opt_help,
    "gff|f=s" => \@opt_files,
    "output|outfile|out|o=s" => \$outfile))

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

if ( ! @opt_files or (@opt_files and ($#opt_files < 1) ) ){
    pod2usage( {
           -message => "$header\nAt least 2 files are mandatory:\n --gff file1 --gff file2\n\n",
           -verbose => 0,
           -exitval => 1 } );
}

######################
# Manage output file #
my $gffout;
if ($outfile) {
  $outfile=~ s/.gff//g;
open(my $fh, '>', $outfile.".gff") or die "Could not open file '$outfile' $!";
  $gffout= Bio::Tools::GFF->new(-fh => $fh, -gff_version => 3 );
}
else{
  $gffout = Bio::Tools::GFF->new(-fh => \*STDOUT, -gff_version => 3);
}


                #####################
                #     MAIN          #
                #####################


######################
### Parse GFF input #

my $file1 = shift @opt_files;
my ($hash_omniscient, $hash_mRNAGeneLink) = slurp_gff3_file_JD({ input => $file1
                                                              });
print ("$file1 GFF3 file parsed\n");
info_omniscient($hash_omniscient);

#Add the features of the other file in the first omniscient. It takes care of name to not have duplicates
foreach my $next_file (@opt_files){
  my ($hash_omniscient2, $hash_mRNAGeneLink2) = slurp_gff3_file_JD({ input => $next_file
                                                              });
  print ("$next_file GFF3 file parsed\n");
  info_omniscient($hash_omniscient2);

  #merge annotation taking care of Uniq name. Does not look if mRNA are identic or so one, it will be handle later.
  merge_omniscients($hash_omniscient, $hash_omniscient2);
  print ("\n$next_file added we now have:\n");
  info_omniscient($hash_omniscient);
}

# Now all the feature are in the same omniscient
# We have to check the omniscient to merge overlaping genes together and remove the identical ones
($hash_omniscient, $hash_mRNAGeneLink) = slurp_gff3_file_JD({ input => $hash_omniscient,
                                                                 merge_loci => 1
                                                               });
print ("\nfinal result:\n");
info_omniscient($hash_omniscient);

########
# Print results
print_omniscient($hash_omniscient, $gffout);

__END__

=head1 NAME

agat_sp_merge_annotations.pl

=head1 DESCRIPTION

This script merge different gff annotation files in one.
It uses the Omniscient parser that takes care of duplicated names and fixes other oddities met in those files.

=head1 SYNOPSIS

    agat_sp_merge_annotations.pl --gff infile1 --gff infile2 --out outFile
    agat_sp_merge_annotations.pl --help

=head1 OPTIONS

=over 8

=item B<--gff> or B<-f>

Input GTF/GFF file(s). You can specify as much file you want like so: -f file1 -f file2 -f file3

=item  B<--out>, B<--output> or B<-o>

Output gff3 file where the gene incriminated will be write.

=item B<--help> or B<-h>

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
