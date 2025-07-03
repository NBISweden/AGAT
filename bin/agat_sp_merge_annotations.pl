#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long;
use Pod::Usage;
use List::MoreUtils qw(uniq);
use AGAT::AGAT;
use File::Spec;
use File::Glob ':glob';

my $header = get_agat_header();
my $config;
my $cpu;
my $outfile = undef;
my @opt_files;
my $file2 = undef;
my $opt_help= 0;

if ( !GetOptions(
    'c|config=s'             => \$config,
                    'thread|threads|cpu|cpus|core|cores|job|jobs=i' => \$cpu,
    "h|help"                 => \$opt_help,
    "gff|f=s"                => \@opt_files,
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

my @expanded_files;
foreach my $file_or_dir (@opt_files) {
    if (-d $file_or_dir) {
        push @expanded_files, bsd_glob(File::Spec->catfile($file_or_dir, '*.{gff,gtf}'));
    } else {
        push @expanded_files, $file_or_dir;
    }
}
@opt_files = @expanded_files;

if ( ! @opt_files or (@opt_files and ($#opt_files < 1) ) ){
    pod2usage( {
           -message => "$header\nAt least 2 files are mandatory:\n --gff file1 --gff file2\n\n",
           -verbose => 0,
           -exitval => 1 } );
}

# --- Manage config ---
initialize_agat({ config_file_in => $config, input => $opt_files[0] });
$CONFIG->{cpu} = $cpu if defined($cpu);

######################
# Manage output file #
my $gffout = prepare_gffout( $outfile);

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>     MAIN     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

######################
### Parse GFF input #

my $file1 = shift @opt_files;
my ($hash_omniscient) = slurp_gff3_file_JD({ input => $file1 });
info_omniscient($hash_omniscient);

#Add the features of the other file in the first omniscient. It takes care of name to not have duplicates
foreach my $next_file (@opt_files){
  my ($hash_omniscient2) = slurp_gff3_file_JD({ input => $next_file });

  info_omniscient($hash_omniscient2);

  #merge annotation is taking care of Uniq name. Does not look if mRNA are identic or so one, it will be handle later.
  merge_omniscients($hash_omniscient, $hash_omniscient2);
  print ("\nTotal raw data of files together:\n");
  info_omniscient($hash_omniscient);
}

# Now all the feature are in the same omniscient
# We have to check the omniscient to merge overlaping genes together. Identical isoforms will be removed
print ("\nNow merging overlaping loci, and removing identical isoforms:\n");
merge_overlap_loci($hash_omniscient);


print ("\nfinal result:\n");
info_omniscient($hash_omniscient);

########
# Print results
print_omniscient( {omniscient => $hash_omniscient, output => $gffout} );

__END__

=head1 NAME

agat_sp_merge_annotations.pl

=head1 DESCRIPTION

This script merge different gff annotation files in one.
It uses the AGAT parser that takes care of duplicated names and fixes other oddities met in those files.

=head1 SYNOPSIS

    agat_sp_merge_annotations.pl --gff infile1 --gff infile2 --out outFile
    agat_sp_merge_annotations.pl --gff /path/to/folder/with/gff --out outFile
    agat_sp_merge_annotations.pl --help

=head1 OPTIONS

=over 8

=item B<--gff> or B<-f>

Input GTF/GFF file(s). You can specify a folder containing GFF3 files with the format .gff or GTF files with .gtf format . You can also specify as much file you want like so: -f file1 -f file2 -f file3

=item  B<--out>, B<--output> or B<-o>

Output gff3 file where the gene incriminated will be write.

=item B<-thread>, B<threads>, B<cpu>, B<cpus>, B<core>, B<cores>, B<job> or B<jobs>

Integer - Number of parallel processes to use for file input parsing (via forking).

=item B<-c> or B<--config>

String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any, 
otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

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
