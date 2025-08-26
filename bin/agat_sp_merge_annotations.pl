#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Pod::Usage;
use List::MoreUtils qw(uniq);
use AGAT::AGAT;
use File::Spec;
use File::Glob ':glob';

my $header = get_agat_header();
my ( $opt, $usage, $config ) = AGAT::AGAT::describe_script_options(
    $header,
    [ 'gff|f=s@', 'Input GTF/GFF files or directories', { required => 1 } ],
);

my @opt_files = @{ $opt->gff };
my $outfile   = $opt->out;
my $opt_verbose = $config->{verbose};

my $log;
if ( my $log_name = $config->{log_path} ) {
    open( $log, '>', $log_name ) or die "Can not open $log_name for printing: $!";
}
dual_print( $log, $header, 3);

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

######################
# Manage output file #
my $gffout = prepare_gffout($config, $outfile);

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>     MAIN     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

######################
### Parse GFF input #

my $file1 = shift @opt_files;
my ($hash_omniscient, $hash_mRNAGeneLink) = slurp_gff3_file_JD({ input => $file1,
                                                                 config => $config
                                                              });
dual_print($log, "$file1 GFF3 file parsed\n");
info_omniscient($hash_omniscient, $log, $opt_verbose);

#Add the features of the other file in the first omniscient. It takes care of name to not have duplicates
foreach my $next_file (@opt_files){
  my ($hash_omniscient2, $hash_mRNAGeneLink2) = slurp_gff3_file_JD({ input => $next_file,
	                                                                   config => $config
                                                                  });
  dual_print($log, "$next_file GFF3 file parsed\n");
  info_omniscient($hash_omniscient2, $log, $opt_verbose);

  #merge annotation is taking care of Uniq name. Does not look if mRNA are identic or so one, it will be handle later.
  merge_omniscients($hash_omniscient, $hash_omniscient2);
  dual_print($log, "\nTotal raw data of files together:\n");
  info_omniscient($hash_omniscient, $log, $opt_verbose);
}

# Now all the feature are in the same omniscient
# We have to check the omniscient to merge overlaping genes together. Identical isoforms will be removed
dual_print($log, "\nNow merging overlaping loci, and removing identical isoforms:\n");
merge_overlap_loci(undef, $hash_omniscient, $hash_mRNAGeneLink, undef);


dual_print($log, "\nfinal result:\n");
info_omniscient($hash_omniscient, $log, $opt_verbose);

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
