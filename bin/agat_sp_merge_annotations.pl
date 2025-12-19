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

start_script();
my $header = get_agat_header();
# -------------------------------- LOAD OPTIONS --------------------------------

my $outfile = undef;
my @opt_files;
my $file2 = undef;
my $opt_help= 0;

# OPTION MANAGEMENT: partition @ARGV into shared vs script options via library
my ($shared_argv, $script_argv) = split_argv_shared_vs_script(\@ARGV);

# Parse script-specific options from its own list
my $script_parser = Getopt::Long::Parser->new;
$script_parser->configure('bundling','no_auto_abbrev');
if ( ! $script_parser->getoptionsfromarray(
    $script_argv,
    'h|help!'                 => \$opt_help,
    'gff|f=s'                 => \@opt_files,
    'output|out|o=s'          => \$outfile,
)) {
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
my $shared_opts = parse_shared_options($shared_argv);
initialize_agat({ config_file_in => ( $shared_opts->{config} ), input => $opt_files[0], shared_opts => $shared_opts });

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
  dual_print1 "\nTotal raw data of files together:\n";
  info_omniscient($hash_omniscient);
}

# Now all the feature are in the same omniscient
# We have to check the omniscient to merge overlaping genes together. Identical isoforms will be removed
dual_print1 "\nNow merging overlaping loci, and removing identical isoforms:\n";
merge_overlap_loci($hash_omniscient);

dual_print1 "\nfinal result:\n";
info_omniscient($hash_omniscient);

########
# Print results
print_omniscient( {omniscient => $hash_omniscient, output => $gffout} );

# --- final messages ---
end_script();

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

=item B<--gff> or B<-f> <file>

Input GTF/GFF file(s). You can specify a folder containing GFF3 files with the format .gff or GTF files with .gtf format . You can also specify as much file you want like so: -f file1 -f file2 -f file3

=item  B<-o>, B<--out> or B<--output> <file>

Output file to create (default GFF3 - see config to modify output format).
If no output file is specified, the output will be written to STDOUT.


=item B<--help> or B<-h>

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
