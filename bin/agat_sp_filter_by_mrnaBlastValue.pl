#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Pod::Usage;
use Getopt::Long;
use Scalar::Util qw(openhandle);
use Time::Piece;
use Time::Seconds;
use URI::Escape;
use AGAT::AGAT;

start_script();
my $header = get_agat_header();
# ---------------------------- OPTIONS ----------------------------
my $outfile = undef;
my $gff     = undef;
my $blast   = undef;
my $opt_help;

# Partition @ARGV into shared vs script options
my ($shared_argv, $script_argv) = split_argv_shared_vs_script(\@ARGV);

# Parse script-specific options
my $script_parser = Getopt::Long::Parser->new;
$script_parser->configure('bundling','no_auto_abbrev');
if ( ! $script_parser->getoptionsfromarray(
    $script_argv,
    'h|help!'    => \$opt_help,
    'gff=s'      => \$gff,
    'blast=s'    => \$blast,
    'outfile=s'  => \$outfile ) )
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

if ( ! (defined($gff)) or !(defined($blast)) ){
    pod2usage( {
           -message => "$header\nAt least 2 parameter is mandatory:\nInput reference gff file (--gff) and Input blast file (--blast)\n\n",
           -verbose => 0,
           -exitval => 1 } );
}

# Parse shared options and initialize AGAT
my ($shared_opts) = parse_shared_options($shared_argv);
initialize_agat({ config_file_in => ( $shared_opts->{config} ), input => $gff, shared_opts => $shared_opts });

# Open Output files #
my $out = prepare_gffout( $outfile );

#### MAIN ####

# Read killlist #
my $killlist = parse_blast($blast);

### Parse GFF input #
my ($hash_omniscient) = slurp_gff3_file_JD({ input => $gff});


# Remove all mRNA specified by the kill-list from their (gene-) parents.
remove_omniscient_elements_from_level2_ID_list ($hash_omniscient, $killlist);

# Write the remaining things to output
print_omniscient( {omniscient => $hash_omniscient, output => $out} );

      #########################
      ######### END ###########
      #########################

end_script();

#######################################################################################################################
        ####################
         #     methods    #
          ################
           ##############
            ############
             ##########
              ########
               ######
                ####
                 ##


sub parse_blast
{
    my @answer;
    my %moreThanOneTest;
    my %duo_answer;
    my %hashAns;
    my $infile = shift;
    my $cpt2   = 0;
    # This is one way to open a file...
    open( my $IN, '<', $infile ) or
      die "FATAL: Can't open BLAST file: $infile for reading.\n$!\n";

    # Streaming the file, line by line
    while (<$IN>) {
        chomp;
        my $line = $_;

        my @elements = split( "\t", $line );

        my ( $query, $target, $score ) = @elements[ 0 .. 2 ];

        # Matches that we need to remove
        if ( $query ne $target and $score > 80.0 ) {
            ####### <<<<<<<<<<<<<<<<<<<< HERE THE BLAST VALUE CONSIDERED
            my $id      = "$query$target";
            my $idInver = "$target$query";

            if ( ( !exists( $hashAns{$id} ) ) && ( !exists( $hashAns{$idInver} ) ) )
            {    # avoid redundance info
                $hashAns{$id}++;
                $cpt2++;

                # keep the 2 ids We will then remove one randomly
                $duo_answer{$target} = [ $target, $query ];

                $moreThanOneTest{$target}++;
                $moreThanOneTest{$query}++
                  ;    # Allows to detect mRNA present more than 1 times
                       # (In this case they will be selected in priority
                       # during step 3)
            }

        }
    } ## end while (<$IN>)

    # We should close the file to make sure that the transaction
    # finishes cleanly.
    close($IN);

    #print "$cpt2\n";

    # Detect case to remove absolutely to select in a tuple this one if
    # the other we can keep it
    my %caseToAvoid;
    my $cpt = 0;
    foreach my $key ( keys %moreThanOneTest ) {

        if ( $moreThanOneTest{$key} > 1 ) {
            $caseToAvoid{$key}++;
            my $valueUnEsc = uri_unescape($key);
            push (@answer, $valueUnEsc);    # name from blast must be unescape
            $cpt++;
        }
    }
    #print "We will remove $cpt\n";

    ## Step3
    my $cptCount = 0;
    my $removed  = 0;

    # We will keep one of the tuple
    foreach my $key ( keys %duo_answer ) {
        my ( $val1, $val2 ) = @{ $duo_answer{$key} };
        if ( ( !exists( $caseToAvoid{$val1} ) ) and
             ( !exists( $caseToAvoid{$val2} ) ) )
        {    # case remove one randomly
            my $valueUnEsc = uri_unescape($val1);
            push (@answer, $valueUnEsc);    # name from blast must be unescape
            $cptCount++;
        }
    }

    #print "We will removed $cptCount more.\n";
    my $nbremove = @answer;
    dual_print1 "$nbremove gene will be removed !\n";

    return \@answer;
} ## end sub parse_blast

__END__

=head1 NAME

agat_sp_filter_by_mrnaBlastValue.pl

=head1 DESCRIPTION

The script aims to remove from a gff file all the sequence that have a similarity
over THRESHOLD with another sequence (will keep only one).
This is typically useful when creating a list of mRNA to use to train abinitio gene finder.
A reciprocal blast of the sequences need to have been performed prior
to the use of this script in order to get the blastp input file.

=head1 SYNOPSIS

    agat_sp_filter_by_mrnaBlastValue.pl --gff infile.gff --blast blastfile --outfile outFile
    agat_sp_filter_by_mrnaBlastValue.pl --help

=head1 OPTIONS

=over 8

=item B<--gff>

Input GTF/GFF file.

=item B<--blast>

The list of the all-vs-all blast file (outfmt 6, blastp)

=item  B<--outfile>

The name of the output file. By default the output is the standard output.


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
