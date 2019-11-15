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
use Agat::Omniscient;

my $header = get_agat_header();
my $outfile = undef;
my $gff     = undef;
my $blast   = undef;
my $help;

if ( !GetOptions(   "help"      => \$help,
                    "gff=s"     => \$gff,
                    "blast=s"   => \$blast,
                    "outfile=s" => \$outfile ))

{
    pod2usage( { -message => 'Failed to parse command line',
                 -verbose => 1,
                 -exitval => 1 } );
}

# Print Help and exit
if ($help) {
    pod2usage( { -verbose => 2,
                 -exitval => 2,
                 -message => "$header\n" } );
}

if ( ! (defined($gff)) or !(defined($blast)) ){
    pod2usage( {
           -message => "$header\nAt least 2 parameter is mandatory:\nInput reference gff file (--gff) and Input blast file (--blast)\n\n",
           -verbose => 0,
           -exitval => 1 } );
}


# Open Input gff3 file #
my $ref_in = Bio::Tools::GFF->new(-file =>$gff , -gff_version => 3);

# Open Output files #
my $out;
if ($outfile) {
    open(my $fh, '>', $outfile) or die "Could not open file '$outfile' $!";
    $out= Bio::Tools::GFF->new(-fh => $fh, -gff_version => 3 );
}
else {
     $out = Bio::Tools::GFF->new(-fh => \*STDOUT, -gff_version => 3);
}


#### MAIN ####

# Read killlist #
my $killlist = parse_blast($blast);

### Parse GFF input #
print ("Parse file $gff\n");
my ($hash_omniscient, $hash_mRNAGeneLink) = slurp_gff3_file_JD({ input => $gff
                                                              });
print ("$gff file parsed\n");

# Remove all mRNA specified by the kill-list from their (gene-) parents.
remove_omniscient_elements_from_level2_ID_list ($hash_omniscient, $killlist);

# Write the remaining things to output
print_omniscient($hash_omniscient, $out); #print gene modified in file

      #########################
      ######### END ###########
      #########################
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
    print "$nbremove gene will be removed !\n";

    return \@answer;
} ## end sub parse_blast


=head1 NAME

gff_filter_by_mrnaBlastValue.pl
ancient name gff_filter_by_mrna_id.pl gff_filter_by_mrnaBlastValue.pl
The script aims to remove from a gff file all the sequence that have a similarity over THRESHOLD with another sequence (will keep only one).
This is typically useful when creating a list of mRNA to use to train abinitio gene finder.
A reciprocal blast of the sequences need to have been performed prior to the use of this script in order to get the blastp input file.

=head1 SYNOPSIS

    ./gff_filter_by_mrnaBlastValue.pl --gff=infile.gff3 --blast blastfile --outfile outFile
    ./gff_filter_by_mrnaBlastValue.pl --help

=head1 OPTIONS

=over 8

=item B<--gff>

Input GFF3 file correponding to gene build.

=item B<--blast>

The list of the all-vs-all blast file (outfmt 6, blastp)

=item  B<--outfile>

The name of the output file. By default the output is the standard output.

=item B<--help> or B<-h>

Display this helpful text.

=back

=cut
