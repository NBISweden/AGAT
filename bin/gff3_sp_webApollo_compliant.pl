#!/usr/bin/env perl

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use Agat::Omniscient;
use Bio::Tools::GFF;

my $header = get_agat_header();
my $start_run = time();
my $opt_gfffile;
my $opt_output;
my $opt_help = 0;

# OPTION MANAGMENT
if ( !GetOptions( 'g|gff=s' => \$opt_gfffile,
                  'o|output=s'      => \$opt_output,

                  'h|help!'         => \$opt_help ) )
{
    pod2usage( { -message => 'Failed to parse command line',
                 -verbose => 1,
                 -exitval => 1 } );
}

if ($opt_help) {
    pod2usage( { -verbose => 2,
                 -exitval => 2 } );
}

if (! defined($opt_gfffile) ){
    pod2usage( {
           -message => "\nAt least 1 parameter is mandatory:\nInput reference gff file (-g).\n\n".
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

                #####################
                #     MAIN          #
                #####################

######################
### Parse GFF input #
my ($hash_omniscient, $hash_mRNAGeneLink) = slurp_gff3_file_JD({ input => $opt_gfffile
                                                              });
print ("GFF3 file parsed\n");

########
# Transform thing needed for webapollo.
webapollo_compliant($hash_omniscient);

#############
# Print result
print_omniscient($hash_omniscient, $gffout); #print gene modified

my $end_run = time();
my $run_time = $end_run - $start_run;
print "Job done in $run_time seconds\n";
__END__

=head1 NAME

gff3_webApollo_compliant.pl -
This script aim to remove useless/problematic information for webapollo, change some feeaturee type to avoid problem whem loading them into webapollo, and optimize some attribute for a nice displaying.
=head1 SYNOPSIS

    ./gff3_webApollo_compliant.pl -g infile.gff [ -o outfile ]
    ./gff3_webApollo_compliant.pl --help

=head1 OPTIONS

=over 8

=item B<-g>, B<--gff> or B<-ref>

Input GFF3 file that will be read (and sorted)

=item B<-o> or B<--output>

Output GFF file.  If no output file is specified, the output will be
written to STDOUT.

=item B<-h> or B<--help>

Display this helpful text.

=back

=cut
