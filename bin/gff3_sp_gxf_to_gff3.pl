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
my $opt_merge = undef;
my $opt_comonTag=undef;
my $opt_verbose = 1;
my $opt_no_check = undef;
my $opt_output;
my $opt_expose_feature_levels = undef;
my $opt_help = 0;
my $opt_version_input = undef;
my $opt_version_output = 3;

# OPTION MANAGMENT
my @copyARGV=@ARGV;
if ( !GetOptions( 'g|gff=s'         => \$opt_gfffile,
                  'c|ct=s'          => \$opt_comonTag,
                  'v=i'             => \$opt_verbose,
                  'o|output=s'      => \$opt_output,
                  'efl|expose!'      => \$opt_expose_feature_levels,
                  'nc|no_check!'      => \$opt_no_check,
                  'gff_version_input|gvi=f'   => \$opt_version_input,
                  'gff_version_output|gvo=f'   => \$opt_version_output,
                  'ml|merge_loci!'     => \$opt_merge,
                  'h|help!'         => \$opt_help ) )
{
    pod2usage( { -message => 'Failed to parse command line',
                 -verbose => 1,
                 -exitval => 1 } );
}

if ($opt_help) {
    pod2usage( { -verbose => 2,
                 -exitval => 2,
                 -message => "$header\n" } );
}

if (! defined($opt_gfffile) ){
    pod2usage( {
           -message => "$header\nAt least 1 parameter is mandatory:\nInput reference gff file (-g).\n\n".
           "Ouptut is optional. Look at the help documentation to know more.\n",
           -verbose => 0,
           -exitval => 1 } );
}

#print perl version
print "-------------------------------------------------------------------------------\n";
print "This script is being run by perl ".$^V."\n";
print "Bioperl location being used: ".substr($INC{"Bio/Tools/GFF.pm"}, 0 , -12)."\n";
print "-------------------------------------------------------------------------------\n";
#############################
# check version input value #
check_version($opt_version_input);
#check_version($opt_version_output);

######################
# Manage output file #
my $gffout;
if ($opt_output) {
  open(my $fh, '>', $opt_output) or die "Could not open file '$opt_output' $!";
  $gffout= Bio::Tools::GFF->new(-fh => $fh, -gff_version => $opt_version_output );
  }
else{
  $gffout = Bio::Tools::GFF->new(-fh => \*STDOUT, -gff_version => $opt_version_output);
}

                #####################
                #     MAIN          #
                #####################

######################
### Parse GFF input #
my ($hash_omniscient, $hash_mRNAGeneLink) = slurp_gff3_file_JD({
                                                               input => $opt_gfffile,
                                                               locus_tag => $opt_comonTag,
                                                               gff_version => $opt_version_input,
                                                               verbose => $opt_verbose,
                                                               merge_loci => $opt_merge,
                                                               no_check => $opt_no_check,
                                                               expose_feature_levels => $opt_expose_feature_levels
                                                               });
print ("GFF3 file parsed\n");

###
# Print result

print_omniscient($hash_omniscient, $gffout); #print gene modified

my $end_run = time();
my $run_time = $end_run - $start_run;
print "usage: $0 @copyARGV\n";
print "Job done in $run_time seconds\n";


sub check_version{
  my ($version) = @_;
  if($version and ($version != 1 and $version != 2 and $version != 3)){
    print "Gff version accepted is 1,2 or 3. $version is not a correct value.\n";
    exit;
  }
}

__END__
=head1 NAME

gff3_sp_gxf_to_gff3.pl -
This script read and print a gff file. It will be read by GFF3::Omniscient parser that will look for duplicate features, duplicate IDs and will print the features sorted.
The result is written to the specified output file, or to STDOUT.

=head1 SYNOPSIS

    ./gff3_sp_gxf_to_gff3.pl -g infile.gff [ -o outfile ]
    ./gff3_sp_gxf_to_gff3.pl --help

=head1 OPTIONS

=over 8

=item B<-g>, B<--gff> or B<-ref>

Input GFF3 file that will be read (and sorted)

=item B<-c> or B<--ct>

When the gff file provided is not correcly formated and features are linked to each other by a comon tag (by default locus_tag), this tag can be provided to parse the file correctly.

=item B<--efl> or B<--expose>

If you want to see, add or modified the feature relationships you will have to use this option.
It will copy past in you working directory the json files used to define the relation between feature types and their level organisation.
Typical level organisation: Level1 => gene; Level2 => mRNA; level3 => exon,cds,utrs
If you get warning from the Omniscient parser that a feature relationship is not defined, you can provide information about it within the exposed json files.
Indeed, if the json files exists in your working directory, they will be used by default.

=item B<--ml> or B<--merge_loci>

Merge loci parameter, default deactivated. You turn on the parameter if you want to merge loci into one locus when they overlap.
(at CDS level for mRNA, at exon level for other level2 features. Strand has to be the same). Prokaryote can have overlaping loci so it should not use it for prokaryote annotation.
In eukaryote, loci rarely overlap. Overlaps could be due to error in the file, mRNA can be merged under the same parent gene if you acticate the option.

=item B<-v>

Verbose option. To modify vefbosity. Default 1. 0 is quiet, 2 and 3 are increasing verbosity.

=item B<--nc> or B<--no_check>

To deacticate all check that can be performed by the parser (e.g fixing UTR, exon, coordinates etc...)

=item B<-o> or B<--output>

Output GFF file.  If no output file is specified, the output will be
written to STDOUT.

=item B<--gvi> or B<--gff_version_input>

If you don't want to use the autodection of the gff/gft version you give as input, you can force the tool to use the parser of the gff version you decide to use: 1,2,2.5 or 3. Remind: 2.5 is suposed to be gtf.

=item B<--gvo> or B<--gff_version_output>

If you don't want to use the autodection of the gff/gft version you give as input, you can force the tool to use the parser of the gff version you decide to use: 1,2,2.5 or 3. Remind: 2.5 is suposed to be gtf.

=item B<-h> or B<--help>

Display this helpful text.

=back

=cut
