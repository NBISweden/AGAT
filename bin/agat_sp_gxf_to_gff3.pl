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

# Print Help and exit
if ($opt_help) {
    pod2usage( { -verbose => 99,
                 -exitval => 0,
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

agat_sp_gxf_to_gff3.pl

=head1 DESCRIPTION

This script fixes and/or standardizes a gtf/gff file into full sorted gff3.
It will be read by the Omniscient parser, that will first detect automtically
which GFF parser to use from bioperl (GFF1,GFF2,GFF3) to read the lines.
Then the Omniscient parser removes duplicate features, fixes duplicated IDs,
add missing ID and/or Parent attributes, deflates factorized attributes
(attributes with several parents are duplicated with uniq ID), add missing features
when possible (e.g add exon if only CDS described, add UTR if CDS and exon described),
fix feature locations (e.g check exon is embeded in th parent features mRNA, gene), etc...
All AGAT's scripts with the _sp_ prefix use the same parser, before to perform suplement tasks.
whith the script you can tuned the Omniscient parser behaviour. I.e you can decide
to merge loci that have an overlap at their CDS features ( Only one top feature
is kept (gene), and the mRNA features become isoforms). This is not activated by
default in case you are working on a prokaryote annotation that often have overlaping
loci.
The Omniscient parser defines relationship between features using 3 levels.
e.g Level1=gene; Level2=mRNA,tRNA; Level3=exon,cds,utr.
The feature type information is stored within the 3rd column of a GTF/GFF file.
Which level a feature is part of is crucial for the parser. This information
is stored by default in a json file coming with the tool. We have implemented the
most common feature types met in gff. If a feature type in your file is not imprelement
the parser will not handle it and inform you. You could easily inform the parser how
to handle it (level1,level2 or level3) by adding it in the corresponding file. How to
access the json files? Easy just use the --expose option and the json files will appear in
the workoing folder. If they are present, the Omniscient parser use the json files
from the working direcrtory by default.


Omniscient parser phylosophy:

 Parse by Parent/child relationship
   ELSE Parse by a comomn tag  (an attribute value shared by feature that must be grouped together.
   By default we are using locus_tag and gene_id as locus tag, but you can specify the one of your choice
     ELSE Parse sequentially (features are grouped in a bucket, and the bucket change at each level2 feature met, and bucket(s) are linked to the first l1 top feature met)

=head1 SYNOPSIS

    agat_sp_gxf_to_gff3.pl -g infile.gff [ -o outfile ]
    agat_sp_gxf_to_gff3.pl --help

=head1 OPTIONS

=over 8

=item B<-g>, B<--gff> or B<-ref>

Input GTF/GFF file.

=item B<-c> or B<--ct>

When the features doesn't have Parent/ID relationships, the parser will try to group
features using a common/shared attribute (i.e. a locus tag.). By default locus_tag and gene_id.
You can provide another specific common/shared attribute using this option.

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

Verbose option. To modify vefbosity. Default is 1. 0 is quiet, 2 and 3 are increasing verbosity.

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
