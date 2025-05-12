#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long;
use Pod::Usage;
use List::MoreUtils qw(uniq);
use AGAT::AGAT;

my $header = get_agat_header();
my $config;
my $cpu;
my $start_run = time();
my $opt_output = undef;
my @opt_files;
my $ref = undef;
my $size_min = 0;
my $opt_help= undef;

# OPTION MANAGMENT
my @copyARGV=@ARGV;
if ( !GetOptions(
    'c|config=s'               => \$config,
                    'thread|threads|cpu|cpus|core|cores|job|jobs=i' => \$cpu,
    "h|help" => \$opt_help,
    "ref|r|i=s" => \$ref,
    "add|a=s" => \@opt_files,
    "size_min|s=i" => \$size_min,
    "output|outfile|out|o=s" => \$opt_output))

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

if (! $ref or ! @opt_files ){
    pod2usage( {
           -message => "$header\nAt least 2 files are mandatory:\n --ref file1 --add file2\n\n",
           -verbose => 0,
           -exitval => 2 } );
}

# --- Manage config ---
initialize_agat({ config_file_in => $config, input => $ref });
$CONFIG->{cpu} = $cpu if defined($cpu);

######################
# Manage output file #
my $gffout = prepare_gffout( $opt_output );

                #####################
                #     MAIN          #
                #####################

######################
### Parse GFF input #

my ($hash_omniscient, $hash_mRNAGeneLink) = slurp_gff3_file_JD({ input => $ref });

info_omniscient($hash_omniscient);

#Add the features of the other file in the first omniscient. It takes care of name to not have duplicates
foreach my $next_file (@opt_files){
  my ($hash_omniscient2, $hash_mRNAGeneLink2) = slurp_gff3_file_JD({ input => $next_file,
	                                                                   config => $config
                                                                });
  print ("$next_file GFF3 file parsed\n");
  info_omniscient($hash_omniscient2);

  # Quick stat hash before complement
  my %quick_stat1;
  foreach my $level ( ('level1', 'level2') ){
    foreach  my $tag (keys %{$hash_omniscient->{$level}}) {
      my $nb_tag = keys %{$hash_omniscient->{$level}{$tag}};
      $quick_stat1{$level}{$tag} = $nb_tag;
    }
  }
  
  ####### COMPLEMENT #######
  complement_omniscients($hash_omniscient, $hash_omniscient2, $size_min); # deal with identical ID by renaming them
  print ("\nComplement done !\n");


 #RESUME COMPLEMENT
  my $complemented=undef;
  # Quick stat hash after complement
  my %quick_stat2;
  foreach my $level ( ('level1', 'level2') ){
    foreach  my $tag (keys %{$hash_omniscient->{$level}}) {
      my $nb_tag = keys %{$hash_omniscient->{$level}{$tag}};
      $quick_stat2{$level}{$tag} = $nb_tag;
    }
  }

  #About tag from hash1 added which exist in hash2
  foreach my $level ( ('level1', 'level2') ){
    foreach my $tag (keys %{$quick_stat1{$level}}){
      if ($quick_stat1{$level}{$tag} != $quick_stat2{$level}{$tag} ){
        print "We added ".($quick_stat2{$level}{$tag}-$quick_stat1{$level}{$tag})." $tag(s)\n";
        $complemented=1;
      }
    }
  }
  #About tag from hash2 added which dont exist in hash1
  foreach my $level ( ('level1', 'level2') ){
    foreach my $tag (keys %{$quick_stat2{$level}}){
      if (! exists $quick_stat1{$level}{$tag} ){
        print "We added ".$quick_stat2{$level}{$tag}." $tag(s)\n";
        $complemented=1;
      }
    }
  }
  #If nothing added
  if(! $complemented){
    print "\nNothing has been added\n";
  }
  else{
    print "\nNow the data contains:\n";
    info_omniscient($hash_omniscient);
  }
}

########
# Print results
print_omniscient( {omniscient => $hash_omniscient, output => $gffout} );
#END
print "usage: $0 @copyARGV\n";
my $end_run = time();
my $run_time = $end_run - $start_run;
print "Job done in $run_time seconds\n";
__END__

=head1 NAME

agat_sp_complement_annotations.pl

=head1 DESCRIPTION

The script allows to complement a reference annotation with other annotations.
A l1 feature from the addfile.gff that does not overlap a l1 feature from the reference annotation will be added.
A l1 feature from the addfile.gff without a CDS that overlaps a l1 feature with a CDS from the reference annotation will be added.
A l1 feature from the addfile.gff with a CDS that overlaps a l1 feature without a CDS from the reference annotation will be added.
A l1 feature from the addfile.gff with a CDS that overlaps a l1 feature with a CDS from the reference annotation will be added only if the CDSs don't overlap.
A l1 feature from the addfile.gff without a CDS that overlaps a l1 feature without a CDS from the reference annotation will be added only if none of the l3 features overlap.
/!\ It is sufficiant that only one isoform is overlapping to prevent the whole gene (l1 feature) from the addfile.gff to be added in the output.

=head1 SYNOPSIS

    agat_sp_complement_annotations.pl --ref annotation_ref.gff --add addfile1.gff --add addfile2.gff --out outFile
    agat_sp_complement_annotations.pl --help

=head1 OPTIONS

=over 8

=item B<--ref>,  B<-r> or B<-i>

Input GTF/GFF file used as reference.

=item B<--add> or B<-a>

Annotation(s) file you would like to use to complement the reference annotation. You can specify as much file you want like so: -a addfile1 -a addfile2 -a addfile3
/!\ The order you provide these files matter. Once the reference file has been complemented by file1, this new annotation becomes the new reference that will be complemented by file2 etc.
/!\ The result with -a addfile1 -a addfile2 will differ to the result from -a addfile2 -a addfile1. So, be aware of what you want if you use several addfiles.

=item  B<--size_min> or B<-s>

Option to keep the non-overlping gene only if the CDS size (in nucleotide) is over the minimum size defined. Default = 0 that means all of them are kept.

=item  B<--out>, B<--output>, B<--outfile> or B<-o>

Output gff3 containing the reference annotation with all the non-overlapping newly added genes from addfiles.gff.

=item B<-thread>, B<threads>, B<cpu>, B<cpus>, B<core>, B<cores>, B<job> or B<jobs>

Integer — Number of parallel processes to use for file input parsing (via forking).

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
