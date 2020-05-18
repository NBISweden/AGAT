#!/usr/bin/env perl

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use AGAT::Omniscient;
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

######################
# Manage output file #

if ($opt_output) {
  $opt_output=~ s/.gff//g;
  }
else{
  print "Default output name: split_result\n";
  $opt_output="split_result";
}

if (-d $opt_output){
  print "The output directory choosen already exists. Please give me another Name.\n";exit();
}
mkdir $opt_output;

                #####################
                #     MAIN          #
                #####################

######################
### Parse GFF input #
my ($hash_omniscient, $hash_mRNAGeneLink) = slurp_gff3_file_JD({ input => $opt_gfffile
                                                              });
print ("GFF3 file parsed\n");


my $topfeatures = get_feature_type_by_agat_value($hash_omniscient, 'level1', 'topfeature');
my $standalones = get_feature_type_by_agat_value($hash_omniscient, 'level1', 'standalone');
my %handlers;
my $gffout;
#################
# == LEVEL 1 == #
foreach my $tag_l1 (keys %{$hash_omniscient->{'level1'}}){ # primary_tag_key_level1 = gene or repeat etc...
	# deal with topfeatures and standalone feature type
	if ( exists_keys ($topfeatures, ($tag_l1) ) or  exists_keys ($standalones, ($tag_l1) ) ){
		open(my $fh, '>', $opt_output."/".$tag_l1.".gff") or die "Could not open file '$tag_l1' $!";
		$gffout= Bio::Tools::GFF->new(-fh => $fh, -gff_version => 3 );
		foreach my $key_l1 (keys %{$hash_omniscient->{'level1'}{$tag_l1}}){
			$gffout->write_feature($hash_omniscient->{'level1'}{$tag_l1}{$key_l1});
		}
		$gffout->close();
	}
	# deal with everything that is not topfeatures or standalone feature
	else{
		foreach my $key_l1 (keys %{$hash_omniscient->{'level1'}{$tag_l1}}){

	    #################
	    # == LEVEL 2 == #
	    my $level1_printed=undef;
	    foreach my $tag_l2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...

	      if ( exists_keys ($hash_omniscient, ('level2', $tag_l2, $key_l1) ) ){
	        foreach my $feature_level2 ( @{$hash_omniscient->{'level2'}{$tag_l2}{$key_l1}}) {
	          #manage handler
	          if(! exists_keys ( \%handlers, ($tag_l2) ) ) {
	            open(my $fh, '>', $opt_output."/".$tag_l2.".gff") or die "Could not open file '$tag_l2' $!";
	            $gffout= Bio::Tools::GFF->new(-fh => $fh, -gff_version => 3 );
	            $handlers{$tag_l2}=$gffout;
	          }
	          $gffout = $handlers{$tag_l2};

	          #################
	          # == LEVEL 1 == #
	          if(! $level1_printed){
	            $gffout->write_feature($hash_omniscient->{'level1'}{$tag_l1}{$key_l1}); # print feature
	            $level1_printed=1;
	          }

	          #################
	          # == LEVEL 2 == #
	          $gffout->write_feature($feature_level2);

	          #################
	          # == LEVEL 3 == #
	          my $level2_ID = lc($feature_level2->_tag_value('ID'));

	          ###########
	          # Before tss
	          if ( exists_keys($hash_omniscient,('level3','tss',$level2_ID)) ){
	            foreach my $feature_level3 ( @{$hash_omniscient->{'level3'}{'tss'}{$level2_ID}}) {
	              $gffout->write_feature($feature_level3);
	            }
	          }

	          ######
	          # FIRST EXON
	          if ( exists_keys($hash_omniscient,('level3','exon',$level2_ID)) ){
	            foreach my $feature_level3 ( @{$hash_omniscient->{'level3'}{'exon'}{$level2_ID}}) {
	              $gffout->write_feature($feature_level3);
	            }
	          }
	          ###########
	          # SECOND CDS
	          if ( exists_keys($hash_omniscient,('level3','cds',$level2_ID)) ){
	            foreach my $feature_level3 ( @{$hash_omniscient->{'level3'}{'cds'}{$level2_ID}}) {
	              $gffout->write_feature($feature_level3);
	            }
	          }

	          ###########
	          # Last tts
	          if ( exists_keys($hash_omniscient,('level3','tts',$level2_ID)) ){
	            foreach my $feature_level3 ( @{$hash_omniscient->{'level3'}{'tts'}{$level2_ID}}) {
	              $gffout->write_feature($feature_level3);
	            }
	          }

	          ###########
	          # The rest
	          foreach my $primary_tag_key_level3 (keys %{$hash_omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...
	            if( ($primary_tag_key_level3 ne 'cds') and ($primary_tag_key_level3 ne 'exon') and ($primary_tag_key_level3 ne 'tss') and ($primary_tag_key_level3 ne 'tts')){
	              if ( exists_keys ($hash_omniscient, ('level3', $primary_tag_key_level3, $level2_ID) ) ){
	                foreach my $feature_level3 ( @{$hash_omniscient->{'level3'}{$primary_tag_key_level3}{$level2_ID}}) {
	                  $gffout->write_feature($feature_level3);
	                }
	              }
	            }
	          }
	        }
	      }
	    }
	  }
	}
}

#Close all FH opened
foreach my $key (keys %handlers){
	$handlers{$key}->close();
}

my $end_run = time();
my $run_time = $end_run - $start_run;
print "Job done in $run_time seconds\n";
__END__

=head1 NAME

agat_sp_separate_by_record_type.pl

=head1 DESCRIPTION

The script will separate the features from the gff input file into different files according to
the record type. A record represent all features linked collectively by Parent/ID relationships.
(e.g gene + mrna + exon + cds + utr of a locus).

a) When the record contains Level2 feature, the record type is the Level2 feature type (e.g tRNA,mRNA,ncRNA etc...)
b) Some features do not have children (top and standalone level1 features) e.g. location,region,chromosome.
In such case the record type is the level1 feature type.

=head1 SYNOPSIS

    agat_sp_separate_by_record_type.pl -g infile.gff [ -o outfolder ]
    agat_sp_separate_by_record_type.pl --help

=head1 OPTIONS

=over 8

=item B<-g>, B<--gff> or B<-ref>

Input GTF/GFF file.

=item B<-o> or B<--output>

Output folder.  If no output folder provided, the default name will be <split_result>.

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
