#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Clone 'clone';
use Getopt::Long;
use Pod::Usage;
use List::MoreUtils qw(uniq);
use Bio::Tools::GFF;
use Agat::Omniscient;

my $header = get_agat_header();
my $gff = undef;
my $opt_help= 0;
my $primaryTag=undef;
my $attributes=undef;
my $opt_merge = undef;
my $opt_comonTag=undef;
my $opt_output=undef;
my $add = undef;
my $cp = undef;

if ( !GetOptions(
    "help|h"          => \$opt_help,
    'c|ct=s'          => \$opt_comonTag,
    "gff|f=s"         => \$gff,
    'ml|merge_loci!'  => \$opt_merge,
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

if ( ! (defined($gff)) ){
    pod2usage( {
           -message => "$header\nAt least 1 parameter is mandatory:\nInput reference gff file (--gff) \n\n",
           -verbose => 0,
           -exitval => 2 } );
}

# Manage Output
my $ostream     = IO::File->new();
if(defined($opt_output))
{
$ostream->open( $opt_output, 'w' ) or
  croak(
     sprintf( "Can not open '%s' for reading: %s", $opt_output, $! ) );
}
else{
  $ostream->fdopen( fileno(STDOUT), 'w' ) or
      croak( sprintf( "Can not open STDOUT for writing: %s", $! ) );
}


                #####################
                #     MAIN          #
                #####################

my $cpt_tag=1;
my %tag_to_number;
my %number_to_tag;
my $content;
######################
### Parse GFF input #
my ($hash_omniscient, $hash_mRNAGeneLink) = slurp_gff3_file_JD({
                                                               input => $gff,
                                                               locus_tag => $opt_comonTag,
                                                               merge_loci => $opt_merge
                                                               });
print ("GFF3 file parsed\n");


foreach my $tag_l1 (keys %{$hash_omniscient->{'level1'}}){
  foreach my $id_l1 (keys %{$hash_omniscient->{'level1'}{$tag_l1}}){

    my $feature_l1=$hash_omniscient->{'level1'}{$tag_l1}{$id_l1};

    manage_attributes($feature_l1);

    #################
    # == LEVEL 2 == #
    #################
    foreach my $tag_l2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...

      if ( exists ($hash_omniscient->{'level2'}{$tag_l2}{$id_l1} ) ){
        foreach my $feature_l2 ( @{$hash_omniscient->{'level2'}{$tag_l2}{$id_l1}}) {

          manage_attributes($feature_l2);
          #################
          # == LEVEL 3 == #
          #################
          my $level2_ID = lc($feature_l2->_tag_value('ID'));

          foreach my $tag_l3 (keys %{$hash_omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...
            if ( exists ($hash_omniscient->{'level3'}{$tag_l3}{$level2_ID} ) ){
              foreach my $feature_l3 ( @{$hash_omniscient->{'level3'}{$tag_l3}{$level2_ID}}) {
                manage_attributes($feature_l3);
              }
            }
          }
        }
      }
    }
  }
}


print $ostream "seq_id\tsource_tag\tprimary_tag\tstart\tend\tscore\tstrand\tframe";
foreach my $key (sort { $a <=> $b } keys %number_to_tag){
   print $ostream "\t".$number_to_tag{$key};
}
print $ostream "\n".$content;
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

sub  manage_attributes{
  my  ($feature)=@_;
  $content .= $feature->seq_id."\t".$feature->source_tag."\t".$feature->primary_tag."\t".$feature->start."\t".$feature->end."\t".$feature->score."\t".$feature->strand."\t".$feature->frame;
  my @tag_list = $feature->get_all_tags();
  my %tag_hash;
  foreach my $tag (@tag_list) {
    $tag_hash{$tag}++;
  }

  foreach my $key (sort { $a <=> $b } keys %number_to_tag){
    my $sorted_tag = $number_to_tag{$key};
    if( exists_keys(\%tag_hash,($sorted_tag))){
      my @values = $feature->get_tag_values($sorted_tag);
      $content .= "\t".join(", ", @values);
      delete $tag_hash{$sorted_tag};
    }
  }

  foreach my $tag ( keys %tag_hash ) {
    $tag_to_number{$tag} = $cpt_tag;
    $number_to_tag{$cpt_tag} = $tag;
    my @values = $feature->get_tag_values($tag);
    $content .= "\t".join(", ", @values);
    $cpt_tag++;
  }
  $content .=  "\n";
}




__END__

=head1 NAME

agat_sp_to_tabulated.pl

=head1 DESCRIPTION

The script aims to convert gtf/gff file into tabulated file.
Attribute's tags from the 9th column become column titles.

=head1 SYNOPSIS

    agat_sp_to_tabulated.pl -gff file.gff [ -o outfile ]
    agat_sp_to_tabulated.pl --help

=head1 OPTIONS

=over 8

=item B<--gff> or B<-f>

Input GTF/GFF file.

=item B<-c> or B<--ct>

When the features doesn't have Parent/ID relationships, the parser will try to group
features using a common/shared attribute (i.e. a locus tag.). By default locus_tag and gene_id.
You can provide another specific common/shared attribute using this option.

=item B<--ml> or B<--merge_loci>

Merge loci parameter, default deactivated. You turn on the parameter if you want to merge loci into one locus when they overlap.
(at CDS level for mRNA, at exon level for other level2 features. Strand has to be the same). Prokaryote can have overlaping loci so it should not use it for prokaryote annotation.
In eukaryote, loci rarely overlap. Overlaps could be due to error in the file, mRNA can be merged under the same parent gene if you acticate the option.

=item B<-o> , B<--output> , B<--out> or B<--outfile>

Output GFF file.  If no output file is specified, the output will be
written to STDOUT.

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
