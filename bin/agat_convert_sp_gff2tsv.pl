#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Clone 'clone';
use Sort::Naturally;
use Getopt::Long;
use Pod::Usage;
use List::MoreUtils qw(uniq);
use AGAT::AGAT;

my $header = get_agat_header();
my $config;
my $gff = undef;
my $opt_help= 0;
my $primaryTag=undef;
my $attributes=undef;
my $opt_output=undef;
my $add = undef;
my $cp = undef;

if ( !GetOptions(
    "h|help"          => \$opt_help,
    "gff|f=s"         => \$gff,
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
           -exitval => 1 } );
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

######################
### Parse GFF input #
my ($hash_omniscient, $hash_mRNAGeneLink) = slurp_gff3_file_JD({ input => $gff,
                                                                 config => $config });
print ("GFF3 file parsed\n");

# ---- List attributes ----
my $attribute_bucket = get_all_attributes($hash_omniscient);
my $nb_attributes= keys %{$attribute_bucket};

# ---- print header specifing columns ----
print $ostream "seq_id\tsource_tag\tprimary_tag\tstart\tend\tscore\tstrand\tframe";
foreach my $key ( sort { "\L$a" cmp "\L$b" } keys %{$attribute_bucket} ){
	 print $ostream "\t".$key;
}
print $ostream "\n";

# ---- Go through features and print tsv values ----
# sort by seq id
my ( $hash_sortBySeq, $hash_sortBySeq_stdf,  $hash_sortBySeq_topf) = collect_l1_info_sorted_by_seqid_and_location($hash_omniscient);

#Read by seqId to sort properly the output by seq ID
foreach my $seqid ( sort { (($a =~ /(\d+)$/)[0] || 0) <=> (($b =~ /(\d+)$/)[0] || 0) } keys %{$hash_sortBySeq}){ # loop over all the feature level1

	#################
	# == LEVEL 1 == #
	#################
	#write_top_features($gffout, $seqid, $hash_sortBySeq_topf, $hash_omniscient);

	foreach my $locationid ( sort { ncmp ($a, $b) } keys %{$hash_sortBySeq->{$seqid} } ){

		my $tag_l1 = $hash_sortBySeq->{$seqid}{$locationid}{'tag'};
		my $id_l1 = $hash_sortBySeq->{$seqid}{$locationid}{'id'};

    my $feature_l1=$hash_omniscient->{'level1'}{$tag_l1}{$id_l1};
    print_tsv_feature($feature_l1,$attribute_bucket);

    #################
    # == LEVEL 2 == #
    #################
    foreach my $tag_l2 (sort keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...

      if ( exists ($hash_omniscient->{'level2'}{$tag_l2}{$id_l1} ) ){
        foreach my $feature_l2 ( @{$hash_omniscient->{'level2'}{$tag_l2}{$id_l1}}) {

          print_tsv_feature($feature_l2,$attribute_bucket);
          #################
          # == LEVEL 3 == #
          #################
          my $level2_ID = lc($feature_l2->_tag_value('ID'));
					my @l3_done;

					if ( exists_keys($hash_omniscient,('level3','tss',$level2_ID)) ){
						foreach my $feature_l3 ( @{$hash_omniscient->{'level3'}{'tss'}{$level2_ID}}) {
							print_tsv_feature($feature_l3,$attribute_bucket);
							push @l3_done, 'tss';
						}
					}
					if ( exists_keys($hash_omniscient,('level3','exon',$level2_ID)) ){
						foreach my $feature_l3 ( @{$hash_omniscient->{'level3'}{'exon'}{$level2_ID}}) {
							print_tsv_feature($feature_l3,$attribute_bucket);
							push @l3_done, 'exon';
						}
					}
					if ( exists_keys($hash_omniscient,('level3','cds',$level2_ID)) ){
						foreach my $feature_l3 ( @{$hash_omniscient->{'level3'}{'cds'}{$level2_ID}}) {
							print_tsv_feature($feature_l3,$attribute_bucket);
							push @l3_done, 'cds';
						}
					}
					if ( exists_keys($hash_omniscient,('level3','tts',$level2_ID)) ){
						foreach my $feature_l3 ( @{$hash_omniscient->{'level3'}{'tts'}{$level2_ID}}) {
							print_tsv_feature($feature_l3,$attribute_bucket);
							push @l3_done, 'tts';
						}
					}

          foreach my $tag_l3 (sort keys %{$hash_omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...
						if (! grep { $_ eq $tag_l3 } @l3_done){
							if ( exists_keys( $hash_omniscient, ('level3', $tag_l3, $level2_ID) ) ){
	              foreach my $feature_l3 ( @{$hash_omniscient->{'level3'}{$tag_l3}{$level2_ID}}) {
	                print_tsv_feature($feature_l3,$attribute_bucket);
	              }
							}
            }
          }
        }
      }
    }
  }
}

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

sub get_all_attributes{
	my ($hash_omniscient)=@_;

	my $attribute_bucket;
	# == LEVEL 1 == #
	foreach my $tag_l1 ( keys %{$hash_omniscient->{'level1'}} ){ # primary_tag_l1 = gene or repeat etc...
		foreach my $id_l1 ( keys %{$hash_omniscient->{'level1'}{$tag_l1}} ) { #sort by position
			my $feature_l1 = $hash_omniscient->{'level1'}{$tag_l1}{$id_l1};

			my @tag_list = $feature_l1->get_all_tags();
			$attribute_bucket->{$_}++ for (@tag_list);
			# == LEVEL 2 == #
			foreach my $tag_l2 (sort keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...

				if ( exists_keys ($hash_omniscient, ('level2', $tag_l2, $id_l1) ) ){
					foreach my $feature_l2 ( @{$hash_omniscient->{'level2'}{$tag_l2}{$id_l1}}) {

						my @tag_list = $feature_l2->get_all_tags();
						$attribute_bucket->{$_}++ for (@tag_list);

						# == LEVEL 3 == #
						my $level2_ID = lc($feature_l2->_tag_value('ID'));
						foreach my $tag_l3 (sort keys %{$hash_omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...

							if ( exists_keys( $hash_omniscient, ('level3', $tag_l3, $level2_ID ) ) ){
								foreach my $feature_l3 ( @{$hash_omniscient->{'level3'}{$tag_l3}{$level2_ID} } ) {
									my @tag_list = $feature_l3->get_all_tags();
									$attribute_bucket->{$_}++ for (@tag_list);
								}
							}
						}
					}
				}
			}
		}
	}
	return $attribute_bucket;
}

sub print_tsv_feature{
  my  ($feature, $attribute_bucket)=@_;

	# get 8 first columns
  my $score = ($feature->score) ? $feature->score : ".";
  my $tsv_line = $feature->seq_id."\t".$feature->source_tag."\t".$feature->primary_tag."\t".$feature->start."\t".$feature->end."\t".$score."\t".$feature->strand."\t".$feature->frame;

	# get attributes of the 9th column
	my @tag_list = $feature->get_all_tags();
  my %tag_hash;
	$tag_hash{$_}++ for (@tag_list);

	# print in same order than header (N/A if attribute does not exist)
	foreach my $key_bucket ( sort { "\L$a" cmp "\L$b" } keys %{$attribute_bucket} ){
		if ( exists_keys(\%tag_hash, ($key_bucket) ) ){
			my @values = $feature->get_tag_values($key_bucket);
			$tsv_line .= "\t".join(", ", @values);
		}
		else{
			$tsv_line .= "\tN/A";# attribute do not exists in this feature. Set N/A
		}
	}

	$tsv_line .= "\n";

	print $ostream $tsv_line;
}




__END__

=head1 NAME

agat_convert_sp_gff2tsv.pl

=head1 DESCRIPTION

The script aims to convert gtf/gff file into tabulated file.
Attribute's tags from the 9th column become column titles.

=head1 SYNOPSIS

    agat_convert_sp_gff2tsv.pl -gff file.gff [ -o outfile ]
    agat_convert_sp_gff2tsv.pl --help

=head1 OPTIONS

=over 8

=item B<--gff> or B<-f>

Input GTF/GFF file.

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
