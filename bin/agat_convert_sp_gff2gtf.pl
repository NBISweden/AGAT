#!/usr/bin/env perl

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use Sort::Naturally;
use Bio::Tools::GFF;
use AGAT::Omniscient;

# for skipping data that may be represented elsewhere; currently, this is
# only the score
my %SKIPPED_TAGS = map { $_ => 1 } qw(score); # BIOPERL FIX.

my $header = get_agat_header();
my $outfile = undef;
my $gff = undef;
my $gtf_version = 3;
my $verbose = undef;
my $relax = undef;
my $help;
my @GTF3 = ("gene", "transcript", "exon", "CDS", "Selenocysteine", "start_codon", "stop_codon", "three_prime_utr", "five_prime_utr");
my @GTF2_5 = ("gene", "transcript", "exon", "CDS", "UTR", "start_codon", "stop_codon", "Selenocysteine");
my @GTF2_2 = ("CDS", "start_codon", "stop_codon", "5UTR", "3UTR", "inter", "inter_CNS", "intron_CNS", "exon");
my @GTF2_1 = ("CDS", "start_codon", "stop_codon", "exon", "5UTR", "3UTR");
my @GTF2 = ("CDS", "start_codon", "stop_codon", "exon");
my @GTF1 = ("CDS", "start_codon", "stop_codon", "exon", "intron");
my %GTF_LIST = (
	1 => \@GTF3,
	2 => \@GTF2,
	"2.1" => \@GTF2_1,
	"2.2" => \@GTF2_2,
	"2.5" => \@GTF2_5,
	3 => \@GTF3
);


if( !GetOptions(
    "help" => \$help,
    "gff|in=s" => \$gff,
		"gtf_version=s" => \$gtf_version,
		"verbose|v!" => \$verbose,
		"relax!" => \$relax,
    "outfile|output|o|out|gtf=s" => \$outfile))
{
    pod2usage( { -message => "Failed to parse command line.",
                 -verbose => 1,
                 -exitval => 1 } );
}

# Print Help and exit
if ($help) {
    pod2usage( { -message => "$header\n",
                 -verbose => 99,
                 -exitval => 0 } );
}

if ( ! (defined($gff)) ){
    pod2usage( {
           -message => "$header\nAt least 1 parameters is mandatory:\nInput gff file (--gff).\n\n",
           -verbose => 0,
           -exitval => 1 } );
}

## Manage output file
my $gtf_out;
if ($outfile) {
open(my $fh, '>', $outfile) or die "Could not open file '$outfile' $!";
  $gtf_out = Bio::Tools::GFF->new(-fh => $fh, -gff_version => 2.5);
}
else{
  $gtf_out = Bio::Tools::GFF->new(-fh => \*STDOUT, -gff_version => 2.5);
}


# check GTF versions
my @gtf_version_list = (1, 2, 2.1, 2.2, 2.5, 3);
my %gtf_version_hash = map { $_ => 1 } @gtf_version_list;
if(! exists_keys (\%gtf_version_hash, ("$gtf_version") ) ) {
	print "$gtf_version is not a valid GTF version. Please choose one among this list: @gtf_version_list\n"; exit;
}
print "converting to GTF$gtf_version\n";

######################
### Parse GFF input #
### Read gff input file.
my ($hash_omniscient, $hash_mRNAGeneLink) = slurp_gff3_file_JD({ input => $gff
                                                              });

# rebuild gene_id and transcript_id feature;

# sort by seq id
my $hash_sortBySeq = gather_and_sort_l1_by_seq_id($hash_omniscient);

my $gene_id=undef;
#################
# == LEVEL 1 == #
#################
foreach my $seqid (sort { (($a =~ /(\d+)$/)[0] || 0) <=> (($b =~ /(\d+)$/)[0] || 0) } keys %{$hash_sortBySeq}){ # loop over all the feature level1


	foreach my $primary_tag_key_level1 (sort {$a cmp $b} keys %{$hash_omniscient->{'level1'}}){
		foreach my $feature_level1 ( @{$hash_sortBySeq->{$seqid}{$primary_tag_key_level1}} ){
			my $id_tag_key_level1 = lc($feature_level1->_tag_value('ID'));

      # Gene ID level1
      my $gene_id_att=undef;
      if($feature_level1->has_tag('gene_id')){
        $gene_id_att=$feature_level1->_tag_value('gene_id');
      }

      my $transcript_id=undef;
      my $level3_gene_id=undef;
      #################
      # == LEVEL 2 == #
      #################
      foreach my $primary_tag_key_level2 ( sort {$a cmp $b} keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...

        if ( exists_keys ($hash_omniscient, ('level2', $primary_tag_key_level2, $id_tag_key_level1) ) ){
          foreach my $feature_level2 ( sort {$a->start <=> $b->start} @{$hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1}}) {



            # Gene ID level2
            my $gene_id_mrna_att=undef;
            if($feature_level2->has_tag('gene_id')){
              $gene_id_mrna_att=$feature_level2->_tag_value('gene_id');
            }

            my $transcript_id_mrna_att=undef;
            if($feature_level2->has_tag('transcript_id')){
              $transcript_id_mrna_att=$feature_level2->_tag_value('transcript_id');
            }

            # get gff3 feature (ID)
            my $level2_ID = lc($feature_level2->_tag_value('ID'));

            my $level3_transcript_id=undef;
            #################
            # == LEVEL 3 == #
            #################

            ############
            # Go through one time to check if gene_id and transcript_id are present and save them
            foreach my $primary_tag_key_level3 ( sort {$a cmp $b} keys %{$hash_omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...
              if ( exists_keys ($hash_omniscient, ('level3', $primary_tag_key_level3, $level2_ID) ) ){
                foreach my $feature_level3 ( sort { $a->start <=> $b->start } @{$hash_omniscient->{'level3'}{$primary_tag_key_level3}{$level2_ID}}) {

                  #Get level3 gene_id
                  if(! $level3_gene_id){
                    if($feature_level3->has_tag('gene_id')){
                      $level3_gene_id=$feature_level3->_tag_value('gene_id');
                    }
                  }

                  #Get level3 transcript_id
                  if(! $level3_transcript_id){
                    if($feature_level3->has_tag('transcript_id')){
                      $level3_transcript_id=$feature_level3->_tag_value('transcript_id');
                    }
                  }
                  if($level3_gene_id and $level3_transcript_id){last;}
                }
              }
              if($level3_gene_id and $level3_transcript_id){last;}
            }

            #################
            # CHOOSE the gene_id. We take the first from level1 to level3.
            if($gene_id_att){
              $gene_id=$gene_id_att;
            }
            elsif($gene_id_mrna_att){
              $gene_id=$gene_id_mrna_att
            }
            elsif($level3_gene_id){
              $gene_id=$level3_gene_id;
            }
            else{ # We didn't find any gene_id we will the ID of level1 as gene_id.
              $gene_id=$feature_level1->_tag_value('ID');
            }

            #################
            # CHOOSE the transcript_id. We take the first from level2 to level3.
            if($transcript_id_mrna_att){
              $transcript_id=$transcript_id_mrna_att;
            }
            elsif($level3_transcript_id){
              $transcript_id=$level3_transcript_id;
            }
            else{ # We didn't find any gene_id we will the ID of level2 as transcript_id.
              $transcript_id=$feature_level2->_tag_value('ID');
            }

            ##############
            # Second pass of level3 features
            # Add gene_id and transcript_id to level3 feature that don't have this information
            foreach my $primary_tag_key_level3 (keys %{$hash_omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...
              if ( exists_keys ($hash_omniscient, ('level3', $primary_tag_key_level3, $level2_ID) ) ){
                foreach my $feature_level3 ( @{$hash_omniscient->{'level3'}{$primary_tag_key_level3}{$level2_ID}}) {

                  #Check add gene_id
                  if(! $feature_level3->has_tag('gene_id')) {
                    $feature_level3->add_tag_value('gene_id', $gene_id);
                  }
                  elsif($feature_level3->_tag_value('gene_id') ne $gene_id) { #gene_id different, we replace it.
                    warn("We replace the transcript_id ".$feature_level3->_tag_value('gene_id')." by ".$gene_id.". Is it normal ?\n");exit;
                    $feature_level3->add_tag_value('gene_id', $gene_id);
                  }
                  #Check add transcript_id
                  if(! $feature_level3->has_tag('transcript_id')){
                    $feature_level3->add_tag_value('transcript_id', $transcript_id);
                  }
                  elsif($feature_level3->_tag_value('transcript_id') ne $transcript_id){ #transcript_id different, we replace it.
                    warn("We replace the transcript_id ".$feature_level3->_tag_value('transcript_id')." by ".$transcript_id.". Is it normal ?\n");exit;
                    $feature_level3->add_tag_value('transcript_id', $transcript_id);
                  }
                }
              }
            }

            ## add level2 missing information gene_id
            if(! $feature_level2->has_tag('gene_id')) {
               $feature_level2->add_tag_value('gene_id', $gene_id);
            }
            elsif($feature_level2->_tag_value('gene_id') ne $gene_id) { #gene_id different, we replace it.
              warn("We replace the transcript_id ".$feature_level2->_tag_value('gene_id')." by ".$gene_id.". Is it normal ?\n");exit;
              $feature_level2->add_tag_value('gene_id', $gene_id);
            }
            # add level2 missing information transcript_id
            if(! $feature_level2->has_tag('transcript_id')){
              $feature_level2->add_tag_value('transcript_id', $transcript_id);
            }
            elsif($feature_level2->_tag_value('transcript_id') ne $transcript_id){ #gene_id transcript_id, we replace it.
              warn("We replace the transcript_id ".$feature_level2->_tag_value('transcript_id')." by ".$transcript_id.". Is it normal ?\n");exit;
              $feature_level2->add_tag_value('transcript_id', $transcript_id);
            }
          }
        }
      }
      ## add level1 missing information gene_id
      if(! $feature_level1->has_tag('gene_id')) {
        $feature_level1->add_tag_value('gene_id', $gene_id);
      }
      elsif($feature_level1->_tag_value('gene_id') ne $gene_id) { #gene_id different, we replace it.
        warn("We replace the transcript_id ".$feature_level1->_tag_value('gene_id')." by ".$gene_id.". Is it normal ?\n");exit;
        $feature_level1->add_tag_value('gene_id', $gene_id);
      }
    }
  }
}

if (! $relax){
	# convert correct
	convert_feature_type($hash_omniscient, $gtf_version);

}

# print results
print_omniscient_filter($hash_omniscient, $gtf_version, $relax, $gtf_out, $relax);

print "Bye Bye\n";

# ---------------------------------

# convert feature type to correct one expected.
# All l1 will become gene type excepted for topfeature and standalone features
# that will be discarded.
sub convert_feature_type{
	my ($hash_omniscient, $gtf_version)=@_;

	my $topfeatures = get_feature_type_by_agat_value($hash_omniscient, 'level1', 'topfeature');
	my $standalones = get_feature_type_by_agat_value($hash_omniscient, 'level1', 'standalone');

	# all l1 are gene now
	foreach my $tag_l1 ( keys %{$hash_omniscient->{'level1'}}){
		if(exists_keys($topfeatures,($tag_l1))){ print "throw $tag_l1\n" if $verbose; next; }
		if(exists_keys($standalones,($tag_l1))){ print "throw $tag_l1\n" if $verbose; next; }

		foreach my $id_l1 ( keys %{$hash_omniscient->{'level1'}{$tag_l1}}){
			if (lc($tag_l1) ne "gene"){
				print "convert $tag_l1 to gene feature\n" if $verbose;
				my $feature = $hash_omniscient->{'level1'}{$tag_l1}{$id_l1};
				create_or_replace_tag( $feature, "original_biotype", $tag_l1 );
				$feature->primary_tag("gene");
			}

			# all l2 are transcript now
			foreach my $tag_l2 ( keys %{$hash_omniscient->{'level2'}}){
				if ( exists ($hash_omniscient->{'level2'}{$tag_l2}{$id_l1} ) ){
					foreach my $feature_level2 ( @{$hash_omniscient->{'level2'}{$tag_l2}{$id_l1}}) {
						if (lc($tag_l2) ne "transcript"){

							create_or_replace_tag($feature_level2, "original_biotype", $tag_l2 );
							$feature_level2->primary_tag("transcript");
						}

						# We do not check for GTF1 and GTF2 UTR specifically because we will remove them later
						if($gtf_version == "1" or $gtf_version == "2"){next;}

						# get gff3 feature (ID)
						my $level2_ID = lc($feature_level2->_tag_value('ID'));

						foreach my $tag_l3_lc (keys %{$hash_omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...

							if ( $tag_l3_lc =~ "utr"){
								if ( exists_keys ($hash_omniscient, ('level3', $tag_l3_lc, $level2_ID) ) ){
									foreach my $feature_level3 ( @{$hash_omniscient->{'level3'}{$tag_l3_lc}{$level2_ID}}) {
										my $tag_l3 = $feature_level3->primary_tag;
										if ($gtf_version eq "2.5"){
											if ($tag_l3 ne "UTR"){
												create_or_replace_tag($feature_level3, "original_biotype", $tag_l3 );
												$feature_level3->primary_tag("UTR");
											}
										}
										else{
											if ( ($tag_l3 =~ "5") or ( $tag_l3_lc =~ "five") ){
												if ($gtf_version eq "3"){
													if ($tag_l3 ne "five_prime_utr"){
														create_or_replace_tag($feature_level3, "original_biotype", $tag_l3 );
														$feature_level3->primary_tag("five_prime_utr");
													}
												}
												else{
													if ($tag_l3 ne "5UTR"){
														create_or_replace_tag($feature_level3, "original_biotype", $tag_l3 );
														$feature_level3->primary_tag("5UTR");
													}
												}
											}
											if ( ($tag_l3 =~ "3") or ( $tag_l3_lc =~ "three") ){
												if ($gtf_version eq "3"){
													if ($tag_l3 ne "three_prime_utr"){
														create_or_replace_tag($feature_level3, "original_biotype", $tag_l3 );
														$feature_level3->primary_tag("three_prime_utr");
													}
												}
												else{
													if ($tag_l3 ne "3UTR"){
														create_or_replace_tag($feature_level3, "original_biotype", $tag_l3 );
														$feature_level3->primary_tag("3UTR");
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
			}
		}
	}
}


# filter feature type to remove not expected ones
sub print_omniscient_filter{
	my ($hash_omniscient, $gtf_version, $relax, $gffout)=@_;

	# --------- deal GTF version --------------
	my $list_ok = $GTF_LIST{$gtf_version};
	$_=lc for @$list_ok;
	my %hash_ok = map { $_ => 1 } @$list_ok;

	# --------- deal with header --------------
  write_headers_gtf($hash_omniscient, $gffout, $gtf_version, $relax);

	# sort by seq id
	my ( $hash_sortBySeq, $hash_sortBySeq_stdf,  $hash_sortBySeq_topf) = collect_l1_info_sorted_by_seqid_and_location($hash_omniscient);

	# Read by seqId to sort properly the output by seq ID
	# sort { (($a =~ /(\d+)$/)[0] || 0) <=> (($b =~ /(\d+)$/)[0] || 0) will provide sorting like that: contig contig1 contig2 contig3 contig10 contig11 contig22 contig100 contig101
	foreach my $seqid (sort { (($a =~ /(\d+)$/)[0] || 0) <=> (($b =~ /(\d+)$/)[0] || 0) } keys %{$hash_sortBySeq}){ # loop over all the feature level1

		# ----- LEVEL 1 -----
		write_top_features_gtf($gffout, $seqid, $hash_sortBySeq_topf, $hash_omniscient, \%hash_ok, $relax);

		foreach my $locationid ( sort { ncmp ($a, $b) } keys %{$hash_sortBySeq->{$seqid} } ){
			my $primary_tag_l1 = $hash_sortBySeq->{$seqid}{$locationid}{'tag'};
			my $id_tag_key_level1 = $hash_sortBySeq->{$seqid}{$locationid}{'id'};

			my $feature_l1 = $hash_omniscient->{'level1'}{$primary_tag_l1}{$id_tag_key_level1};
			my $primary_tag_l1_gtf = lc($feature_l1->primary_tag());

			if(exists_keys (\%hash_ok, ( $primary_tag_l1_gtf) ) or $relax ) {
				write_feature_JD(	$gffout, $feature_l1 ); # print feature
			}

			# ----- LEVEL 2 -----
			foreach my $primary_tag_l2 (sort {$a cmp $b} keys %{$hash_omniscient->{'level2'}}){ # primary_tag_l2 = mrna or mirna or ncrna or trna etc...
				if ( exists_keys( $hash_omniscient, ('level2', $primary_tag_l2, $id_tag_key_level1) ) ){
					foreach my $feature_level2 ( sort { ncmp ($a->start.$a->end.$a->_tag_value('ID'), $b->start.$b->end.$b->_tag_value('ID') ) } @{$hash_omniscient->{'level2'}{$primary_tag_l2}{$id_tag_key_level1}}) {

						my $primary_tag_l2_gtf = lc($feature_level2->primary_tag());
						if(exists_keys (\%hash_ok, ($primary_tag_l2_gtf) )  or $relax ) {
							write_feature_JD($gffout, $feature_level2); # print feature
						}

						# ----- LEVEL 3 -----
						my $level2_ID = lc($feature_level2->_tag_value('ID'));
						my @l3_done;

						######
						# FIRST EXON
						if ( exists_keys( $hash_omniscient, ('level3', 'exon', $level2_ID) ) ){
							foreach my $feature_level3 ( sort {$a->start <=> $b->start} @{$hash_omniscient->{'level3'}{'exon'}{$level2_ID}}) {
								write_feature_JD($gffout, $feature_level3);
								push @l3_done, 'exon';
							}
						}
						###########
						# SECOND CDS
						if ( exists_keys( $hash_omniscient, ('level3', 'cds', $level2_ID) ) ){
							foreach my $feature_level3 ( sort {$a->start <=> $b->start} @{$hash_omniscient->{'level3'}{'cds'}{$level2_ID}}) {
								write_feature_JD($gffout, $feature_level3);
								push @l3_done, 'cds';
							}
						}

						############
						# THEN ALL THE REST
						foreach my $primary_tag_l3 (sort {$a cmp $b} keys %{$hash_omniscient->{'level3'}}){ # primary_tag_l3 = cds or exon or start_codon or utr etc...
							if (! grep { $_ eq $primary_tag_l3 } @l3_done){
								if ( exists_keys( $hash_omniscient, ('level3', $primary_tag_l3, $level2_ID) ) ){
									my $primary_tag_l3_gtf = lc($hash_omniscient->{'level3'}{$primary_tag_l3}{$level2_ID}->[0]->primary_tag() );

									if(exists_keys (\%hash_ok, ( $primary_tag_l3_gtf) )  or $relax ) {
										foreach my $feature_level3 ( sort {$a->start <=> $b->start} @{$hash_omniscient->{'level3'}{$primary_tag_l3}{$level2_ID}}) {
											write_feature_JD($gffout, $feature_level3);
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
}

# @Purpose: Print the headers when first time we access the fh
# @input: 2 =>  ref omniscient, gff fh
# @output none => none
sub write_headers_gtf{
  	my ($hash_omniscient, $gffout, $gtf_version, $relax ) = @_  ;

  # If first time we write in this fh
  if ($gffout->{'_first'} ){

    # we write the very top header describing gff version
		if ($relax){
			$gffout->_print("##gtf-version X\n");
			my $gffXtra=$gffout->{"_filehandle"}; #to add extra lines to gff!!
			print $gffXtra "# GFF-like GTF i.e. not checked against any GTF specification. Conversion based on GFF input, standardised by AGAT.";
		}else{
    	$gffout->_print("##gtf-version ".$gtf_version."\n");
		}

    # Now we inject the header catched when parsing input file
    if (exists_keys( $hash_omniscient, ('other', 'header') ) ){
      my $gffXtra=$gffout->{"_filehandle"}; #to add extra lines to gff!!
      foreach my $header_line ( @{$hash_omniscient->{'other'}{'header'} } ) {
        print $gffXtra $header_line;
      }
    }
    # to avoid to write again the header later in the file
    $gffout->{'_first'} = 0;
  }
}

sub write_top_features_gtf{

  my ( $gffout, $seqid, $hash_sortBySeq_topf, $hash_omniscient, $hash_ok, $relax ) = @_;

  if ( exists_keys( $hash_sortBySeq_topf, ($seqid) ) ){

    foreach my $locationid ( sort { ncmp ($a, $b) } keys %{$hash_sortBySeq_topf->{$seqid}} ){
			my $tag_l1 = $hash_sortBySeq_topf->{$seqid}{$locationid}{'tag'};
			my $id_l1 = $hash_sortBySeq_topf->{$seqid}{$locationid}{'id'};
			my $feature_l1 = $hash_omniscient->{'level1'}{$tag_l1}{$id_l1};

			if(exists_keys ($hash_ok, ( $tag_l1) ) or $relax ) {
				write_feature_JD($gffout, $feature_l1); # print feature
			}
    }
  }
}

sub write_feature_JD {
    my ($self, @features) = @_;
    return unless @features;
    if( $self->{'_first'} && $self->gff_version() == 3 ) {
        $self->_print("##gff-version 3\n");
    }
    $self->{'_first'} = 0;
    foreach my $feature ( @features ) {
        $self->_print(_gff25_string_JD($self, $feature)."\n");
    }
}

=head2 _gff25_string_JD

 Title   : _gff25_string
 Usage   : $gffstr = $gffio->_gff2_string
 Function: To get a format of GFF that is peculiar to Gbrowse/Bio::DB::GFF
 Example : 9th column: ID "gene-1"; Name "name 1" name2;
 Returns : A GFF2.5-formatted string representation of the SeqFeature
 Args    : A Bio::SeqFeatureI implementing object to be GFF2.5-stringified
 Comments: GFF2.5 is suposed to be similar as GTF (with semicolon at the end).
=cut

sub _gff25_string_JD {
    my ($gff, $origfeat) = @_;
    my $feat;
    if ($origfeat->isa('Bio::SeqFeature::FeaturePair')){
        $feat = $origfeat->feature2;
    } else {
        $feat = $origfeat;
    }
    my ($str1, $str2,$score,$frame,$name,$strand);

    if( $feat->can('score') ) {
        $score = $feat->score();
    }
    $score = '.' unless defined $score;

    if( $feat->can('frame') ) {
        $frame = $feat->frame();
    }
    $frame = '.' unless defined $frame;

    $strand = $feat->strand();
    if(! $strand) {
        $strand = ".";
    } elsif( $strand == 1 ) {
        $strand = '+';
    } elsif ( $feat->strand == -1 ) {
        $strand = '-';
    }

    if( $feat->can('seqname') ) {
        $name = $feat->seq_id();
    }
    $name ||= 'SEQ';

    $str1 = join("\t",
                 $name,
                 $feat->source_tag(),
                 $feat->primary_tag(),
                 $feat->start(),
                 $feat->end(),
                 $score,
                 $strand,
                 $frame);

    my @all_tags = $feat->all_tags;
    my @group; my @firstgroup;

    if (@all_tags) {   # only play this game if it is worth playing...
        foreach my $tag ( @all_tags ) {
            next if exists $SKIPPED_TAGS{$tag};
            my @v;
            foreach my $value ( $feat->get_tag_values($tag) ) {
                unless( defined $value && length($value) ) {
                    $value = '""';
                } else{ # quote all type of values
                    $value =~ s/\t/\\t/g; # substitute tab and newline
                    # characters
                    $value =~ s/\n/\\n/g; # to their UNIX equivalents
                    $value = '"' . $value . '"';
                }
                push @v, $value;
            }
            $v[$#v] =~ s/\s+$//; #remove left space of the last value
            if (($tag eq 'gene_id') || ($tag eq 'transcript_id')){ # hopefully we won't get both...
                push @firstgroup, "$tag ".join(" ", @v);
            } else {
                push @group, "$tag ".join(" ", @v);
            }
        }
    }
		@firstgroup = sort @firstgroup if @firstgroup;
    $str2 = join('; ', (@firstgroup, @group));
    $str2 = $str2.";";
    # Add Target information for Feature Pairs
    if( ! $feat->has_tag('Target') && # This is a bad hack IMHO
        ! $feat->has_tag('Group') &&
        $origfeat->isa('Bio::SeqFeature::FeaturePair') ) {
        $str2 = sprintf("Target %s ; tstart %d ; tend %d", $origfeat->feature1->seq_id,
                        ( $origfeat->feature1->strand < 0 ?
                          ( $origfeat->feature1->end,
                            $origfeat->feature1->start) :
                          ( $origfeat->feature1->start,
                            $origfeat->feature1->end)
                        )) . ($str2?" ; ".$str2:""); # need to put the target info before other tag/value pairs - mw
    }
    return $str1 . "\t".  $str2;
}

__END__

=head1 NAME

agat_convert_sp_gff2gtf.pl

=head1 DESCRIPTION

The script aims to convert any GTF/GFF file into a proper GTF file.
Full information about the format can be found here: https://github.com/NBISweden/GAAS/blob/master/annotation/knowledge/gxf.md
You can choose among 6 different GTF types (1, 2, 2.1, 2.2, 2.5, 3).
Depending the version selected the script will filter out the features that are not accepted.
For GTF2.5 and 3, every level1 feature (e.g nc_gene pseudogene) will be converted into
gene feature and every level2 feature (e.g mRNA ncRNA) will be converted into
transcript feature.
You can even produce a GFF-like GTF using the --relax option. It allows to keep all
original feature types (3rd column).

To be fully GTF compliant all feature have a gene_id and a transcript_id attribute.
The gene_id	is unique identifier for the genomic source of the transcript, which is
used to group transcripts into genes.
The transcript_id	is a unique identifier for the predicted transcript,
which is used to group features into transcripts.

=head1 SYNOPSIS

    agat_convert_sp_gff2gtf.pl --gff infile.gtf [ -o outfile ]
    agat_convert_sp_gff2gtf -h

=head1 OPTIONS

=over 8

=item B<--gff> or B<--in>

Input GFF file that will be read

=item B<--gtf_version>
version of the GTF output. Default 3 (for GTF3)

GTF3 (9 feature types accepted): gene, transcript, exon, CDS, Selenocysteine, start_codon, stop_codon, three_prime_utr and five_prime_utr

GTF2.5 (8 feature types accepted): gene, transcript, exon, CDS, UTR, start_codon, stop_codon, Selenocysteine

GTF2.2 (9 feature types accepted): CDS, start_codon, stop_codon, 5UTR, 3UTR, inter, inter_CNS, intron_CNS and exon

GTF2.1 (6 feature types accepted): CDS, start_codon, stop_codon, exon, 5UTR, 3UTR

GTF2 (4 feature types accepted): CDS, start_codon, stop_codon, exon

GTF1 (5 feature types accepted): 	CDS, start_codon, stop_codon, exon, intron

=item B<--relax>

Relax option avoid to apply strict GTF format specification. All feature type will be kept.
No modification e.g. mRNA to transcript.
No filtering i.e. feature type not accepted by GTF format are kept.
gene_id and transcript_id attributes will be added, and the attributes will follow the
GTF formating.

=item B<-o> , B<--output> , B<--out> , B<--outfile> or B<--gtf>

Output GTF file. If no output file is specified, the output will be
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
