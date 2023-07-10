#!/usr/bin/env perl

use strict;
use warnings;
use Pod::Usage;
use Clone 'clone';
use Getopt::Long;
use Sort::Naturally;
use Bio::SeqIO;
use Bio::DB::Fasta;
use AGAT::AGAT;

my $header = get_agat_header();
my $config;
my $opt_plus_strand = undef;
my $start_run = time();
my $codonTable=1;

my $opt_gfffile;
my $opt_fastafile;
my $opt_output;
my $opt_AA=undef;
my $opt_help = 0;
my $opt_full=undef;
my $opt_split=undef;
my $opt_merge=undef;
my $opt_extremity_only=undef;
my $opt_upstreamRegion=undef;
my $opt_downRegion=undef;
my $opt_cdna=undef;
my $opt_mrna=undef;
my $opt_OFS=undef;
my $opt_type = 'cds';
my $opt_keep_attributes = undef;
my $opt_keep_parent_attributes = undef;
my $opt_cleanFinalStop=undef;
my $opt_cleanInternalStop=undef;
my $opt_remove_orf_offset = undef;
my $opt_quiet = undef;
my $opt_plus_strand_only = undef;
my $opt_revcomp=undef;
my $opt_verbose=undef;
my $opt_alternative_start_codon = undef;

# OPTION MANAGMENT
my @copyARGV=@ARGV;
if ( !GetOptions( 'alternative_start_codon|asc!' => \$opt_alternative_start_codon,
                  'cdna!'                        => \$opt_cdna,
                  'cfs|clean_final_stop!'        => \$opt_cleanFinalStop,
                  'cis|clean_internal_stop!'     => \$opt_cleanInternalStop,
                  'do|3|three|down|downstream=i' => \$opt_downRegion,
                  'eo!'                          => \$opt_extremity_only,
                  'f|fa|fasta=s'                 => \$opt_fastafile,
                  'full!'                        => \$opt_full,
                  'g|gff=s'                      => \$opt_gfffile,
                  'c|config=s'               => \$config,
                  'h|help!'                      => \$opt_help,
                  'merge!'                       => \$opt_merge,
                  'mrna|transcript!'             => \$opt_mrna,
                  'ofs=s'                        => \$opt_OFS,
                  'o|output=s'                   => \$opt_output,
                  'plus_strand_only!'            => \$opt_plus_strand_only,
                  'p|protein|aa!'                => \$opt_AA,
                  'q|quiet!'                     => \$opt_quiet,
                  'remove_orf_offset|roo!'       => \$opt_remove_orf_offset,
                  'revcomp!'                     => \$opt_revcomp,
                  'split!'                       => \$opt_split,
                  'keep_attributes!'             => \$opt_keep_attributes,
                  'keep_parent_attributes!'      => \$opt_keep_parent_attributes,
                  'table|codon|ct=i'             => \$codonTable,
                  't|type=s'                     => \$opt_type,
                  'up|5|five|upstream=i'         => \$opt_upstreamRegion,
                  'verbose|v!'                   => \$opt_verbose ) )
{
    pod2usage( { -message => "$header\nFailed to parse command line",
                 -verbose => 1,
                 -exitval => 1 } );
}


# Print Help and exit
if ($opt_help) {
    pod2usage( { -verbose => 99,
                 -exitval => 0,
                 -message => "$header\n" } );
}

if ( (! (defined($opt_gfffile)) ) or (! (defined($opt_fastafile)) ) ){
    pod2usage( {
           -message => "\nAt least 2 parametes are mandatory:\nInput reference gff file (-g);  Input reference fasta file (-f)\n\n".
           "Output is optional. Look at the help documentation to know more.\n",
           -verbose => 0,
           -exitval => 2 } );
}

# --- Manage config ---
$config = get_agat_config({config_file_in => $config});

# activate warnings limit
my %warnings;
activate_warning_limit(\%warnings, 10);

# shortcut for cdna
if($opt_cdna){$opt_type="exon"; $opt_merge=1; $opt_revcomp=1;}
# shortcut for mrna/transcript
if($opt_mrna){$opt_type="exon"; $opt_merge=1;}

if( $opt_full   and $opt_split){print "Options --full and --split cannot be used concomitantly.\n"; exit;}
if( $opt_full   and $opt_merge){print "Options --full and --merge cannot be used concomitantly.\n"; exit;}
if( $opt_split   and $opt_merge){print "Options --split and --merge cannot be used concomitantly.\n"; exit;}

my $ostream;
if ($opt_output) {
  open(my $fh, '>', $opt_output) or die "Could not open file '$opt_output' $!";
  $ostream= Bio::SeqIO->new(-fh => $fh, -format => 'Fasta' );
}
else{
  $ostream = Bio::SeqIO->new(-fh => \*STDOUT, -format => 'Fasta');
}

print "We will extract the $opt_type sequences.\n";
$opt_type=lc($opt_type);

# deal with codon table
if($codonTable<0 and $codonTable>25){
  print "$codonTable codon table is not a correct value. It should be between 0 and 25 (0,23 and 25 can be problematic !)\n";
}

# deal with OFS
my $OFS=" ";
if($opt_OFS){
  $OFS = $opt_OFS;
}

# if parental option activated with must be sure that the default one is activated too.
if ($opt_keep_parent_attributes){
	$opt_keep_attributes = 1;
}

##### MAIN ####
#### read gff file and save info in memory
######################
### Parse GFF input #
print "Reading file $opt_gfffile\n";
my ($hash_omniscient, $hash_mRNAGeneLink) = slurp_gff3_file_JD({ input => $opt_gfffile,
                                                                 config => $config
                                                              });
print "Parsing Finished\n";
### END Parse GFF input #
#########################

# extract level
my $hash_level = $hash_omniscient->{'other'}{'level'};
my $hash_l1_grouped = group_l1features_from_omniscient($hash_omniscient);

#### read fasta
my $nbFastaSeq=0;
my $db = Bio::DB::Fasta->new($opt_fastafile);
my @ids      = $db->get_all_primary_ids;
my %allIDs; # save ID in lower case to avoid cast problems
foreach my $id (@ids ){$allIDs{lc($id)}=$id;}


print ("Fasta file parsed\n");
# ----------------------------------- LEVEL 1 ----------------------------------
foreach my $seqname (sort { (($a =~ /(\d+)$/)[0] || 0) <=> (($b =~ /(\d+)$/)[0] || 0) } keys %{$hash_l1_grouped}) {

  foreach my $feature_l1 ( sort { ncmp ($a->start.$a->end.$a->_tag_value('ID'), $b->start.$b->end.$b->_tag_value('ID') ) } @{$hash_l1_grouped->{$seqname}}) {

    my $id_l1=$feature_l1->_tag_value('ID');
    my $name=undef;

    if( $opt_type eq lc($feature_l1->primary_tag()) or $opt_type eq "l1" or $opt_type eq "level1" ){

      #Handle Header
      my $id_seq = clean_string($id_l1);
      my $description.=clean_tag("seq_id=").clean_string($seqname).$OFS.clean_tag("type=").clean_string($opt_type);
			if($opt_keep_attributes){
				print "Extract attributes level1\n" if ($opt_verbose);
				my $attributes = extract_attributes($feature_l1);
				$description.=$OFS.$attributes;
			}

      my @ListSeq=($feature_l1);
      extract_sequences(\@ListSeq, $db, $id_seq, $description, $opt_full, $opt_upstreamRegion, $opt_downRegion, $opt_split, $opt_extremity_only, 'level1');
    }

# ----------------------------------- LEVEL 2 ----------------------------------
    foreach my $ptag_l2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...

      if (exists_keys ( $hash_omniscient, ('level2', $ptag_l2, lc($id_l1)) ) ){
        foreach my $feature_l2 ( @{$hash_omniscient->{'level2'}{$ptag_l2}{lc($id_l1)}}) {

          #For Header
          my $id_l2  = $feature_l2->_tag_value('ID');

          #Handle Header
          my $id_seq = clean_string($id_l2);
          my $description=clean_tag("gene=").clean_string($id_l1).$OFS.clean_tag("seq_id=").clean_string($seqname).$OFS.clean_tag("type=").clean_string($opt_type);

          if( $opt_type eq $ptag_l2 or $opt_type eq "l2" or $opt_type eq "level2" ){

						if($opt_keep_attributes ){
							print "Extract attributes level2\n" if ($opt_verbose);
							my @List_l1=($feature_l1);
							my $attributes = extract_attributes( $feature_l2, \@List_l1 );
							$description.=$OFS.$attributes;
						}

            my @ListSeq=($feature_l2);
            extract_sequences(\@ListSeq, $db, $id_seq, $description, $opt_full, $opt_upstreamRegion, $opt_downRegion, $opt_split, $opt_extremity_only, 'level2');
          }

# ----------------------------------- LEVEL 3 ----------------------------------
          # $opt_keep_attributes will extracted within the extract_sequences function
          foreach my $ptag_l3 (keys %{$hash_omniscient->{'level3'}}){
            if ( exists_keys ( $hash_omniscient, ('level3', $ptag_l3, lc($id_l2)) ) ){
              if( $opt_type eq $ptag_l3 or $opt_type eq "l3" or $opt_type eq "level3" ){
                my @list_parent_features = ($feature_l1, $feature_l2);
                extract_sequences(\@{$hash_omniscient->{'level3'}{$ptag_l3}{lc($id_l2)}}, $db, $id_seq, $description, $opt_full, $opt_upstreamRegion, $opt_downRegion, $opt_split, $opt_extremity_only, 'level3', \@list_parent_features);
              }
            }
          }
        }
      }
    }
  }
}

#END
print "usage: $0 @copyARGV\n";

if($opt_upstreamRegion and $opt_downRegion){
  print "$nbFastaSeq $opt_type converted in fasta with $opt_upstreamRegion upstream nucleotides and $opt_downRegion downstream nucleotides.\n";
}
elsif($opt_upstreamRegion){
  print "$nbFastaSeq $opt_type converted in fasta with $opt_upstreamRegion upstream nucleotides.\n";
}
elsif($opt_downRegion){
  print "$nbFastaSeq $opt_type converted in fasta with $opt_downRegion downstream nucleotides.\n";
}
else{
  print "$nbFastaSeq $opt_type converted in fasta.\n";
}

my $end_run = time();
my $run_time = $end_run - $start_run;
print "Job done in $run_time seconds\n";

#######################################################################################################################
        ####################
         #     METHODS    #
          ################
           ##############
            ############
             ##########
              ########
               ######
                ####
                 ##

# $lpa is all parent features
sub extract_attributes{
	my ($feature, $lpa) = @_;

	my $hash_attributes;
	if ($opt_keep_parent_attributes){
		my @combined;
		if(ref($feature) eq 'ARRAY'){
			@combined = (@$feature, @$lpa);
		} else {
			@combined = @$lpa;
			push @combined, $feature;
		}
		$hash_attributes = create_attribute_hash(\@combined);
	} else {
		$hash_attributes = create_attribute_hash($feature);
	}
	my $string = stringify_attribute_hash($hash_attributes);

	return $string;
}

#stringify my info HASH
sub stringify_attribute_hash{
	my ($info) = @_;

	my $string;
	foreach my $key (sort keys %{$info}){
		my $values_string=undef;
		foreach my $value (sort keys %{$info->{$key}}){
			if(! $values_string){
				$values_string=$value;
			}
			else{
				$values_string .=",".$value
			}
		}

		if(! $string){
			$string=clean_tag($key."=").clean_string($values_string);
		}
		else{
			$string .= $OFS.clean_tag($key."=").clean_string($values_string);
		}
	}
	return $string;
}

# read all attribute of a feature or a feature list a make a hash of tag => value
sub create_attribute_hash{
	my ($feature) = @_;

	# it a list of feature
	my %info;
	if(ref($feature) eq 'ARRAY'){
		foreach my $ft (@{$feature}){
			my @list_tag = $ft->get_all_tags();
			foreach my $tag (@list_tag){
				my @values = $ft->get_tag_values($tag);
				foreach my $value (@values){
					$info{$tag}{$value}++;
				}
			}
		}
	}
	# it is a single feature
	else{
		my @list_tag = $feature->get_all_tags();
		foreach my $tag (@list_tag){
			my @values = $feature->get_tag_values($tag);
			foreach my $value (@values){
				$info{$tag}{$value}++;
			}
		}
	}
	return \%info;
}

# cleanify the string to be sure to be compliant with OFS
sub clean_string{
  my ($string) = @_;

  my $replaceBy = "_";
  if($OFS eq "_"){$replaceBy = "-";}

      if($string =~ m/\Q$OFS/){
        if ($OFS eq " "){
          warn "The string <$string> contains spaces while is is used as Output Field Separator (OFS) to create fasta header, so we have quoted it (\"string\").\n".
          "If you want to keep the string/header intact, please choose another OFS using the option --ofs\n" if ! $opt_quiet;
          $string="\"".$string."\"";
        }
        else{
          warn "The fasta header has been modified !! Indeed, the string <$string> contains the Output Field Separator (OFS) <$OFS> used to build the header, so we replace it by <$replaceBy>.".
          "If you want to keep the string/header intact, please choose another OFS using the option --ofs\n" if ! $opt_quiet;
          eval "\$string =~ tr/\Q$OFS\E/\Q$replaceBy\E/";
        }
      }
  return $string
}

sub clean_tag{
  my ($string) = @_;

  my $replaceBy = "_";
  if($OFS eq "="){$replaceBy = ":";}

      if($string =~ m/\Q$OFS/){
        eval "\$string =~ tr/\Q$OFS\E/\Q$replaceBy\E/";
      }
  return $string
}

# Main function to extract sequences
# lpa = list parent attributes
sub extract_sequences{
  my($feature_list, $db, $id_seq, $description, $opt_full, $opt_upstreamRegion, $opt_downRegion, $opt_split, $opt_extremity_only, $level, $lpa )=@_;

  #sort the list
  my @sortedList = sort {$a->start <=> $b->start} @$feature_list;
  my $seq_id = $sortedList[0]->seq_id;
  #set strand, check if need to be reverse complement
  my $minus = undef;
  if($sortedList[0]->strand eq "-1" or $sortedList[0]->strand eq "-"){ $minus = 1; }
  my $phase = 0;

  # ------ Full sequence with introns ------
  if($opt_full){
    my $start = $sortedList[0]->start;
    my $end = $sortedList[$#sortedList]->end;
		# Deal with phase for CDS not starting at 0 when we want to do translation
		if( $opt_remove_orf_offset and  lc($sortedList[0]->primary_tag) eq "cds" ){
			if($minus and $sortedList[$#sortedList]->frame != 0){
				$end = $sortedList[$#sortedList]->end - $sortedList[$#sortedList]->frame;
        $phase = $sortedList[$#sortedList]->frame;
			}
			elsif (! $minus and $sortedList[0]->frame != 0){
				$start = $sortedList[0]->start + $sortedList[0]->frame;
        $phase = $sortedList[0]->frame;
			}
		}

    my $info = ""; my $right_piece = ""; my $left_piece = ""; my $sequence = "";

    # take and append the left piece if asked for
    if ( ( $opt_upstreamRegion and ! $minus ) or ( $opt_downRegion and $minus ) ){
      ($left_piece, $info) = get_left_extremity($db, $seq_id, $opt_upstreamRegion, $opt_downRegion, $minus, $start, $end, $info);
      $phase = 0; # as we add non cds sequence at the beginning we set the phase to 0
    }

    # take and append the right piece if asked for
    if ( ( $opt_downRegion and !$minus ) or ( $opt_upstreamRegion and $minus ) ){
      ($right_piece, $info) = get_right_extremity($db, $seq_id, $opt_upstreamRegion, $opt_downRegion, $minus, $start, $end, $info);
    }

    # append only extremities
    if($opt_extremity_only){
      $sequence = $left_piece.$right_piece;
      $phase = 0; # as we add non cds sequence at the beginning we set the phase to 0
    }
    else{ # append extremity to main sequence even if empty
      $sequence = get_sequence($db, $seq_id, $start, $end);
      $sequence = $left_piece.$sequence.$right_piece;
    }

		# catch attributes for Level3
		if($opt_keep_attributes and $level eq 'level3' ){ #update header's id information
			print "Extract attributes level3 full\n" if ($opt_verbose);
			my $attributes = extract_attributes(\@sortedList, $lpa);
			$description.=$OFS.$attributes;
		}

    # create object
    my $seqObj = create_seqObj($sequence, $id_seq, $description, $minus, $info);
    # print object
    print_seqObj($ostream, $seqObj, $opt_AA, $codonTable, $phase);
  }
  # --------------------------------------


   # ------ all pieces independently ------
  elsif($opt_split){

    foreach my $feature ( @sortedList ){
      my $start = $feature->start;
      my $end = $feature->end;
			# Deal with phase for CDS not starting at 0 when we want to do translation
			if( $opt_remove_orf_offset and  lc($feature->primary_tag) eq "cds" ){
				if($minus and $feature->frame != 0){
					$end = $feature->end - $feature->frame;
          $phase = $feature->frame;
				}
				elsif (! $minus and $feature->frame != 0){
					$start =$feature->start + $feature->frame;
          $phase = $feature->frame;
				}
			}

      my $info = ""; my $right_piece = ""; my $left_piece = ""; my $sequence = "";

      # take and append the left piece if asked for
      if ( ( $opt_upstreamRegion and ! $minus ) or ( $opt_downRegion and $minus ) ){
        ($left_piece, $info) = get_left_extremity($db, $seq_id, $opt_upstreamRegion, $opt_downRegion, $minus, $start, $end, $info);
        $phase = 0; # as we add non cds sequence at the beginning we set the phase to 0
      }

      # take and append the right piece if asked for
      if ( ( $opt_downRegion and !$minus ) or ( $opt_upstreamRegion and $minus ) ){
        ($right_piece, $info) = get_right_extremity($db, $seq_id, $opt_upstreamRegion, $opt_downRegion, $minus, $start, $end, $info);
      }

      # append only extremities
      if($opt_extremity_only){
        $sequence = $left_piece.$right_piece;
        $phase = 0; # as we add non cds sequence at the beginning we set the phase to 0
      }
      else{  # append extremity to main sequence even if empty
        $sequence = get_sequence($db, $seq_id, $start, $end);
        $sequence = $left_piece.$sequence.$right_piece;
      }

      my $seqObj = undef;
      if($level eq 'level3' ){ #update header's id information
        my $id_l3  = $feature->_tag_value('ID');
        my $updated_description="transcript=".$id_seq.$OFS.$description;

				# catch attributes for Level3
				if( $opt_keep_attributes ){ #update header's id information
					print "Extract attributes level3 split\n" if ($opt_verbose);
					my $attributes = extract_attributes($feature, $lpa);
					$updated_description.=$OFS.$attributes;
				}

        #create object
        $seqObj = create_seqObj($sequence, $id_l3, $updated_description, $minus, $info);
      }
      else{
        $seqObj = create_seqObj($sequence, $id_seq, $description, $minus, $info);
      }

      #print object
      print_seqObj($ostream, $seqObj, $opt_AA, $codonTable, $phase);
    }
  }
  # --------------------------------------


   # ------ feature natural ------
  else{

    my $feature_type = $sortedList[0]->primary_tag;

    # ------ SPREADED feature need to be collapsed else only if merge option activated ------
    if( exists_keys($hash_level,'spread',lc($feature_type) ) or ( $opt_merge ) ){

    	my $sequence="";my $info = "";

      # create sequence part 1
      foreach my $feature ( @sortedList ){
         $sequence .= get_sequence($db, $feature->seq_id, $feature->start, $feature->end);
      }
  		# Deal with phase for CDS not starting at 0 when we want to do translation or if asked by activation of $opt_remove_orf_offset parameter
  		if( ($opt_remove_orf_offset or $opt_AA) and  lc($sortedList[0]->primary_tag) eq "cds" ){

  			if($minus){
  				if ( $sortedList[$#sortedList]->frame eq "." ){
  					warn_no_phase();
  				}
  			 	elsif ( $sortedList[$#sortedList]->frame != 0 ){
  					$sequence = substr $sequence, 0, -$sortedList[$#sortedList]->frame; # remove offset end
            $phase = $sortedList[$#sortedList]->frame;
          }
  			}
  			else { # ! minus
  				if ( $sortedList[0]->frame  eq "." ){
  					warn_no_phase();
  				}
  				elsif( $sortedList[0]->frame != 0 ){
  					$sequence = substr $sequence, $sortedList[0]->frame; # remove offset start
            $phase = $sortedList[0]->frame;
  				}
  			}
  		}

       # update sequence with extremities if option
       if($opt_upstreamRegion or $opt_downRegion){
         my $start = $sortedList[0]->start;
         my $end = $sortedList[$#sortedList]->end;
         my $right_piece = ""; my $left_piece = "";

         # take and append the left piece if asked for
         if ( ( $opt_upstreamRegion and ! $minus ) or ( $opt_downRegion and $minus ) ){
           ($left_piece, $info) = get_left_extremity($db, $seq_id, $opt_upstreamRegion, $opt_downRegion, $minus, $start, $end, $info);
           $phase = 0; # as we add non cds sequence at the beginning we set the phase to 0
         }

         # take and append the right piece if asked for
         if ( ( $opt_downRegion and !$minus ) or ( $opt_upstreamRegion and $minus ) ){
           ($right_piece, $info) = get_right_extremity($db, $seq_id, $opt_upstreamRegion, $opt_downRegion, $minus, $start, $end, $info);
         }

         # append only extremities
         if($opt_extremity_only){
           $sequence = $left_piece.$right_piece;
           $phase = 0; # as we add non cds sequence at the beginning we set the phase to 0
         }
         else{  # append extremity to main sequence even if empty
           $sequence = $left_piece.$sequence.$right_piece;
         }
      }

			# catch attributes for Level3
			if($opt_keep_attributes and $level eq 'level3' ){ #update header's id information
				print "Extract attributes level3 natural spread merged\n" if ($opt_verbose);
				my $attributes = extract_attributes(\@sortedList, $lpa);
				$description.=$OFS.$attributes;
			}

      #create object
      my $seqObj = create_seqObj($sequence, $id_seq, $description, $minus, $info);
      #print object
      print_seqObj($ostream, $seqObj, $opt_AA, $codonTable, $phase);
    }

    # ---- Non spreaded feature extract them one by one
    else{
      foreach my $feature ( @sortedList ){
        my $start = $feature->start;
        my $end = $feature->end;
        my $info = ""; my $right_piece = ""; my $left_piece = ""; my $sequence = "";

        # take and append the left piece if asked for
        if ( ( $opt_upstreamRegion and ! $minus ) or ( $opt_downRegion and $minus ) ){
          ($left_piece, $info) = get_left_extremity($db, $seq_id, $opt_upstreamRegion, $opt_downRegion, $minus, $start, $end, $info);
        }

        # take and append the right piece if asked for
        if ( ( $opt_downRegion and !$minus ) or ( $opt_upstreamRegion and $minus ) ){
          ($right_piece, $info) = get_right_extremity($db, $seq_id, $opt_upstreamRegion, $opt_downRegion, $minus, $start, $end, $info);
        }

        # append only extremities
        if($opt_extremity_only){
          $sequence = $left_piece.$right_piece;
        }
        else{  # append extremity to main sequence even if empty
          $sequence = get_sequence($db, $seq_id, $start, $end);
          $sequence = $left_piece.$sequence.$right_piece;
        }

        my $seqObj = undef;
        if($level eq 'level3' ){ #update header's id information
          my $id_l3  = $feature->_tag_value('ID');
          my $updated_description="transcript=".$id_seq.$OFS.$description;

					# catch attributes for Level3
					if( $opt_keep_attributes ){ #update header's id information
						print "Extract attributes level3 natural not spread or spread not merged\n" if ($opt_verbose);
						my $attributes = extract_attributes($feature, $lpa);
						$updated_description.=$OFS.$attributes;
					}

					#create object
          $seqObj = create_seqObj($sequence, $id_l3, $updated_description, $minus, $info);
        }
        else{
          $seqObj = create_seqObj($sequence, $id_seq, $description, $minus, $info);
        }

        #print object
        print_seqObj($ostream, $seqObj, $opt_AA, $codonTable, $phase);
      }
    }
  }
}

# patch to avoid different warn due to different lines.
sub warn_no_phase{
	warn "No phase is specify in the CDS. We will assume it start in phase 0.";
}

# Get left extremity regardless if it is 5' or 3'
sub get_left_extremity{

  my ($db, $seq_id, $opt_upstreamRegion, $opt_downRegion, $minus, $start, $end, $info)=@_;

  if( $info ne ""){$info.=$OFS;}

  my $new_start = undef;

  if ( $minus ){
    $new_start = $start-$opt_downRegion;
    # add info left it is 3'
    if($new_start < 0){$info.=clean_tag("3'extra=").($start-1)."nt";}
    else{$info.=clean_tag("3'extra=").$opt_downRegion."nt";}
  }
  else{

    $new_start=$start-$opt_upstreamRegion;
    # add info left it is 5'
    if($new_start < 0){$info.=clean_tag("5'extra=").($start-1)."nt";}
    else{$info.=clean_tag("5'extra=").$opt_upstreamRegion."nt";}
  }

  # extract the chunck
  my $sequence = "";
  if ($new_start > $start){ # Deal with neagtive value for $opt_upstreamRegion, $opt_downRegion (e.g when trying to extract the start and stop codons from a CDS or splice sites of intron feature)
    $sequence = get_sequence($db, $seq_id, $start, $new_start-1);
  }
  else{ # Majority of cases, positive value for $opt_upstreamRegion, $opt_downRegion
    $sequence = get_sequence($db, $seq_id, $new_start, $start-1);
  }

  return $sequence, $info;
}


# Get right extremity regardless if it is 5' or 3'
sub get_right_extremity{
  my ($db, $seq_id, $opt_upstreamRegion, $opt_downRegion, $minus, $start, $end, $info)=@_;

  if( $info ne ""){$info.=$OFS;}

  my $new_end= undef;

  if ( $minus ){
    $new_end = $end+$opt_upstreamRegion;
    if($end > $db->length($seq_id) ){ $info.=clean_tag("5'extra=").($db->length($seq_id)-$end)."nt" ;}
    else{$info.=clean_tag("5'extra=").$opt_upstreamRegion."nt";}
  }
  else{
    $new_end = $end+$opt_downRegion;
    # add info right it is 3'
    if($new_end > $db->length($seq_id) ){$info.=clean_tag("3'extra=").$db->length($seq_id)-$end."nt" ;}
    else{$info.=clean_tag("3'extra=").$opt_downRegion."nt";}
  }

  # extract the chunck
  my $sequence = "";
  if ($new_end < $end){ # Deal with neagtive value for $opt_upstreamRegion, $opt_downRegion (e.g when trying to extract the start and stop codons from a CDS or splice sites of intron feature)
    $sequence = get_sequence($db, $seq_id, $new_end+1, $end);
  }
  else{ # Majority of cases, positive value for $opt_upstreamRegion, $opt_downRegion
    $sequence = get_sequence($db, $seq_id, $end+1, $new_end);
  }

  return $sequence, $info;
}


#
sub create_seqObj{
  my ($sequence, $id_seq, $description, $minus, $info)=@_;

  my $seqObj  = Bio::Seq->new( '-format' => 'fasta' , -seq => $sequence);

  #check if need to be reverse complemented (minus strand features)
  if(! $opt_plus_strand_only){ # with option plus strand you avoid to reverse complement
    $seqObj=$seqObj->revcom if $minus;
  }
  # complement asked
  if($opt_revcomp){
    $seqObj=$seqObj->revcom;
  }

  # build description
  if($info){
    $description.=$OFS.$info;
  }

  # fill object with id and description
  $seqObj->id($id_seq);
  $seqObj->description($description);

  return $seqObj;
}


# extract the sequence from the DB
sub  get_sequence{
  my  ($db, $seq_id, $start, $end) = @_;

  my $sequence="";
  my $seq_id_correct = undef;
  if ( exists_keys ( \%allIDs, (lc($seq_id)) ) ){

    $seq_id_correct = $allIDs{lc($seq_id)};

    $sequence = $db->subseq($seq_id_correct, $start, $end);

    if($sequence eq ""){
      warn "Problem ! no sequence extracted for - $seq_id !\n";  exit;
    }
    if( length($sequence) != abs($end-$start+1) ){
      my $wholeSeq = $db->subseq($seq_id_correct);
      $wholeSeq = length($wholeSeq);
      warn "Problem ! The size of the extracted sequence ".length($sequence)." is different than the specified span: ".abs($end-$start+1).
      ".\nThat often occurs when the fasta file does not correspond to the annotation file. Or the index file comes from another fasta file which had the same name and haven't been removed.\n".
      "As last possibility your gff contains location errors (Already encountered for a Maker annotation)\n",
      "Supplement information: seq_id=$seq_id ; seq_id_correct=$seq_id_correct ; start=$start ; end=$end ; $seq_id sequence length: $wholeSeq )\n";
    }
  }
  else{
    warn "Problem ! ID $seq_id not found !\n";
  }

  return $sequence;
}

# Print the sequence object
sub print_seqObj{
  my($ostream, $seqObj, $opt_AA, $codonTable, $phase) = @_;


  if($opt_AA){ #translate if asked

			if ( length($seqObj->seq()) < 3 ){warn "Sequence to translate for ".$seqObj->id()." < 3 nucleotides! Skipped...\n"; return; }

      my $transObj = $seqObj->translate(-CODONTABLE_ID => $codonTable);

      # translate alternative start codon by a M
      my $start_codon = substr($seqObj->seq(),0,3); # get start codon
      my $myCodonTable  = Bio::Tools::CodonTable->new( -id => $codonTable );
      my $first_AA = substr($transObj->seq(),0,1);
      if ($phase == 0 and $opt_alternative_start_codon and $myCodonTable->is_start_codon($start_codon)){
        if($first_AA ne "M"){ # if the start codon was not a M while it is a valid start codon we have to replace it by a methionine
          my $translated_seq = substr($transObj->seq(),1); # removing first AA
          $transObj->seq("M".$translated_seq);  # adding M as first AA
          print "Replacing valid alternative start codon (AA=$first_AA) by a methionine (AA=M) for ".$seqObj->id().".\n";
        }
      }

      if($opt_cleanFinalStop and $opt_cleanInternalStop){ #this case is needed to be able to remove two final stop codon in a raw when the both options are activated.
        my $lastChar = substr $transObj->seq(),-1,1;
        my $cleanedSeq=$transObj->seq();
        if ($lastChar eq "*"){ # if last char is a stop we remove it
          chop $cleanedSeq;
        }
        $cleanedSeq =~ tr/*/X/; #X = Any / unknown Amino Acid
        $transObj->seq($cleanedSeq);
      }

      elsif($opt_cleanFinalStop){
		    my $lastChar = substr $transObj->seq(),-1,1;

        if ($lastChar eq "*"){ # if last char is a stop we remove it
		      my $cleanedSeq=$transObj->seq();
		      chop $cleanedSeq;
		      $transObj->seq($cleanedSeq);
		    }
      }

      elsif($opt_cleanInternalStop){
        my $lastChar = substr $transObj->seq(),-1,1;

        my $seqMinus1=$transObj->seq();
        chop $seqMinus1;
        $seqMinus1 =~ tr/*/X/; #X = Any / unknown Amino Acid
        my $cleanedSeq=$seqMinus1.$lastChar;
        $transObj->seq($cleanedSeq);
      }

			if ( length($transObj->seq()) == 0 ){warn "Translated sequence empty for ".$transObj->id().". The nucleotide sequence was only coding for stop codon(s). The translated stop codon(s) '*' has/have been removed due to stop codon(s)'s cleaning  parameter(s) activated. Skipped...\n"; return; }

      $ostream->write_seq($transObj);
   }
  else{
    $ostream->write_seq($seqObj);
  }

	$nbFastaSeq++;
}


__END__

=head1 NAME

agat_sp_extract_sequences.pl

=head1 DESCRIPTION

This script extracts sequences in fasta format according to features described
in a gff file. You can extract the fasta of any type of feature. The feature
type is defined within the 3rd column in the gff file.
The result is written to the specified output file, or to STDOUT.

Features spanning several locations (e.g. UTR, CDS), are extracted chunk by chunk
and merged to create the biological feature. If you wish to extract each chunck independently,
please refer to the --split parameter. To see the list of features that may span over several locations
within AGAT run: agat levels --expose
and then look at section called spread in the feature_levels.yaml file.

The headers are formated like that:
>ID gene=gene_ID seq_id=Chromosome_ID type=cds 5'extra=VALUE tag=value tag=value

The ID is the identifier of the feature (ID attribute in the 9th column.
If missing it is created by AGAT)
The gene value will be the ID of the level1 feature (the top feature of the record)
The seq_id value is the value from 1st column within the gff.
The type value holds the information of the feature type extracted (3rd column in the gff).
5'extra or 3'extra is optional, it holds the information of extra nucleotides
removed or added when using the downstream and/or upstream parameter.
The tag=value are optional and represent attributes from the gff (See --keep_attributes parameter).

The OFS of all values can be modified excepted for the ID (see --ofs parameter).
In such case the tool gives a warning.

Some examples:
To extract the coding regions:
agat_sp_extract_sequences.pl -g infile.gff -f infile.fasta
or
agat_sp_extract_sequences.pl -g infile.gff -f infile.fasta -t cds
To extract and translate the coding regions:
agat_sp_extract_sequences.pl -g infile.gff -f infile.fasta -p
or
agat_sp_extract_sequences.pl -g infile.gff -f infile.fasta -t cds -p
To extract the mRNA (biological definition UTR+CDS):
agat_sp_extract_sequences.pl -g infile.gff -f infile.fasta -t exon --merge
To extract each cds chunk independently:
agat_sp_extract_sequences.pl -g infile.gff -f infile.fasta -t cds --split
To extract 5'UTR with introns:
agat_sp_extract_sequences.pl -g infile.gff -f infile.fasta -t "5'UTR" --full
To extract 100nt upstream region of a gene:
agat_sp_extract_sequences.pl -g infile.gff -f infile.fasta -t gene --upstream 100


    agat_sp_extract_sequences.pl -g infile.gff -f infile.fasta  [ -o outfile ]
    agat_sp_extract_sequences.pl --help

=head1 OPTIONS

=over 8

=item B<--alternative_start_codon> or B<--asc>

Bolean - When activated it can affect the translation of the start codon.
Indeed alternative start codons exist, and are translated by the cells'machinery
by a Methionine (M). By default AGAT translates the first codon as other codons by the
corresponding AA. If you wish to translate the first codon by a M when it is a valid
alternative start codon, activate this parameter.
If the sequence you try to translate is a CDS (or starting by a CDS), the phase
is checked and the alternative start codon is accepted only if the phase is 0.

=item B<--cdna>

Boolean - This extract the cdna sequence (i.e reverse complement of the mRNA:
transcribed sequence (devoid of introns, but containing untranslated exons,
then reverse complemented). It corresponds to extract the exons sequences,
merge them, and reverse complement the sequence (--type exon --merge --revcomp).

=item B<--clean_final_stop> or B<--cfs>

Boolean - The Clean Final Stop option allows removing the translation of the
final stop codons that is represented by the <*> character.
This character can be disturbing for many programs (e.g interproscan)

=item B<--clean_internal_stop> or B<--cis>

Boolean - The Clean Internal Stop option allows replacing the translation of the
stop codons present among the sequence that is represented by the <*> character
by <X>. Indeed the <*> character can be disturbing for many programs
(e.g interproscan)

=item B<--codon>, B<--table> or B<--ct>

Integer - Allow to choose the codon table for the translation. [default 1]

=item B<--do>, B<-3>, B<--three>, B<-down> or B<-downstream>

Integer - It will take that number of nucleotide in more at the 3' extremity.
/!\ You must activate the option "--full" if you with to extract only the most downstream part of certain feature (exon,cds,utr)
otherwise you will extract each downstream parts of the subfeatures (e.g many cds parts may be needed to shape a cds in its whole).

=item B<--eo>

Boolean - Called 'extremity only', this option will extract only the adjacent parts of a feature.
This option has to be activated with -u and/or -p option.
/!\ using -u and -p together builds a chimeric sequence which will be the concatenation of the left and right extremities of a feature.

=item B<-f> or B<--fasta>

String - Input fasta file.

=item B<--full>

Boolean - This option allows dealing with feature that may span over several locations
like cds or exon, in order to extract the full sequence from the start extremity
of the first chunck to the end extremity of the last chunk.
The use of that option with '--type exon' will extract the pre-mRNA sequence (i.e with introns).
Use of that option on cds will give the pre-mRNA without the untraslated regions (UTRs).
(To extract an mRNA as it is defined biologicaly you need to use the
`-t exon` option with the --merge option)

=item B<-g>, B<--gff> or B<-ref>

String - Input GTF/GFF file.

=item B<-h> or B<--help>

Boolean - Display this helpful text.

=item B<--keep_attributes>

Boolean - The value of the attribute tags will be extracted from the feature type
specified by the option --type and stored in the fasta header.

=item B<--keep_parent_attributes>

Boolean - The value of the attribute tags will be extracted from the feature type
specified by the option --type along with those from parent features and stored in the fasta header.

=item B<--merge>

Boolean - By default, only features that span several locations (e.g. cds and utr can
span over several exons) are merged together. In order to merge other
type of features (e.g. exon) you must activate this parameter.

=item B<--mrna>

Boolean - This extract the mrna sequence (i.e transcribed sequence (devoid of
introns, but containing untranslated exons)). It corresponds to extract the exon
sequences and merge them (--type exon --merge).

=item B<--ofs>

String - Output Fields Separator for the description field of the header of the
created fasta sequences. By default it's a space < > but can be modified by any String or
character using this option.

=item B<-o> or B<--output>

String - Output fasta file.  If no output file is specified, the output will be
written to STDOUT.

=item B<--plus_strand_only>

Boolean - By default the extrated feature sequences from a minus strand is
reverse complemented. Activating this option you will always get sequence from plus
strand ( not reverse complemented).
You can get the opposite (minus strand only) by using --plus_strand_only --revcomp

=item B<-p>, B<--protein> or B<--aa>

Boolean - Will translate the extracted sequence in Amino acid.
By default the codon table used is the 1 (Standard).
See --table parameter for more options.

=item B<--remove_orf_offset> or B<--roo>

Boolean - CDS can start with a phase different from 0 when a gene model is fragmented.
When asking for protein translation this (start) offset is trimmed out automatically.
But when you extract CDS dna sequences, this  (start) offset is not removed by default.
To remove it activate this option. If --up or --do option are used too, the (start) offset
is trimmed first, then is added the piece of sequence asked for.

=item B<--revcomp>

Boolean - To reverse complement the extracted sequence [default - False].
By default the extrated feature sequences from a minus strand is
reverse complemented. Consequently, for minus strand features that option will
extract the sequences from plus strand from left to right.

=item B<--split>

Boolean -  By default, all features that span several locations
(e.g. cds and utr can span over several exons) are merge together to shape
the biological feature (e.g. several cds chuncks are merged to create the CDS
 in its whole).
If you wish to extract all the chuncks independently activate this option.

=item B<-t> or B<--type>

String - Define the feature you want to extract the sequence from.
Default 'cds'.
Most common choice are: gene,mrna,exon,cds,trna,three_prime_utr,five_prime_utr.
When you choose exon (or cds,utr,etc.), all the exons of a same parent feature
are attached together before to extract the sequence. If you wish to extract each
exon of an mRNA independently, see option --split.
/!\ `-t mRNA` will extract the features labeled as "mRNA" and corresponds to the cdna*
because it contains the introns if any. It does not actually extract the mRNAs as
it is defined biologicaly. To extract the mRNA as defined biologicaly you must use `-t exon`.
*Not a real cdna because it is not reversed

=item B<--up>, B<-5>, B<--five> or B<-upstream>

Integer - It will take that number of nucleotide in more at the 5' extremity.
/!\ You must activate the option "--full" if you wish to extract only the most
upstream part of certain features (exon,cds,utr)
otherwise you will extract each upstream parts of the subfeatures
(e.g many cds parts may be needed to shape a cds in its whole).

=item B<-c> or B<--config>

String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any, 
otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
The --config option gives yo the possibility to use your own AGAT config file (located elsewhere or named differently).

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
