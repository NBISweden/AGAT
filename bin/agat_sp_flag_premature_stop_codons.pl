#!/usr/bin/env perl

use strict;
use warnings;
use POSIX qw(strftime);
use File::Basename;
use Sort::Naturally;
use Carp;
use Getopt::Long;
use IO::File;
use Pod::Usage;
use Bio::DB::Fasta;
use AGAT::AGAT;


my $header = get_agat_header();
my $config;
my $opt_file;
my $opt_output;
my $file_fasta;
my $codonTable = 1;
my $opt_help = 0;

my @copyARGV=@ARGV;
if ( !GetOptions(	'gff|ref|reffile=s' => \$opt_file,
                	'o|out|output=s'    => \$opt_output,
					"fasta|fa|f=s"      => \$file_fasta,
					"table|codon|ct=i"  => \$codonTable,
                 	'c|config=s'        => \$config,
                 	'h|help!'           => \$opt_help ) )
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

if ( !$opt_file or !$file_fasta) {
    pod2usage( {
           -message => "$header\nMust specify at least 2 parameters:\n".
					 "Reference data gff3 file (--gff)\n".
					 "Reference data fasta file (--gff)\n",
           -verbose => 0,
           -exitval => 1 } );
}

# --- Manage config ---
initialize_agat({ config_file_in => $config, input => $opt_file });

# --- Check codon table
$codonTable = get_proper_codon_table($codonTable);

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    PARAMS    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

my $ostreamReport_file;
if (defined($opt_output) ) {
  my ($filename,$path,$ext) = fileparse($opt_output,qr/\.[^.]*/);

  $ostreamReport_file = $path.$filename."_report.txt";
}

my $gffout = prepare_gffout( $opt_output );
my $ostreamReport = prepare_fileout($ostreamReport_file);

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    EXTRA     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# Print info
my $string1 = strftime "%m/%d/%Y at %Hh%Mm%Ss", localtime;
$string1 .= "\n\nusage: $0 @copyARGV\n\n";

print $ostreamReport $string1;
if($opt_output){print $string1;}

# activate warnings limit
my %warnings;
activate_warning_limit(\%warnings, 10);

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>     MAIN     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

######################
### Parse GFF input #
print "Reading ".$opt_file,"\n";
my ($hash_omniscient, $hash_mRNAGeneLink) = slurp_gff3_file_JD({ input => $opt_file });
### END Parse GFF input #
#########################

####################
# index the genome #
my $db_fasta = Bio::DB::Fasta->new($file_fasta);
# save ID in lower case to avoid cast problems
my %allIDs;
my @ids_db_fasta     = $db_fasta->get_all_primary_ids;
foreach my $id (@ids_db_fasta ){$allIDs{lc($id)}=$id;}
print ("Fasta file parsed\n");


my $nb_cases=0;
my $nb_cases_l1=0;
# sort by seq id
my ( $hash_sortBySeq, $hash_sortBySeq_std, $hash_sortBySeq_topf ) = collect_l1_info_sorted_by_seqid_and_location($hash_omniscient);

# Read by seqId to sort properly the output by seq ID
foreach my $seqid (sort { (($a =~ /(\d+)$/)[0] || 0) <=> (($b =~ /(\d+)$/)[0] || 0) } keys %{$hash_sortBySeq}){ # loop over all the feature level1

	foreach my $locationid ( sort { ncmp ($a, $b) } keys %{$hash_sortBySeq->{$seqid} } ){

		my $primary_tag_l1 = $hash_sortBySeq->{$seqid}{$locationid}{'tag'};
		my $id_l1 = $hash_sortBySeq->{$seqid}{$locationid}{'id'};
		my $feature_l1 = $hash_omniscient->{'level1'}{$primary_tag_l1}{$id_l1};

		my $nb_this_l2_pseudo=0;
		foreach my $tag_l2 (keys %{$hash_omniscient->{'level2'}}){
			if (exists_keys ( $hash_omniscient,('level2',$tag_l2,$id_l1) ) ){
				foreach my $feature_l2 ( @{$hash_omniscient->{'level2'}{$tag_l2}{lc($id_l1)}}) {
					my $level2_ID = lc($feature_l2->_tag_value('ID'));

		      if ( exists_keys($hash_omniscient,('level3','cds', $level2_ID)) ){
		        my @sortedList = sort {$a->start <=> $b->start} @{$hash_omniscient->{'level3'}{'cds'}{lc($level2_ID)}};
						my $sequence = "";

						# check strand
						my $minus = undef;
						if($sortedList[0]->strand eq "-1" or $sortedList[0]->strand eq "-"){ $minus = 1; }

						# create sequence
						foreach my $feature ( @sortedList ){
							 $sequence .= get_sequence($db_fasta, $feature->seq_id, $feature->start, $feature->end);
						}

						# Deal with offset
						my $start_position=0;
						if($minus){
							$start_position = $sortedList[$#sortedList]->end();

							if ( $sortedList[$#sortedList]->frame eq "." ){
								warn_no_phase();
							}
							elsif ( $sortedList[$#sortedList]->frame != 0 ){
								$sequence = substr $sequence, 0, -$sortedList[$#sortedList]->frame; # remove offset end
								$start_position -= $sortedList[$#sortedList]->frame;
							}
						}
						else { # ! minus
							$start_position = $sortedList[0]->start();

							if ( $sortedList[0]->frame  eq "." ){
								warn_no_phase();
							}
							elsif( $sortedList[0]->frame != 0 ){
								$sequence = substr $sequence, $sortedList[0]->frame; # remove offset start
								$start_position += $sortedList[0]->frame;
							}
						}

						#create sequence object
						 my $seqObj  = Bio::Seq->new( '-format' => 'fasta' , -seq => $sequence);

						 #check if need to be reverse complement
						 $seqObj=$seqObj->revcom if $minus;
						if ( length($seqObj->seq()) < 3 ){warn "Sequence to translate for ".$seqObj->id()." < 3 nucleotides! Skipped...\n"; return; }

						# translate
						my $transObj = $seqObj->translate(-CODONTABLE_ID => $codonTable);

						# get protein sequence
						my $seqMinus1=$transObj->seq();
						# remove last character
						chop $seqMinus1;
						# count if there is any stop codon
						my $count = $seqMinus1 =~ tr/*//;

						if ($count){

							my @positions;
							my $position = $start_position;
							foreach my $char (split //, $seqMinus1) {
								if($char eq "*"){
									push @positions, $position;
								}
								$position = $minus ? $position-3 : $position+3; # follow nt position codon per codon
							}
							$feature_l2->add_tag_value('pseudo', @positions);
							my $arrSize = @positions;
							my $toprint = "We flag the $tag_l2 $level2_ID that contained $arrSize premature stop codons\n";
							print $ostreamReport $toprint; print $toprint;
						  $nb_cases++;
							$nb_this_l2_pseudo++;
						}
		      }
		    }
				if($nb_this_l2_pseudo eq $#{$hash_omniscient->{'level2'}{$tag_l2}{lc($id_l1)}} + 1){
					$feature_l1->add_tag_value('pseudo', "yes");
					$nb_cases_l1++;
				}
				elsif($nb_this_l2_pseudo){
					print "Not all isoforms of $id_l1 are pseudogenes, so we do not flag the gene.\n";
				}
			}
	  }
	}
}

my $toprint = "We found $nb_cases cases where mRNAs contain premature stop codons. They have been flagged as pseudogene.\n".
							"$nb_cases_l1 genes have been flagged as pseudogene.\n";
print $ostreamReport $toprint;
if($opt_output){print $toprint;}

print_omniscient( {omniscient => $hash_omniscient, output => $gffout} );

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
     warn "Problem ! The size of the sequence extracted ".length($sequence)." is different than the specified span: ".abs($end-$start+1).
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

# warning
sub warn_no_phase{
	warn "No phase is specify in the CDS. We will assume it start in phase 0.";
}
__END__


=head1 NAME

agat_sp_flag_premature_stop_codons.pl

=head1 DESCRIPTION

The script flags the mRNAs containing premature stop codons.
It will add the attribute "pseudo" and the value will be the positions of all premature stop codons.
Gene are flagged as pseudogene only if all the isoforms are pseudogenes. The attribute
will also be "pseudo" but will not contains any location.

=head1 SYNOPSIS

    agat_sp_flag_premature_stop_codons.pl --gff infile.gff --fasta infile.fa --out outfile
    agat_sp_flag_premature_stop_codons.pl --help

=head1 OPTIONS

=over 8

=item B<--gff>, B<--ref> or B<-reffile>

Input GTF/GFF file.

=item B<-f>, B<--fa> or B<--fasta>

Imput fasta file.

=item B<--ct>, B<--codon> or B<--table>

Codon table to use. [default 1]

=item  B<--out>, B<--output> or B<-o>

Output gff3 file where the result will be printed.

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
