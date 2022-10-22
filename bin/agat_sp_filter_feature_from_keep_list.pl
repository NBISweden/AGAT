#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use POSIX qw(strftime);
use Pod::Usage;
use IO::File;
use AGAT::AGAT;

my $header = get_agat_header();
my $config = get_agat_config();
my $primaryTag=undef;
my $opt_output= undef;
my $opt_keep_list = undef;
my $opt_attribute = 'ID';
my $opt_gff = undef;
my $opt_verbose = undef;
my $opt_help;

# OPTION MANAGMENT
my @copyARGV=@ARGV;
if ( !GetOptions( 'f|ref|reffile|gff=s' => \$opt_gff,
                  'kl|keep_list=s'      => \$opt_keep_list,
                  "p|type|l=s"          => \$primaryTag,
                  'o|output=s'          => \$opt_output,
                  'a|attribute=s'       => \$opt_attribute,
                  'v|verbose!'          => \$opt_verbose,
                  'h|help!'             => \$opt_help ) )
{
    pod2usage( { -message => 'Failed to parse command line',
                 -verbose => 1,
                 -exitval => 1 } );
}

if ($opt_help) {
    pod2usage( { -verbose => 99,
                 -exitval => 0,
                 -message => "$header\n" } );
}

if ( ! $opt_gff or ! $opt_keep_list ){
    pod2usage( {
           -message => "$header\nAt least 2 parameters are mandatory:\n1) Input reference gff file: --gff\n".
           "2) A keep list (one value per line): --keep_list\n\n",
           -verbose => 0,
           -exitval => 2 } );
}

###############
# Manage Output

## FOR GFF FILE
my $gffout_ok_file ;
my $ostreamReport_file;

if ($opt_output) {
  my ($outfile,$path,$ext) = fileparse($opt_output,qr/\.[^.]*/);

  # set file names
  $gffout_ok_file = $path.$outfile.$ext;
  $ostreamReport_file = $path.$outfile."_report.txt";
}
my $gffout_ok = prepare_gffout($config, $gffout_ok_file);
my $ostreamReport = prepare_fileout($ostreamReport_file);

# Manage $primaryTag
my @ptagList;
my $print_feature_string;
if(! $primaryTag or $primaryTag eq "all"){
  $print_feature_string = "all features";
  push(@ptagList, "all");
}
elsif($primaryTag =~/^level[123]$/){
  $print_feature_string .= "$primaryTag features ";
  push(@ptagList, $primaryTag);
}
else{
   @ptagList= split(/,/, $primaryTag);
   foreach my $tag (@ptagList){
      if($tag =~/^level[123]$/){
        $print_feature_string .= "$primaryTag features ";
      }
      else{
        $print_feature_string .= "$tag feature ";
      }
   }
}

# Manage keep List
my %keep_hash;
open my $in_keep_list, "<:encoding(utf8)", $opt_keep_list or die "$opt_keep_list: $!";
while (my $line = <$in_keep_list>) {
    chomp $line;
    $line =~ s/^\s+|\s+$//g; #removing leading and trailing white spaces
    $keep_hash{lc($line)}++;
}
my $nb_to_keep = keys %keep_hash;

# start with some interesting information
my $stringPrint = strftime "%m/%d/%Y at %Hh%Mm%Ss", localtime;
$stringPrint .= "\nusage: $0 @copyARGV\n";
$stringPrint .= "We will keep the records that have $print_feature_string sharing the value of the $opt_attribute attribute with the keep list.\n";
$stringPrint .= "The keep list contains $nb_to_keep uniq IDs\n";

if ($opt_output){
  print $ostreamReport $stringPrint;
  print $stringPrint;
}
else{ print $stringPrint; }
                          #######################
# >>>>>>>>>>>>>>>>>>>>>>>>#        MAIN         #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                          #######################
my %all_cases = ('l1' => 0, 'l2' => 0, 'l3' => 0, 'all' => 0);
######################
### Parse GFF input #
my ($hash_omniscient, $hash_mRNAGeneLink) =  slurp_gff3_file_JD({ input => $opt_gff,
                                                                  config => $config
                                                                });
print("Parsing Finished\n");
### END Parse GFF input #
#########################
# sort by seq id
my $hash_sortBySeq = gather_and_sort_l1_by_seq_id($hash_omniscient);

my $keepit=undef;
my @keeplist=();
#################
# == LEVEL 1 == #
#################
foreach my $seqid (sort { (($a =~ /(\d+)$/)[0] || 0) <=> (($b =~ /(\d+)$/)[0] || 0) } keys %{$hash_sortBySeq}){ # loop over all the feature level1

	foreach my $tag_l1 (sort {$a cmp $b} keys %{$hash_omniscient->{'level1'}}){
		foreach my $feature_l1 ( @{$hash_sortBySeq->{$seqid}{$tag_l1}} ){
			my $id_l1 = lc($feature_l1->_tag_value('ID'));

		  $keepit = check_feature($feature_l1, 'level1', \@ptagList);
			# we can remove feature L1 now because we are looping over $hash_sortBySeq not $hash_omniscient
		  if ($keepit){
		  	push @keeplist, $id_l1;
		   	next;
		  }

		  #################
		  # == LEVEL 2 == #
		  #################
			my @list_l2_to_remove;
			my @list_l3_to_remove;
		  foreach my $tag_l2 (sort keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...

		    if ( exists_keys( $hash_omniscient, ('level2', $tag_l2, $id_l1) ) ){
			    my @list_fl2 = @{$hash_omniscient->{'level2'}{$tag_l2}{$id_l1}};
	        foreach my $feature_l2 ( @list_fl2 ) {

            $keepit = check_feature($feature_l2,'level2', \@ptagList);
				    if ($keepit){
				      push @keeplist, $id_l1;
				      last;
				    }
				    #################
				    # == LEVEL 3 == #
				    #################
			      my $id_l2 = lc($feature_l2->_tag_value('ID'));

				    foreach my $tag_l3 (sort keys %{$hash_omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...
				      if ( exists_keys( $hash_omniscient, ('level3', $tag_l3, $id_l2) ) ){
				       	my @list_fl3 = @{$hash_omniscient->{'level3'}{$tag_l3}{$id_l2}};
				       	foreach my $feature_l3 ( @list_fl3 ) {

					        $keepit = check_feature($feature_l3, 'level3', \@ptagList);
					        if ($keepit){
					          push @keeplist, $id_l1;
					          last;
					        }
					      }
					      last if $keepit; #We found one l3 we can skip the rest of l3
					    }
					  }
					  last if $keepit; #We found one l3 we can skip the rest of l2
				  }
				  last if $keepit; #We found one l2 we can skip the rest of l2
		    }
			}
		}
	}
}

# create omniscient with only selected recoreds
my $hash_kept = subsample_omniscient_from_level1_id_list_delete($hash_omniscient, \@keeplist);
print_omniscient( {omniscient => $hash_kept, output => $gffout_ok} );#print gene modified in file
$stringPrint = ($#keeplist+1)." records kept!\n";
if ($opt_output){
  print $ostreamReport $stringPrint;
  print $stringPrint;
} else{ print $stringPrint; }

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

sub check_feature{
  my  ($feature, $level, $ptagList)=@_;

  my $keepit=undef;
  my $primary_tag=$feature->primary_tag;
	if($feature->has_tag($opt_attribute)){

	  # check primary tag (feature type) to handle
	  foreach my $ptag (@$ptagList){

	    if($ptag eq "all"){
	      $keepit = 1 if( exists_keys(\%keep_hash, (lc($feature->_tag_value($opt_attribute)))));
	    }
	    elsif(lc($ptag) eq $level){
	      $keepit = 1 if( exists_keys(\%keep_hash, (lc($feature->_tag_value($opt_attribute)))));
	    }
	    elsif(lc($ptag) eq lc($primary_tag) ){
	      $keepit = 1 if( exists_keys(\%keep_hash, (lc($feature->_tag_value($opt_attribute)))));
	    }
	  }
	}
	else{
		warn "No attribute $opt_attribute found for the following feature:\n".$feature->gff_string."\n";
	}
  return $keepit;
}

__END__

=head1 NAME

agat_sp_filter_feature_from_keep_list.pl

=head1 DESCRIPTION

The script aims to keep records based on a keeplist.
The default behaviour is to look at the features's ID. If the feature has an ID
(case insensitive) listed among the keeplist it will be kept along with all
related features (the whole record is kept. A record repsent all features linked
 by relationship e.g. gene+transcript+exon+cds of a same locus).

=head1 SYNOPSIS

    agat_sp_filter_feature_from_keep_list.pl --gff infile.gff --keep_list file.txt  [ --output outfile ]
    agat_sp_filter_feature_from_keep_list.pl --help

=head1 OPTIONS

=over 8

=item B<-f>, B<--reffile>, B<--gff>  or B<-ref>

Input GFF3 file that will be read

=item B<-p>,  B<--type> or  B<-l>

primary tag option, case insensitive, list. Allow to specied the feature types that will be handled.
You can specified a specific feature by given its primary tag name (column 3) as: cds, Gene, MrNa
You can specify directly all the feature of a particular level:
      level2=mRNA,ncRNA,tRNA,etc
      level3=CDS,exon,UTR,etc
By default all feature are taking into account. fill the option by the value "all" will have the same behaviour.


=item B<--kl> or B<--keep_list>

Keep list. One value per line.

=item  B<-a> or B<--attribute>

Attribute tag to specify the attribute to analyse. Case sensitive. Default: ID


=item B<-o> or B<--output>

Output GFF file.  If no output file is specified, the output will be
written to STDOUT.

=item B<-v>

Verbose option for debugging purpose.

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
