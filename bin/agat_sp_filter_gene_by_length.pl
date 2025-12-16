#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Pod::Usage;
use IO::File;
use AGAT::AGAT;

start_script();
my $header = get_agat_header();
# -----------------------------------------------------------------------------------------------
my $opt_test=">";
my $opt_output= undef;
my $opt_size = 100;
my $opt_gff = undef;
my $opt_help;

# ---------------------------- OPTIONS ----------------------------
# Partition @ARGV into shared vs script options
my ($shared_argv, $script_argv) = split_argv_shared_vs_script(\@ARGV);

# Parse script-specific options
my $script_parser = Getopt::Long::Parser->new;
$script_parser->configure('bundling','no_auto_abbrev');
if ( ! $script_parser->getoptionsfromarray(
  $script_argv,
                  'f|ref|reffile|gff=s' => \$opt_gff,
                  't|test=s'            => \$opt_test,
                  "s|size=i"            => \$opt_size,
                  'o|output=s'          => \$opt_output,
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

if ( ! $opt_gff ){
    pod2usage( {
           -message => "$header\nAt least 1 parameter is mandatory:\n1) Input reference gff file: --gff\n\n",
           -verbose => 0,
           -exitval => 2 } );
}

# --- Manage config ---
# Parse shared options and initialize AGAT
my ($shared_opts) = parse_shared_options($shared_argv);
initialize_agat({ config_file_in => ( $shared_opts->{config} ), input => $opt_gff, shared_opts => $shared_opts });

# -----------------------------------------------------------------------------------------------

###############
# Manage Output

## FOR GFF FILE
my $gffout_ok_file ;
my $gffout_notok_file ;
my $ostreamReport_file ;

if ($opt_output) {
  my ($outfile,$path,$ext) = fileparse($opt_output,qr/\.[^.]*/);

  # set file names
  $gffout_ok_file = $path.$outfile.$ext;
  $gffout_notok_file = $path.$outfile."_remaining".$ext;
  $ostreamReport_file = $path.$outfile."_report.txt";
}

my $gffout_ok = prepare_gffout( $gffout_ok_file );
my $gffout_notok = prepare_gffout( $gffout_notok_file );
my $ostreamReport = prepare_fileout( $ostreamReport_file );

#Manage test option
if($opt_test ne "<" and $opt_test ne ">" and $opt_test ne "<=" and $opt_test ne ">=" and $opt_test ne "="){
  die "The test to apply is Wrong: $opt_test.\nWe want something among this list: <,>,<=,>= or =.";
}

# start with some interesting information
my $stringPrint = "We will select l1 feature (e.g. gene) that have length $opt_test $opt_size bp.\n";

print $ostreamReport $stringPrint if ($opt_output);
dual_print1 $stringPrint;

                          #######################
# >>>>>>>>>>>>>>>>>>>>>>>>#        MAIN         #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                          #######################

######################
### Parse GFF input #
my ($hash_omniscient) =  slurp_gff3_file_JD({ input => $opt_gff });
### END Parse GFF input #
#########################
# sort by seq id
my $hash_sortBySeq = gather_and_sort_l1_by_seq_id($hash_omniscient);

my @listok;
my @listNotOk;
#################
# == LEVEL 1 == #
#################
foreach my $seqid (sort { (($a =~ /(\d+)$/)[0] || 0) <=> (($b =~ /(\d+)$/)[0] || 0) } keys %{$hash_sortBySeq}){ # loop over all the feature level1

	foreach my $tag_l1 (sort {$a cmp $b} keys %{$hash_omniscient->{'level1'}}){
		foreach my $feature_l1 ( @{$hash_sortBySeq->{$seqid}{$tag_l1}} ){
			my $id_l1 = lc($feature_l1->_tag_value('ID'));
      my $gene_length=$feature_l1->end()-$feature_l1->start()+1;
      my $successl1 = test_size( $gene_length, $opt_test, $opt_size );
      my $longer_concat_exon=undef;
	    #################
	    # == LEVEL 2 == #
	    #################
	    foreach my $tag_l2 (sort keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
	      if ( exists_keys( $hash_omniscient, ('level2', $tag_l2, $id_l1) ) ){
	        foreach my $feature_l2 ( @{$hash_omniscient->{'level2'}{$tag_l2}{$id_l1}} ) {
            my $id_l2 = lc($feature_l2->_tag_value('ID'));
	          #################
	          # == LEVEL 3 == #
	          #################
            foreach my $tag_l3 (sort keys %{$hash_omniscient->{'level3'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
  	          if ( exists_keys( $hash_omniscient, ('level3', 'exon', $id_l2) ) ){
                my $local_size=0;
                foreach my $feature_l3 ( @{$hash_omniscient->{'level3'}{'exon'}{$id_l2}} ) {
                  $local_size += $feature_l3->end()-$feature_l3->start()+1;
                }
                if($longer_concat_exon and $longer_concat_exon<$local_size){
                  $longer_concat_exon = $local_size;
                }
                elsif(! $longer_concat_exon) {
                  $longer_concat_exon = $local_size;
                }
              }
	          }
	        }
	      }
	    }
      # case we had exon (we look at the longest mRNA)
      if($longer_concat_exon){
        dual_print2("$id_l1 does have exon(s). Longest concatenated exons: $longer_concat_exon\n");
        if( test_size( $longer_concat_exon, $opt_test, $opt_size ) ){
          dual_print2("$id_l1 pass the test\n");
          push @listok, $id_l1;
        }
        else{
          dual_print2("$id_l1 do not pass the test\n");
          push @listNotOk, $id_l1;
        }
      }
      else{
        dual_print2("$id_l1 does not have any exon. $tag_l1 size: $gene_length\n");
        # No exon, L1 pass test
        if($successl1){
          dual_print2("$id_l1 pass the test\n");
          push @listok, $id_l1;
        }
        # No exon, L1 do not pass test
        else{
          dual_print2("$id_l1 do not pass the test\n");
          push @listNotOk, $id_l1;
        }
      }
    }
  }
}

# print ok
my $hash_ok = subsample_omniscient_from_level1_id_list_delete($hash_omniscient, \@listok);
print_omniscient( {omniscient => $hash_ok, output => $gffout_ok} );
%{$hash_ok} = ();
# print remaining if an output is provided
if($opt_output){
  my $hash_remaining = subsample_omniscient_from_level1_id_list_delete($hash_omniscient, \@listNotOk);
  print_omniscient( {omniscient => $hash_remaining, output => $gffout_notok} );
  %{$hash_remaining} = ();
}

my $test_success = scalar @listok;
my $test_fail = scalar @listNotOk;

my $stringPrint2 = "$test_success l1 feature (e.g. gene) selected with a length $opt_test $opt_size bp.\n";
$stringPrint2 .= "$test_fail remaining l1 feature (e.g. gene) do not pass the test.\n";
print $ostreamReport $stringPrint2 if ($opt_output);
dual_print1 $stringPrint2;

# --- final messages ---
end_script();

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

sub test_size{
  my ($size, $operator, $nb_ref) = @_;

  if ($operator eq ">"){
    if ($size > $nb_ref){
      return "true";
    }
  }
  if ($operator eq "<"){
    if ($size < $nb_ref){
      return "true";
    }
  }
  if ($operator eq "=" or $operator eq "=="){
    if ($size == $nb_ref){
      return "true";
    }
  }
  if ($operator eq "<="){
    if ($size <= $nb_ref){
      return "true";
    }
  }
  if ($operator eq ">="){
    if ($size >= $nb_ref){
      return "true";
    }
  }
  return undef;
}

__END__

=head1 NAME

agat_sp_filter_gene_by_length.pl

=head1 DESCRIPTION

The script aims to filter level1 feature (e.g. gene, match, etc) by length.
It will create two files. one with the feature passing the length filter,
the other one with the remaining features.
If the level1 feature has exon features, the size is computed by concatenating
the exon together. If the level1 feature has several level2 features (e.g. mRNA)
we apply the test on the longest one (the longest concatenated exon set).

Some examples:
Select L1 feature shorter than 1000bp:
agat_sp_filter_gene_by_length.pl --gff infile.gff  --size 1000 --test "<" -o result.gff
Select genes longer than 200bp:
agat_sp_filter_gene_by_length.pl --gff infile.gff --size 200 --test ">" -o result.gff

=head1 SYNOPSIS

    agat_sp_filter_gene_by_length.pl --gff infile.gff --test ">=" --nb 10 [ --output outfile ]
    agat_sp_filter_gene_by_length.pl --help

=head1 OPTIONS

=over 8

=item B<-f>, B<--reffile>, B<--gff>  or B<-ref>

Input GFF3 file that will be read

=item B<-s> or B<--size>

Integer - Gene size in pb [Default 100]

=item B<-t> or B<--test>

Test to apply (>, <, =, >= or <=). If you use one of these two characters >, <,
please do not forget to quote your parameter like that "<=". Else your terminal will complain.
[Default "="]

=item B<-o> or B<--output>

Output GFF file.  If no output file is specified, the output will be
written to STDOUT.

=item B<-h> or B<--help>

Display this helpful text.

=back

=head1 SHARED OPTIONS

Shared options are defined in the AGAT configuration file and can be overridden via the command line for this script only.
Common shared options are listed below; for the full list, please refer to the AGAT agat_config.yaml.

=over 8

=item B<--config>

String - Path to a custom AGAT configuration file.  
By default, AGAT uses `agat_config.yaml` from the working directory if present, otherwise the default file shipped with AGAT
(available locally via `agat config --expose`).

=item B<--cpu>, B<--core>, B<--job> or B<--thread>

Integer - Number of parallel processes to use for file input parsing (via forking).

=item B<-v> or B<--verbose>

Integer - Verbosity, choice are 0,1,2,3,4. 0 is quiet, 1 is normal, 2,3,4 is more verbose. Default 1.

=back

=head1 FEEDBACK

For questions, suggestions, or general discussions about AGAT, please use the AGAT community forum:
https://github.com/NBISweden/AGAT/discussions

=head1 BUG REPORTING

Bug reports should be submitted through the AGAT GitHub issue tracker:
https://github.com/NBISweden/AGAT/issues

=cut

AUTHOR - Jacques Dainat
