#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long;
use Clone 'clone';
use Pod::Usage;
use AGAT::AGAT;

start_script();
my $header = get_agat_header();
# -----------------------------------------------------------------------------------------------
my $PROT_LENGTH = 100;
my $file_fasta=undef;
my $outfile = undef;
my $opt_test = undef;
my $gff = undef;
my $opt_help= 0;

# Partition @ARGV into shared vs script options
my ($shared_argv, $script_argv) = split_argv_shared_vs_script(\@ARGV);

# Parse script-specific options
my $script_parser = Getopt::Long::Parser->new;
$script_parser->configure('bundling','no_auto_abbrev');
if ( ! $script_parser->getoptionsfromarray(
    $script_argv,
    "h|help"                 => \$opt_help,
    "g|gff=s"                => \$gff,
    't|test=s'               => \$opt_test,
    "size|s=i"               => \$PROT_LENGTH,
    "output|out|o=s" => \$outfile))

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
           -message => "$header\nAt least 1 parameter is mandatory:\n Input reference gff file (--gff)\n\n",
           -verbose => 0,
           -exitval => 1 } );
}

my ($shared_opts) = parse_shared_options($shared_argv);
initialize_agat({ config_file_in => $shared_opts->{config}, input => $gff, shared_opts => $shared_opts });

# -----------------------------------------------------------------------------------------------

######################
# Option check

if($opt_test){
  if($opt_test ne "<" and $opt_test ne ">" and $opt_test ne "<=" and $opt_test ne ">=" and $opt_test ne "=" and $opt_test ne "=="){
    die "The test to apply is Wrong: $opt_test.\nWe want something among this list: <,>,<=,>=,== or =.";
  }
}
else{
  $opt_test = ">";
}

# To avoid > < character in output files
my $opt_test_to_print = $opt_test;
$opt_test_to_print =~ s/>/sup/ig;
$opt_test_to_print =~ s/</inf/ig;

######################
# Manage output file #
my $gffout_pass_file;
my $gffout_notpass_file;
if ($outfile) {
  $outfile=~ s/.gff//g;
  $gffout_pass_file = $outfile."_".$opt_test_to_print.$PROT_LENGTH.".gff";
  $gffout_notpass_file = $outfile."_NOT_".$opt_test_to_print.$PROT_LENGTH.".gff";
}

my $gffout_pass = prepare_gffout( $gffout_pass_file);
my $gffout_notpass = prepare_gffout( $gffout_notpass_file);

# print 
dual_print1 "We are filtering the gene with protein size $opt_test $PROT_LENGTH";

                #####################
                #     MAIN          #
                #####################

######################
### Parse GFF input #
my ($hash_omniscient) = slurp_gff3_file_JD({ input => $gff });

# Create an empty omniscient hash to store the discarded features and copy the config in
my %hash_omniscient_discarded;
    $hash_omniscient_discarded{'config'} = clone($hash_omniscient->{'config'});

my $number_mRNA_discarded=0;
my $number_gene_discarded=0;
my $number_gene_affected=0;
foreach my $primary_tag_l1 (keys %{$hash_omniscient->{'level1'}}){ # primary_tag_l1 = gene or repeat etc...
  foreach my $gene_id_l1 (keys %{$hash_omniscient->{'level1'}{$primary_tag_l1}}){
    
    my $gene_feature=$hash_omniscient->{'level1'}{$primary_tag_l1}{$gene_id_l1};
    dual_print2 "Study gene $gene_id_l1";
		my $no_l2=1;# see if standalone or topfeature

    foreach my $primary_tag_l2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
       
      if ( exists_keys( $hash_omniscient, ('level2', $primary_tag_l2, $gene_id_l1) ) ){

        $no_l2 = undef;  
        my $there_is_cds=undef;
        my @l2_to_discard; 
        my @l2_to_keep; 

        foreach my $level2_feature ( @{$hash_omniscient->{'level2'}{$primary_tag_l2}{$gene_id_l1}}) {

          # get level2 id
          my $id_level2 = lc($level2_feature->_tag_value('ID'));

          ##############################
          #If it's a mRNA with a CDS. #
          if( exists_keys( $hash_omniscient, ('level3', 'cds', $id_level2 ) ) ) {
            $there_is_cds="true";

            # Manage the CDS 
            my $cds_size=0;
            foreach my $cds (@{$hash_omniscient->{'level3'}{'cds'}{$id_level2}}){
              $cds_size+= ($cds->end() - $cds->start() + 1);
            }
            $cds_size = ($cds_size - 3) / 3; # Remove the stop codon and divide by 3 to get Amnino acid
            # test the CDS
            if( test_size($cds_size, $PROT_LENGTH, $opt_test) ){
              push @l2_to_keep, $id_level2;
            } else {
              $number_mRNA_discarded++;
              push @l2_to_discard, $id_level2;
            }
          } 
          else {
            $number_mRNA_discarded++;
            push @l2_to_discard, $id_level2;
          }
        }

        # ---------- CASE there is at least one CDS -----------
        if( $there_is_cds ){
          # All transcript discarded
          if( @l2_to_keep == 0){
            dual_print2 "Case all L2 discarded \n";
            $number_gene_discarded++;
            $number_gene_affected++;
            # move L3
            foreach my $level2_ID (@l2_to_discard){ 
              foreach my $primary_tag_l3 (keys %{$hash_omniscient->{'level3'}}){ # 
                if ( exists ($hash_omniscient->{'level3'}{$primary_tag_l3}{$level2_ID} ) ){
                  $hash_omniscient_discarded{'level3'}{$primary_tag_l3}{$level2_ID} = delete $hash_omniscient->{'level3'}{$primary_tag_l3}{$level2_ID} 
                }
              }
						}
            # move L2
            $hash_omniscient_discarded{'level2'}{$primary_tag_l2}{$gene_id_l1} = delete $hash_omniscient->{'level2'}{$primary_tag_l2}{$gene_id_l1};
            # move L1
            $hash_omniscient_discarded{'level1'}{$primary_tag_l1}{$gene_id_l1} = delete $hash_omniscient->{'level1'}{$primary_tag_l1}{$gene_id_l1};
          }
          # Only part of the isoforms have been discarded
          elsif ( @l2_to_discard > 0){
            dual_print2 "Case some L2 discarded \n";
            $number_gene_affected++;
            # handle L3
            
            foreach my $level2_ID (@l2_to_discard){
              foreach my $primary_tag_l3 (keys %{$hash_omniscient->{'level3'}}){ # 
                if ( exists ($hash_omniscient->{'level3'}{$primary_tag_l3}{$level2_ID} ) ){
                  $hash_omniscient_discarded{'level3'}{$primary_tag_l3}{$level2_ID} = delete $hash_omniscient->{'level3'}{$primary_tag_l3}{$level2_ID} 
                }
              }
            }
            # handle L2
            my @new_l2_flist_keep;
            my @new_l2_flist_discard;
            if (exists_keys ($hash_omniscient, ('level2', $primary_tag_l2, $gene_id_l1) ) ){ 
              foreach my $feature_l2 ( @{$hash_omniscient->{'level2'}{$primary_tag_l2}{$gene_id_l1}}) {
                my $id_level2 = lc($feature_l2->_tag_value('ID'));
                if ( grep {$_ eq $id_level2} @l2_to_discard ){
                  push @new_l2_flist_discard, $feature_l2 ;
                } else {
                  push @new_l2_flist_keep, $feature_l2 ;
                }
              }
            }
            $hash_omniscient_discarded{'level2'}{$primary_tag_l2}{$gene_id_l1} = \@new_l2_flist_discard;
            $hash_omniscient->{'level2'}{$primary_tag_l2}{$gene_id_l1} = \@new_l2_flist_keep;
            # handle L1
            $hash_omniscient_discarded{'level1'}{$primary_tag_l1}{$gene_id_l1} = clone( $hash_omniscient->{'level1'}{$primary_tag_l1}{$gene_id_l1} );
          }
        }
        # ---------- CASE there is no CDS -----------
        else{
          dual_print2 "No cds for $gene_id_l1\n";
        }
      }
      # ---------- CASE NO L2 -----------
      if($no_l2){ # case of l1 feature without child
        dual_print2 "No child for $gene_id_l1\n";
      }
    }
  }
}

#resume
print_omniscient( {omniscient => $hash_omniscient, output => $gffout_pass} );
print_omniscient( {omniscient => \%hash_omniscient_discarded, output => $gffout_notpass} );

dual_print1 "\n$number_gene_affected genes have at least one transcript removed.\n";
dual_print1 "$number_gene_discarded genes discarded\n";
dual_print1 "$number_mRNA_discarded transcripts discarded.\n";

# END
end_script();

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

sub test_size{
  my ($size, $PROT_LENGTH, $operator) = @_;

  if ($operator eq ">"){
    if ($size > $PROT_LENGTH){
      return "true";
    }
  }
  if ($operator eq "<"){
    if ($size < $PROT_LENGTH){
      return "true";
    }
  }
  if ($operator eq "=" or $operator eq "=="){
    if ($size == $PROT_LENGTH){
      return "true";
    }
  }
  if ($operator eq "<="){
    if ($size <= $PROT_LENGTH){
       return "true";
    }
  }
  if ($operator eq ">="){
    if ($size >= $PROT_LENGTH){
      return "true";
    }
  }
return undef;
}

__END__

=head1 NAME

agat_sp_filter_by_ORF_size.pl

=head1 DESCRIPTION

The script reads a gff annotation file, and create two output files,
one contains the gene models with ORF passing the test, the other contains the rest.
By default the test is "> 100" that means all gene models that have ORF longer
than 100 Amino acids, will pass the test.
In the case of isoforms, the isoforms that do not pass the test are removed
(If all isoforms are removed, the gene is removed).
A gene with with any transcript having any CDS will be considered as non
coding gene and will not be removed.

=head1 SYNOPSIS

    agat_sp_filter_by_ORF_size.pl --gff infile.gff [ -o outfile ]
    agat_sp_filter_by_ORF_size.pl -h

=head1 OPTIONS

=over 8

=item B<-g> or B<--gff> <file>

Input GTF/GFF file.

=item B<-s> or B<--size> <int>

ORF size to apply the test. Default 100.

=item B<-t> or B<--test> <operator>

Test to apply (> < = >= <=). If you us one of these two character >, <, please don't forget to quote you parameter like that "<=" otherwise your terminal will complain.
By default it will be ">"

=item B<-o>, B<--out> or B<--output> <file>

Output file to create (default GFF3 - see config to modify output format).
If no output file is specified, the output will be written to STDOUT.

=item B<-h> or B<--help>

Display this helpful text.

=back

=head1 SHARED OPTIONS

Shared options are defined in the AGAT configuration file and can be overridden via the command line for this script only.
Common shared options are listed below; for the full list, please refer to the AGAT agat_config.yaml.

=over 8

=item B<--config> <file>

Path to a custom AGAT configuration file.  
By default, AGAT uses `agat_config.yaml` from the working directory if present, otherwise the default file shipped with AGAT
(available locally via `agat config --expose`).

=item B<--cpu>, B<--core>, B<--job> or B<--thread> <int>

Number of parallel processes to use for file input parsing (via forking).

=item B<-v> or B<--verbose> <int>

Verbosity, choice are 0,1,2,3,4. 0 is quiet, 1 is normal, 2,3,4 is more verbose. Default 1.

=back

=head1 FEEDBACK

For questions, suggestions, or general discussions about AGAT, please use the AGAT community forum:
https://github.com/NBISweden/AGAT/discussions

=head1 BUG REPORTING

Bug reports should be submitted through the AGAT GitHub issue tracker:
https://github.com/NBISweden/AGAT/issues

=cut

AUTHOR - Jacques Dainat
