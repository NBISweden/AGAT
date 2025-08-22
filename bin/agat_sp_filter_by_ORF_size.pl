#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use POSIX qw(strftime);
use Getopt::Long::Descriptive;
use Clone 'clone';
use Pod::Usage;
use AGAT::AGAT;

my $start_run = time();
my $header    = get_agat_header();
my @copyARGV  = @ARGV;
my ( $opt, $usage, $config ) = AGAT::AGAT::describe_script_options(
    $header,
    [ 'gff|g|f|ref|reffile=s', 'Input GTF/GFF file', { required => 1 } ],
    [ 'test|t=s',              'Comparison test (<,>,<=,>=,==,=)',
      { default => '>',
        callbacks => {
            valid => sub {
                $_[0] =~ /^(?:<|>|<=|>=|==|=)$/
                  or die 'Test to apply must be one of <, >, <=, >=, == or =';
                return 1;
            }
        }
      }
    ],
    [ 'size|s=i', 'Protein length threshold', { default => 100 } ],
);

my $gff         = $opt->gff;
my $opt_test    = $opt->test;
my $PROT_LENGTH = $opt->size;
my $outfile     = $config->{output};

my $log;
if ( my $log_name = $config->{log_path} ) {
    open( $log, '>', $log_name )
      or die "Can not open $log_name for printing: $!";
    dual_print( $log, $header, 0 );
}

######################
# Option check

if($opt_test){
  if($opt_test ne "<" and $opt_test ne ">" and $opt_test ne "<=" and $opt_test ne ">=" and $opt_test ne "=" and $opt_test ne "=="){ 
    dual_print( $log, "The test to apply is Wrong: $opt_test.\nWe want something among this list: <,>,<=,>=,== or =.\n", 1 );
    exit;
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

my $gffout_pass = prepare_gffout($config, $gffout_pass_file);
my $gffout_notpass = prepare_gffout($config, $gffout_notpass_file);

# print usage performed
my $stringPrint = strftime "%m/%d/%Y at %Hh%Mm%Ss", localtime;
$stringPrint = "Launched the ".$stringPrint."\nusage: $0 @copyARGV\n";
$stringPrint .= "We are filtering the gene with protein size $opt_test $PROT_LENGTH\n";
dual_print( $log, $stringPrint, $config->{verbose} );

                #####################
                #     MAIN          #
                #####################

######################
### Parse GFF input #
  my ($hash_omniscient, $hash_mRNAGeneLink) = slurp_gff3_file_JD({ input => $gff,
                                                                   config => $config
                                                                 });
  dual_print( $log, "GFF3 file parsed\n", $config->{verbose} );

# Create an empty omniscient hash to store the discarded features and copy the config in
my %hash_omniscient_discarded;
    $hash_omniscient_discarded{'config'} = clone($hash_omniscient->{'config'});

my $number_mRNA_discarded=0;
my $number_gene_discarded=0;
my $number_gene_affected=0;
foreach my $primary_tag_l1 (keys %{$hash_omniscient->{'level1'}}){ # primary_tag_l1 = gene or repeat etc...
  foreach my $gene_id_l1 (keys %{$hash_omniscient->{'level1'}{$primary_tag_l1}}){
    
    my $gene_feature=$hash_omniscient->{'level1'}{$primary_tag_l1}{$gene_id_l1};
      dual_print( $log, "Study gene $gene_id_l1\n", $config->{verbose} );
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
            dual_print( $log, "Case all L2 discarded \n", $config->{verbose} );
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
            dual_print( $log, "Case some L2 discarded \n", $config->{verbose} );
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
          dual_print( $log, "No cds for $gene_id_l1\n", $config->{verbose} );
        }
      }
      # ---------- CASE NO L2 -----------
      if($no_l2){ # case of l1 feature without child
        dual_print( $log, "No child for $gene_id_l1\n", $config->{verbose} );
      }
    }
  }
}

#resume
print_omniscient( {omniscient => $hash_omniscient, output => $gffout_pass} );
print_omniscient( {omniscient => \%hash_omniscient_discarded, output => $gffout_notpass} );

dual_print( $log, "\n$number_gene_affected genes have at least one transcript removed.\n", $config->{verbose} );
dual_print( $log, "$number_gene_discarded genes discarded\n", $config->{verbose} );
dual_print( $log, "$number_mRNA_discarded transcripts discarded.\n", $config->{verbose} );

# END
my $end_run = time();
my $run_time = $end_run - $start_run;
dual_print( $log, "Job done in $run_time seconds\n", $config->{verbose} );
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

=item B<-g> or B<--gff>

Input GTF/GFF file.

=item B<-s> or B<--size>

ORF size to apply the test. Default 100.

=item B<-t> or B<--test>
Test to apply (> < = >= <=). If you us one of these two character >, <, please don't forget to quote you parameter like that "<=" otherwise your terminal will complain.
By default it will be ">"

=item B<-v>

Verbose. Useful for debugging purpose. Bolean

=item B<-o> or B<--out> or B<--output> or B<--outfile>

Output GFF file.  If no output file is specified, the output will be
written to STDOUT.

=item B<-c> or B<--config>

String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any, 
otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

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
