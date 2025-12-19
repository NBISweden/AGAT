#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use Carp;
use Getopt::Long;
use IO::File;
use Pod::Usage;
use AGAT::AGAT;

start_script();
my $header = get_agat_header();
# ---------------------------- OPTIONS ----------------------------
my $opt_file;
my $opt_output=undef;
my $Xsize=10;
my $opt_help = 0;

# ---------------------------- OPTIONS ----------------------------
# OPTION MANAGEMENT: partition @ARGV into shared vs script options via library
my ($shared_argv, $script_argv) = split_argv_shared_vs_script(\@ARGV);

# Parse script-specific options from its own list
my $script_parser = Getopt::Long::Parser->new;
$script_parser->configure('bundling','no_auto_abbrev');
if ( ! $script_parser->getoptionsfromarray(
  $script_argv,
  'f|gff|ref|reffile=s' => \$opt_file,
  'o|out|output=s'      => \$opt_output,
  'i|intron_size=i'     => \$Xsize,
  'h|help!'             => \$opt_help ) )
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

if ( ! defined($opt_file) ) {
    pod2usage( {
           -message => "$header\nMust specify at least 1 parameters:\nReference data gff3 file (--gff)\n",
           -verbose => 0,
           -exitval => 1 } );
}

# --- Manage config ---
my $shared_opts = parse_shared_options($shared_argv);
initialize_agat({ config_file_in => ( $shared_opts->{config} ), input => $opt_file, shared_opts => $shared_opts });

# ----------------------------------------------------------------------------

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    PARAMS    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

my $ostreamReport_file;
if (defined($opt_output) ) {
  my ($filename,$path,$ext) = fileparse($opt_output,qr/\.[^.]*/);
  $ostreamReport_file = $path.$filename."_report.txt";
}

my $gffout = prepare_gffout( $opt_output );
my $ostreamReport = prepare_fileout($ostreamReport_file);

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>     MAIN     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

######################
### Parse GFF input #
my ($hash_omniscient) = slurp_gff3_file_JD({ input => $opt_file });

### END Parse GFF input #
#########################

my $nb_cases=0;
my $tag = "pseudo";
######################
### Parse GFF input #
foreach my $tag_l1 (keys %{$hash_omniscient->{'level1'}}){
  foreach my $id_l1 (keys %{$hash_omniscient->{'level1'}{$tag_l1}}){
    my $shortest_intron=10000000000;
    foreach my $tag_l2 (keys %{$hash_omniscient->{'level2'}}){
      if (exists_keys($hash_omniscient,('level2',$tag_l2,$id_l1) ) ){
        foreach my $feature_l2 (@{$hash_omniscient->{'level2'}{$tag_l2}{$id_l1}}){
          my $level2_ID = lc($feature_l2->_tag_value('ID'));

          if ( exists_keys($hash_omniscient,('level3','cds',$level2_ID)) ){
            my $counterL3=-1;
            my $indexLast = $#{$hash_omniscient->{'level3'}{'cds'}{$level2_ID}};
            my @sortedList = sort {$a->start <=> $b->start} @{$hash_omniscient->{'level3'}{'cds'}{$level2_ID}};
            foreach my $feature_l3 ( @sortedList ){
              #count number feature of tag_l3 type
              $counterL3++;
              #Manage Introns## from the second intron to the last (from index 1 to last index of the table sortedList) ## We go inside this loop only if we have more than 1 feature.
              if($counterL3 > 0 and $counterL3 <= $indexLast){
                my $intronSize = $sortedList[$counterL3]->start - $sortedList[$counterL3-1]->end;
                $shortest_intron = $intronSize if($intronSize < $shortest_intron)
              }
            }
          }
        }
      }
    }
    dual_print2 "Shortest intron for $id_l1:".$shortest_intron."\n" if ($shortest_intron != 10000000000);
    if ($shortest_intron < $Xsize){
      dual_print1 "flag the gene $id_l1\n";
      $nb_cases++;

      my $feature_l1 = $hash_omniscient->{'level1'}{$tag_l1}{$id_l1};
      $feature_l1->add_tag_value($tag, $shortest_intron);
      if($feature_l1->has_tag('product') ){
        $feature_l1->add_tag_value('note', $feature_l1->get_tag_values('product'));
        $feature_l1->remove_tag('product');
      }
      foreach my $tag_l2 (keys %{$hash_omniscient->{'level2'}}){
        if (exists_keys ($hash_omniscient, ('level2', $tag_l2, $id_l1) ) ) {
          foreach my $feature_l2 (@{$hash_omniscient->{'level2'}{$tag_l2}{$id_l1}}){
            my $level2_ID = lc($feature_l2->_tag_value('ID'));
            $feature_l2->add_tag_value($tag, $shortest_intron);
            if($feature_l2->has_tag('product') ){
              $feature_l2->add_tag_value('note', $feature_l2->get_tag_values('product'));
              $feature_l2->remove_tag('product');
            }

            foreach my $tag_l3 (keys %{$hash_omniscient->{'level3'}}){
              if ( exists_keys($hash_omniscient, ('level3', $tag_l3, $level2_ID) ) ){
                foreach my $feature_l3 (@{$hash_omniscient->{'level3'}{$tag_l3}{$level2_ID}}){
                  $feature_l3->add_tag_value($tag, $shortest_intron);
                  if($feature_l3->has_tag('product') ){
                    $feature_l3->add_tag_value('note', $feature_l3->get_tag_values('product'));
                    $feature_l3->remove_tag('product');
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

my $toprint = "We found $nb_cases cases where introns were < $Xsize, we flagged them with the attribute $tag. The value of this tag is size of the shortest intron found in this gene.\n";
print $ostreamReport $toprint if($opt_output);
dual_print1 "$toprint";

print_omniscient( {omniscient => $hash_omniscient, output => $gffout} );

# --- final messages ---
end_script();

# ----------------------------------------------------------------------------


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



__END__


=head1 NAME

agat_sp_flag_short_introns_ebi.pl

=head1 DESCRIPTION

The script flags records that contain short introns (default 10bp) within coding sequences (CDS) with the <pseudo> attribute and changes the <product> attribute into a <note> attribute.
This is useful for avoiding ERROR messages when submitting data to the EBI.
(Typical EBI error message: ERROR: Intron usually expected to be at least 10 nt long. Please check the accuracy.)

=head1 SYNOPSIS

    agat_sp_flag_short_introns_ebi.pl --gff infile --out outfile
    agat_sp_flag_short_introns_ebi.pl --help

=head1 OPTIONS

=over 8

=item B<--gff>, B<-f>, B<--ref> or B<-reffile> <file>

Input GTF/GFF file.

=item  B<--intron_size> or B<-i> <int>

Minimum intron size, default 10. All genes with an intron < of this size will be
flagged with the pseudo attribute (the value will be the size of the smallest
intron found within the incriminated gene)

=item  B<-o>, B<--out> or B<--output> <file>

Output file to create (default GFF3 - see config to modify output format).
If no output file is specified, the output will be written to STDOUT.

=item B<--help> or B<-h>

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
