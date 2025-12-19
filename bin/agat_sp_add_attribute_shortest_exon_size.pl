#!/usr/bin/env perl

use strict;
use warnings;
use POSIX qw(strftime);
use File::Basename;
use Carp;
use Getopt::Long;
use IO::File;
use Pod::Usage;
use AGAT::AGAT;

start_script();
my $header = get_agat_header();
# -------------------------------- LOAD OPTIONS --------------------------------
my $opt_file;
my $opt_output = undef;
my $opt_help = 0;

# OPTION MANAGEMENT: partition @ARGV into shared vs script options via library
my ($shared_argv, $script_argv) = split_argv_shared_vs_script(\@ARGV);

# Parse script-specific options from its own list
my $script_parser = Getopt::Long::Parser->new;
$script_parser->configure('bundling', 'no_auto_abbrev');
if ( ! $script_parser->getoptionsfromarray(
    $script_argv,
    'f|gff|ref=s'     => \$opt_file,
    'o|out|output=s'  => \$opt_output,
    'h|help!'         => \$opt_help,
  ) ) {
  pod2usage({
    -message => 'Failed to parse command line',
    -verbose => 1,
    -exitval => 1
  });
}

# Print Help and exit
if ($opt_help) {
  pod2usage({ -verbose => 99,
        -exitval => 0,
        -message => "$header\n" });
}

if ( ! defined($opt_file) ) {
  pod2usage({
    -message => "$header\nMust specify at least 1 parameters:\nReference data gff3 file (--gff)\n",
    -verbose => 0,
    -exitval => 1
  });
}

# Parse shared options
my ($shared_opts) = parse_shared_options($shared_argv);

# --- Manage config ---
initialize_agat({ config_file_in => ($shared_opts->{config}), input => $opt_file, shared_opts => $shared_opts });
# -----------------------------------------------------------------------------------------------

# #######################
# # START Manage Option #
# #######################
my $ostreamReport_filename;
if (defined($opt_output) ) {
  my ($filename,$path,$ext) = fileparse($opt_output,qr/\.[^.]*/);
  $ostreamReport_filename=$path.$filename."_report.txt";
}
my $gffout = prepare_gffout( $opt_output);
my $ostreamReport = prepare_fileout($ostreamReport_filename);

                                                      #######################
                                                      #        MAIN         #
#                     >>>>>>>>>>>>>>>>>>>>>>>>>       #######################       <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

######################
### Parse GFF input #
my ($hash_omniscient) = slurp_gff3_file_JD({ input => $opt_file });
### END Parse GFF input #
#########################

my $nb_cases_l1=0;
my $nb_cases_l2=0;
my $tag = "shortest_exon";
######################
### Parse GFF input #
foreach my $tag_l1 (keys %{$hash_omniscient->{'level1'}}){
  foreach my $id_l1 (keys %{$hash_omniscient->{'level1'}{$tag_l1}}){
    my $feature_l1 = $hash_omniscient->{'level1'}{$tag_l1}{$id_l1};
    my $shortest_exon=undef;
    foreach my $tag_l2 (keys %{$hash_omniscient->{'level2'}}){
      if (exists_keys($hash_omniscient,('level2',$tag_l2,$id_l1) ) ){
        my $shortest_exon_l2=undef;
        foreach my $feature_l2 (@{$hash_omniscient->{'level2'}{$tag_l2}{$id_l1}}){
          my $level2_ID = lc($feature_l2->_tag_value('ID'));

          if ( exists_keys($hash_omniscient,('level3','exon',$level2_ID)) ){

            foreach my $feature_l3 (  @{$hash_omniscient->{'level3'}{'exon'}{$level2_ID}} ){
              my $exonSize = $feature_l3->end - $feature_l3->start+1;
              if(! $shortest_exon_l2){
                $shortest_exon_l2 = $exonSize;
              }
              elsif( $exonSize < $shortest_exon_l2){
                $shortest_exon_l2 = $exonSize;
              }
            }
            $feature_l2->add_tag_value($tag, $shortest_exon_l2);
            $nb_cases_l2++;
            if(! $shortest_exon){
              $shortest_exon = $shortest_exon_l2;
            }
            elsif( $shortest_exon_l2 < $shortest_exon){
              $shortest_exon = $shortest_exon_l2;
            }
          }
        }
      }
    }
    if($shortest_exon){
      $feature_l1->add_tag_value($tag, $shortest_exon);
      $nb_cases_l1++;
    }
  }
}

my $toprint = "$nb_cases_l1 $tag flags/attributes added to level1 features and $nb_cases_l2 $tag flags/attributes added to level2 features. The value of the attribute is size of the shortest exon found.\n";
print $ostreamReport $toprint;
dual_print1 $toprint;
print_omniscient( {omniscient => $hash_omniscient, output => $gffout} );
      #########################
      ######### END ###########
      #########################

# --- final messages ---
end_script($ostreamReport);


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

agat_sp_add_attribute_shortest_exon_size.pl

=head1 DESCRIPTION

The script add the attribute <shortest_exon> to each gene and rna, which will hold the size of the shortest exon in bp.

=head1 SYNOPSIS

    agat_sp_add_attribute_shortest_exon_size.pl --gff infile --out outfile
    agat_sp_add_attribute_shortest_exon_size.pl --help

=head1 OPTIONS

=over 8

=item B<--gff>, B<-f> or B<--ref> <file>

Input GTF/GFF file.

=item B<-o>, B<--out> or B<--output> <file>

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
