#!/usr/bin/env perl

use strict;
use warnings;
use List::MoreUtils  qw(natatime);;
use Carp;
use Getopt::Long;
use Pod::Usage;
use Clone 'clone';
use AGAT::AGAT;

start_script();
my $header = get_agat_header();
# -------------------------------- LOAD OPTIONS --------------------------------
my $intergenicID = 1;
my $opt_file;
my $opt_output = undef;
my $opt_help = 0;
my @copyARGV = @ARGV;

# OPTION MANAGEMENT: partition @ARGV into shared vs script options via library
my ($shared_argv, $script_argv) = split_argv_shared_vs_script(\@ARGV);

# Parse script-specific options from its own list
my $script_parser = Getopt::Long::Parser->new;
$script_parser->configure('bundling', 'no_auto_abbrev');
if ( ! $script_parser->getoptionsfromarray(
    $script_argv,
    'f|gff|ref=s'        => \$opt_file,
    'o|out|output=s'     => \$opt_output,
    'h|help!'            => \$opt_help,
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
    -message => "$header\nMust specify at least 1 parameters:\nReference data GFF/GTF file (--gff)\n",
    -verbose => 0,
    -exitval => 1
  });
}

# Parse shared options
my ($shared_opts) = parse_shared_options($shared_argv);

# --- Manage config ---
initialize_agat({ config_file_in => ($shared_opts->{config}), input => $opt_file, shared_opts => $shared_opts });
# -----------------------------------------------------------------------------------------------

# # START Manage Option #
# #######################
my $gffout = prepare_gffout( $opt_output );

#                         #######################
# >>>>>>>>>>>>>>>>>>>>>>>>#        MAIN         #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#                         #######################

  ######################
  ### Parse GFF input #
  my ($omniscient) = slurp_gff3_file_JD({ input => $opt_file });
  ### END Parse GFF input #
  #########################

if(! exists_keys($omniscient,('level1', "gene") ) ){
  dual_print1 "No gene feature found in $opt_file, intergenic regions cannot be determinded!";
  exit 0;
}

# Gather all Level1 features
my $sortBySeq = gather_and_sort_l1_by_seq_id_for_l1type($omniscient, 'gene');

# --------------------------- COLLECT GENE LOCATIONS -----------------------
dual_print2 "Now colleting the gene locations\n";
my $flattened_locations = {};
foreach my $locusID ( sort keys %{$sortBySeq}){ # tag_l1 = gene or repeat etc...
  # check if gene  exits for this sequence
  if( exists_keys($sortBySeq,($locusID, "gene") ) ){

    # Go through location from left to right ### !! if not empty
    while ( my $gene_feature = shift @{$sortBySeq->{$locusID}{"gene"}} ){

      # Define location l1
      my $current_location_l1 = [$gene_feature->start(), $gene_feature->end()];

      # save appending a list
      push @{$flattened_locations->{$locusID}}, $current_location_l1 ;
    }
  }
}

# --------------------------- FIX OVERLAPPING LOCATIONS -----------------------
# Will merge locations that overlap

dual_print2 "Now flattening the locations\n";
foreach my $locusID (  keys %{$flattened_locations} ){

  my @newlocations;
  my $previous_location = undef;
  foreach my $location ( sort {$a->[0] <=> $b->[0]} @{$flattened_locations->{$locusID}} ){

    # first round
    if (! $previous_location){
        push @newlocations, $location;
        $previous_location = $location;
    }
    # Not first round
    else{
      #  location A  -------------------------
      #  location B           -------------------------
      if ( ($previous_location->[0] <= $location->[1]) and ($previous_location->[1] >= $location->[0])){

        #  location A  -------------------------
        #  location B           -------------------------
        if($previous_location->[1] <= $location->[1]){
          # take back last location pushed in the array
          $previous_location = pop @newlocations;
          $previous_location  = [$previous_location->[0], $location->[1]];
          # push back into the array the previous location that has been modified
          push @newlocations, $previous_location ;
        }
      }
      else{
        push @newlocations, $location ;
        $previous_location = $location;
      }
    }
  }
  @{$flattened_locations->{$locusID}} = @newlocations ;
}

# --------------------------- NOW creating intergenic location -----------------------
dual_print2 "Now creating intergenic regions\n";
my $intergenic_added=0;
# Go through location from left to right ### !! if not empty
foreach my $locusID ( sort keys %{$flattened_locations}){ # tag_l1 = gene or repeat etc...
  my $round = 0;
  my $previous_location = undef;
  my $last_location = undef;

  while ( my $location = shift @{$flattened_locations->{$locusID}} ){
    $round++;
     
    $last_location = 1 if (scalar @{$flattened_locations->{$locusID}} == 0);
    if($round == 1 or $last_location){
      $previous_location=$location;
      next;
    }

    my $uniq_ID = "intergenic_region_".$intergenicID;
    my $intergenic_feature = Bio::SeqFeature::Generic->new( -seq_id => $locusID, 
                                                            -source_tag => "AGAT",
                                                            -primary_tag => 'intergenic_region' , 
                                                            -start => $previous_location->[1]+1,  
                                                            -end => $location->[0]-1, 
                                                            -frame => ".", 
                                                            -strand => "." , 
                                                            -tag => { 'ID' => $uniq_ID  }) ;
    # add feature in omniscient
    $omniscient->{"level1"}{"intergenic_region"}{lc($uniq_ID)}=$intergenic_feature; 
    
    # update local variables
    $intergenic_added++;
    $intergenicID++;
    $previous_location=$location;
  }     
}

# print result
print_omniscient( {omniscient => $omniscient, output => $gffout} );

dual_print1 "$intergenic_added intergenic_region added!\nBye Bye\n";
      #########################
      ######### END ###########
      #########################

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

__END__

=head1 NAME

agat_sp_add_intergenic_regions.pl

=head1 DESCRIPTION

The script aims to add intergenic features (intergenic_region) to gtf/gff file.
The intergenic regions are deduced from gene features (feature type gene from the 3rd column).

=head1 SYNOPSIS

    agat_sp_add_intergenic_regions.pl --gff infile --out outFile
    agat_sp_add_intergenic_regions.pl --help

=head1 OPTIONS

=over 8

=item B<--gff>, B<-f> or B<--ref> <file>

Input GTF/GFF file.

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
