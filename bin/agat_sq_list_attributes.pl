#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Clone 'clone';
use Getopt::Long;
use Pod::Usage;
use IO::File;
use List::MoreUtils qw(uniq);
use AGAT::AGAT;
start_script();
my $header = get_agat_header();
# ---------------------------- OPTIONS ----------------------------
my %handlers;
my $gff = undef;
my $opt_help= 0;
my $primaryTag=undef;
my $outfile=undef;

# OPTION MANAGEMENT: partition @ARGV into shared vs script options via library
my ($shared_argv, $script_argv) = split_argv_shared_vs_script(\@ARGV);

# Parse script-specific options from its own list
my $script_parser = Getopt::Long::Parser->new;
$script_parser->configure('bundling','no_auto_abbrev');
if ( ! $script_parser->getoptionsfromarray(
  $script_argv,
  "h|help"                 => \$opt_help,
  "gff|f=s"                => \$gff,
  "p|t|l=s"                => \$primaryTag,
  "output|out|o=s"         => \$outfile))

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
           -message => "$header\nAt least 1 parameter is mandatory:\nInput reference gff file (--gff) \n\n",
           -verbose => 0,
           -exitval => 2 } );
}

# Parse shared options (CPU, config, etc.)
my ($shared_opts) = parse_shared_options($shared_argv);

# --- Manage config ---
initialize_agat({ config_file_in => $shared_opts->{config}, input => $gff, shared_opts => $shared_opts });

# ------------------------------------------------------------------------------

# Manage $primaryTag
my @ptagList;
if(! $primaryTag or $primaryTag eq "all"){
  dual_print1 "We will work on attributes from all features\n";
  push(@ptagList, "all");
}
else{
   @ptagList= split(/,/, $primaryTag);
   foreach my $tag (@ptagList){
    dual_print1 "We will work on attributes from $tag feature.\n";
   }
}

# Manage input fasta file
my $format = $CONFIG->{force_gff_input_version};
if(! $format ){ $format = select_gff_format($gff); }
my $inputfh = open_maybe_gz($gff);
my $ref_in = AGAT::BioperlGFF->new(-fh => $inputfh, -gff_version => $format);


                #####################
                #     MAIN          #
                #####################

my %all_attributes;
my %attributes_per_level;
######################
### Parse GFF input #

# set progression bar
set_progression_counter( $gff);
my $line_cpt=0;

my $geneName=undef;
while (my $feature = $ref_in->next_feature() ) {
  $line_cpt++;


    manage_attributes($feature, \@ptagList, \%all_attributes, \%attributes_per_level);

  #Display progression
  update_progression_counter($line_cpt);
}

#print "We added $nbNameAdded Name attributes\n";

my $out = prepare_fileout($outfile);

# Print information by feature
my $nbFeat = scalar keys %attributes_per_level;
print $out "\nWe met ".$nbFeat." different feature types.";
foreach my $feature_type ( sort keys %attributes_per_level){
  my $nbAtt = scalar keys %{$attributes_per_level{$feature_type}};
  print $out "\nHere the list of all the attributes tags met for the feature type <".$feature_type."> (".$nbAtt." attributes):\n";
  foreach my $attribute ( sort keys %{$attributes_per_level{$feature_type}} ){
    print $out $attribute."\n";
  }
}

# Print Global information
my $nbAtt = scalar keys %all_attributes;
print $out "\nHere the list of all the attributes tags met (".$nbAtt." attributes):\n";
foreach my $attribute ( sort keys %all_attributes){
  print $out $attribute."\n";
}

# --- final messages ---
end_script();

# ---------------------------- FUNCTIONS ----------------------------

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

sub  manage_attributes{
  my  ($feature, $ptagList, $all_attributes, $attributes_per_level)=@_;

  my $primary_tag=$feature->primary_tag;

  # check primary tag (feature type) to handle
  foreach my $ptag (@$ptagList){

    if($ptag eq "all"){
      tag_from_list($feature,$all_attributes, $attributes_per_level);
    }
    elsif(lc($ptag) eq lc($primary_tag) ){
      tag_from_list($feature,$all_attributes, $attributes_per_level);
    }
  }
}

sub tag_from_list{
  my  ($feature, $all_attributes, $attributes_per_level)=@_;

  foreach my $tag ($feature->get_all_tags) {
      # create handler if needed (on the fly)
      if(! exists_keys( $all_attributes,($tag) ) ) {
        $all_attributes{$tag}++;
      }
      if(! exists_keys ( $attributes_per_level,($feature->primary_tag,$tag) ) ) {
        $attributes_per_level{$feature->primary_tag}{$tag}++;
      }

  }
}


__END__


=head1 NAME

agat_sq_list_attributes.pl

=head1 DESCRIPTION

The script take a gff3 file as input. -
The script give information about attribute tags used within you file.

=head1 SYNOPSIS

    agat_sq_list_attributes.pl -gff file.gff -p level2,cds,exon [ -o outfile ]
    agat_sq_list_attributes.pl --help

=head1 OPTIONS

=over 8

=item B<--gff> or B<-f> <file>

Input GTF/GFF file.

=item B<-p>,  B<-t> or  B<-l> <string>

primary tag option, case insensitive, list. Allow to specied the feature types that will be handled.
You can specified a specific feature by given its primary tag name (column 3) as: cds, Gene, MrNa
You can specify directly all the feature of a particular level:
      level2=mRNA,ncRNA,tRNA,etc
      level3=CDS,exon,UTR,etc
By default all feature are taking in account. fill the option by the value "all" will have the same behaviour.

=item  B<-o>, B<--out> or B<--output> <file>

Output file to create (default GFF3 - see config to modify output format).
If no output file is specified, the output will be written to STDOUT.

=item B<-h> or B<--help>

Display this helpful text.

=back

=head1 SHARED OPTIONS

Shared options are defined in the AGAT configuration file and can be overridden via the command line for this script only.
Common shared options are listed below; for the full list, please refer to the AGAT agat_config.yaml.
Note: For _sq_ scripts, only the following options are supported: verbose, output_format, gff_output_version, gtf_output_version, progress_bar, and tabix.

=over 8

=item B<--config> <file>

Path to a custom AGAT configuration file.  
By default, AGAT uses `agat_config.yaml` from the working directory if present, otherwise the default file shipped with AGAT
(available locally via `agat config --expose`).

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
