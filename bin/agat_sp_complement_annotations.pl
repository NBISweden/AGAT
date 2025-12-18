#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long;
use Pod::Usage;
use List::MoreUtils qw(uniq);
use AGAT::AGAT;

start_script();
my $header = get_agat_header();
# -----------------------------------------------------------------------------------------------
my $config;
my $opt_output = undef;
my @opt_files;
my $ref = undef;
my $size_min = 0;
my $opt_help;

# OPTION MANAGMENT
my @copyARGV=@ARGV;
my ($shared_argv, $script_argv) = split_argv_shared_vs_script(\@ARGV);
my $script_parser = Getopt::Long::Parser->new;
$script_parser->configure('bundling','no_auto_abbrev');
if ( !$script_parser->getoptionsfromarray(
  $script_argv,
  "h|help"                   => \$opt_help,
  "ref|r|i=s"                => \$ref,
  "add|a=s"                  => \@opt_files,
  "size_min|s=i"             => \$size_min,
  "output|out|o=s"   => \$opt_output))
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

if (! $ref or ! @opt_files ){
    pod2usage( {
           -message => "$header\nAt least 2 files are mandatory:\n --ref file1 --add file2\n\n",
           -verbose => 0,
           -exitval => 2 } );
}

my ($shared_opts) = parse_shared_options($shared_argv);
initialize_agat({ config_file_in => $shared_opts->{config}, input => $ref, shared_opts => $shared_opts });
# -----------------------------------------------------------------------------------------------

######################
# Manage output file #
my $gffout = prepare_gffout( $opt_output );

                #####################
                #     MAIN          #
                #####################

######################
### Parse GFF input #

my ($hash_omniscient) = slurp_gff3_file_JD({ input => $ref });

info_omniscient($hash_omniscient);

#Add the features of the other file in the first omniscient. It takes care of name to not have duplicates
foreach my $next_file (@opt_files){
  my ($hash_omniscient2) = slurp_gff3_file_JD({ input => $next_file });
  dual_print1 "$next_file GFF3 file parsed\n";
  info_omniscient($hash_omniscient2);

  # Quick stat hash before complement
  my %quick_stat1;
  foreach my $level ( ('level1', 'level2') ){
    foreach  my $tag (keys %{$hash_omniscient->{$level}}) {
      my $nb_tag = keys %{$hash_omniscient->{$level}{$tag}};
      $quick_stat1{$level}{$tag} = $nb_tag;
    }
  }
  
  ####### COMPLEMENT #######
  complement_omniscients($hash_omniscient, $hash_omniscient2, $size_min); # deal with identical ID by renaming them
  dual_print1 "\nComplement done !\n";


 #RESUME COMPLEMENT
  my $complemented=undef;
  # Quick stat hash after complement
  my %quick_stat2;
  foreach my $level ( ('level1', 'level2') ){
    foreach  my $tag (keys %{$hash_omniscient->{$level}}) {
      my $nb_tag = keys %{$hash_omniscient->{$level}{$tag}};
      $quick_stat2{$level}{$tag} = $nb_tag;
    }
  }

  #About tag from hash1 added which exist in hash2
  foreach my $level ( ('level1', 'level2') ){
    foreach my $tag (keys %{$quick_stat1{$level}}){
      if ($quick_stat1{$level}{$tag} != $quick_stat2{$level}{$tag} ){
        dual_print1 "We added ".($quick_stat2{$level}{$tag}-$quick_stat1{$level}{$tag})." $tag(s)\n";
        $complemented=1;
      }
    }
  }
  #About tag from hash2 added which dont exist in hash1
  foreach my $level ( ('level1', 'level2') ){
    foreach my $tag (keys %{$quick_stat2{$level}}){
      if (! exists $quick_stat1{$level}{$tag} ){
        dual_print1 "We added ".$quick_stat2{$level}{$tag}." $tag(s)\n";
        $complemented=1;
      }
    }
  }
  #If nothing added
  if(! $complemented){
    dual_print1 "\nNothing has been added\n";
  }
  else{
    dual_print1 "\nNow the data contains:\n";
    info_omniscient($hash_omniscient);
  }
}

########
# Print results
print_omniscient( {omniscient => $hash_omniscient, output => $gffout} );

# --- final messages ---
end_script();

# -----------------------------------------------------------------------------------------------
__END__

=head1 NAME

agat_sp_complement_annotations.pl

=head1 DESCRIPTION

The script allows to complement a reference annotation with other annotations.
A l1 feature from the addfile.gff that does not overlap a l1 feature from the reference annotation will be added.
A l1 feature from the addfile.gff without a CDS that overlaps a l1 feature with a CDS from the reference annotation will be added.
A l1 feature from the addfile.gff with a CDS that overlaps a l1 feature without a CDS from the reference annotation will be added.
A l1 feature from the addfile.gff with a CDS that overlaps a l1 feature with a CDS from the reference annotation will be added only if the CDSs don't overlap.
A l1 feature from the addfile.gff without a CDS that overlaps a l1 feature without a CDS from the reference annotation will be added only if none of the l3 features overlap.
/!\ It is sufficiant that only one isoform is overlapping to prevent the whole gene (l1 feature) from the addfile.gff to be added in the output.

=head1 SYNOPSIS

    agat_sp_complement_annotations.pl --ref annotation_ref.gff --add addfile1.gff --add addfile2.gff --out outFile
    agat_sp_complement_annotations.pl --help

=head1 OPTIONS

=over 8

=item B<--ref>,  B<-r> or B<-i> <file>

Input GTF/GFF file used as reference.

=item B<--add> or B<-a> <file>

Annotation(s) file you would like to use to complement the reference annotation. You can specify as much file you want like so: -a addfile1 -a addfile2 -a addfile3
/!\ The order you provide these files matter. Once the reference file has been complemented by file1, this new annotation becomes the new reference that will be complemented by file2 etc.
/!\ The result with -a addfile1 -a addfile2 will differ to the result from -a addfile2 -a addfile1. So, be aware of what you want if you use several addfiles.

=item  B<--size_min> or B<-s> <int>

Option to keep the non-overlping gene only if the CDS size (in nucleotide) is over the minimum size defined. Default = 0 that means all of them are kept.

=item B<--output>, B<--out> or B<-o> <file>

Output file to create (default GFF3 - see config to modify output format).
If no output file is specified, the output will be written to STDOUT.
It will contain the reference annotation with all the non-overlapping newly added genes from addfiles.gff.

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
