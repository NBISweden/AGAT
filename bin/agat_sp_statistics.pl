#!/usr/bin/env perl

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use AGAT::AGAT;

start_script();
my $header = get_agat_header();

# ------------------------------- LOAD OPTIONS --------------------------------
my $gff = undef;
my $opt_output = undef;
my $opt_yaml = undef;
my $opt_percentile = 90;
my $opt_genomeSize = undef;
my $opt_plot = undef;
my $opt_raw = undef;
my $opt_help = 0;

# OPTION MANAGEMENT: split shared vs script options
my ($shared_argv, $script_argv) = split_argv_shared_vs_script(\@ARGV);

# Parse script-specific options
my $script_parser = Getopt::Long::Parser->new;
$script_parser->configure('bundling','no_auto_abbrev');
if ( ! $script_parser->getoptionsfromarray(
    $script_argv,
    'h|help!'       => \$opt_help,
    'o|out|output=s'    => \$opt_output,
    'percentile=i'  => \$opt_percentile,
    'yaml!'         => \$opt_yaml,
    'r|raw!'        => \$opt_raw,
    'd|p!'          => \$opt_plot,
    'g|f|gs=s'      => \$opt_genomeSize,
    'gff|i=s'       => \$gff,
  ) )
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
       -exitval => 1 } );
}

# --- Manage config ---
my $shared_opts = parse_shared_options($shared_argv);
initialize_agat({ config_file_in => ( $shared_opts->{config} ), input => $gff, shared_opts => $shared_opts });

#### IN / OUT
my $out = prepare_fileout($opt_output);
if(defined($opt_yaml)){
  if( defined($opt_output)){
    $opt_yaml = $opt_output.".yaml";
  }
}

# Manage raw data
if($opt_raw){
  if ($opt_output){
    $opt_raw = $opt_output."_raw_data";
  }
  else{
    $opt_raw = "raw_data";

    if (-f $opt_raw){
      die "Cannot create a directory with the name $opt_raw because a file with this name already exists.\n";
    }
    if (-d $opt_raw){
      die "Cannot create a directory with the name $opt_raw because a folder with this name already exists.\n";
    }
  }
}

#Manage plot folder output
if($opt_plot){

  # Check if dependencies for plot are available
  if ( ! may_i_plot() ) {
    $opt_plot = undef;
  }
  else{
    if ($opt_output){
      $opt_plot = $opt_output."_distribution_plots";
    }
    else{
      $opt_plot = "distribution_plots";

      if (-f $opt_plot){
        die "Cannot create a directory with the name $opt_plot because a file with this name already exists.\n";
      }
      if (-d $opt_plot){
        die "The default output directory $opt_plot use to save the distribution plots already exists. Please give me another folder name.\n";
      }
    }
  }
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>     MAIN     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

######################
### Parse GFF input #
my ($hash_omniscient) =  slurp_gff3_file_JD({ input => $gff });
### END Parse GFF input #
#########################

##############
# STATISTICS #
dual_print1 "Compute statistics\n";
print_omniscient_statistics ({ input   => $hash_omniscient,
                                 genome    => $opt_genomeSize,
                                 percentile=> $opt_percentile,
                                 output    => $out,
                                 yaml      => $opt_yaml,
                                 raw       => $opt_raw,
                                 distri    => $opt_plot,
                                 isoform   => 1,
                                 verbose   => $CONFIG->{verbose}
                               });
# END STATISTICS #
##################

# --- final messages ---
end_script();

__END__

=head1 NAME

agat_sp_statistics.pl

=head1 DESCRIPTION

The script provides exhaustive statistics of a gft/gff file.
/!\ If you have isoforms in your file, even if correct, some values calculated
might sounds incoherent: e.g. total length mRNA can be superior than the genome size.
Because all isoforms length is added... It is why by default
we always compute the statistics twice when there are isoforms, once with the
isoforms, once without (In that case we keep the longest isoform per locus).

=head1 SYNOPSIS

    agat_sp_statistics.pl --gff file.gff  [ -o outfile ]
    agat_sp_statistics.pl --help

=head1 OPTIONS

=over 8

=item B<--gff> or B<-i> <file>

Input GTF/GFF file.

=item B<--gs>, B<-f> or B<-g> <file or int>

This option inform about the genome size in oder to compute more statistics. You can give the size in Nucleotide or directly the fasta file.

=item B<-d> or B<-p>

When this option is used, an histogram of distribution of the features will be printed in pdf files in a folder with distribution_plots suffix. (d means distribution, p means plot).

=item  B<-o>, B<--out> or B<--output> <file>

File where will be written the result. If no output file is specified, the output will be written to STDOUT.

=item B<--percentile> <int>

Percentile to compute. Default is 90.

=item B<-r> or B<--raw>

When this option is used, the raw data (same as used to create histogram of distribution of the features) are printed in a dedicated folder with raw_data suffix.

=item B<--yaml>

When this option is activated, a second output will be printed either in STDOUT if no output provided or in <output.yaml> (a .yaml suffix is added to the --output value provided).

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
