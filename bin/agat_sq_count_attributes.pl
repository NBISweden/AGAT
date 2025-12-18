#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Clone 'clone';
use Getopt::Long;
use Pod::Usage;
use List::MoreUtils qw(uniq);
use AGAT::AGAT;
start_script();

my $header = get_agat_header();
# ---------------------------- OPTIONS ----------------------------
my $gff = undef;
my $opt_help= 0;
my $attribute='transcript_id';
my $outfile=undef;
my $cpt_case=0;

# OPTION MANAGEMENT: partition @ARGV into shared vs script options via library
my ($shared_argv, $script_argv) = split_argv_shared_vs_script(\@ARGV);

# Parse script-specific options from its own list
my $script_parser = Getopt::Long::Parser->new;
$script_parser->configure('bundling','no_auto_abbrev');
if ( ! $script_parser->getoptionsfromarray(
  $script_argv,
  "h|help"      => \$opt_help,
  "gff|f=s"     => \$gff,
  "tag|att=s"   => \$attribute,
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

if ( ! $gff ){
    pod2usage( {
           -message => "$header\nAt least 1 parameter is mandatory:\nInput reference gff file (--gff)\n\n",
           -verbose => 0,
           -exitval => 2 } );
}

# Parse shared options (CPU, config, etc.)
my ($shared_opts) = parse_shared_options($shared_argv);

# --- Manage config ---
initialize_agat({ config_file_in => ( $shared_opts->{config} ), input => $gff, shared_opts => $shared_opts });

# ------------------------------------------------------------------------------
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>     MAIN     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

dual_print1 "Looking to $attribute attribute.\n";

# Manage input gff file
my $format = $CONFIG->{force_gff_input_version};
if(! $format ){ $format = select_gff_format($gff); }
my $inputfh = open_maybe_gz($gff);
my $ref_in = AGAT::BioperlGFF->new(-fh => $inputfh, -gff_version => $format);

# set progression bar
set_progression_counter( $gff);
my $line_cpt=0;

my %hash_values;
my $nb_attributes=0;
while (my $feature = $ref_in->next_feature() ) {
  $line_cpt++;

  if ($feature->has_tag($attribute) ){
    $hash_values{$feature->_tag_value($attribute)}++;
    $nb_attributes++;
  }


  #####################
  #Display progression
  update_progression_counter($line_cpt);
}

##Last round
my $result = scalar keys %hash_values;
dual_print1 "$line_cpt features read. Among them, $nb_attributes has the $attribute attribute.\n".
"There is $result unique value within $attribute attribute\n";

# --- final messages ---
end_script();

# ---------------------------- FUNCTIONS ----------------------------

__END__

=head1 NAME

agat_sq_count_attributes.pl

=head1 DESCRIPTION

The script count the number of a choosen attribute and also count the number of
unique value of this attribute.

=head1 SYNOPSIS

    agat_sq_count_attributes.pl --gff file.gff  --att gene_id [ -o outfile ]
    agat_sq_count_attributes.pl --help

=head1 OPTIONS

=over 8

=item B<--gff> or B<-f> <file>

Input GTF/GFF file.

=item B<--tag>, B<--att> <string>
The name of the attribute that will be investigated.

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
