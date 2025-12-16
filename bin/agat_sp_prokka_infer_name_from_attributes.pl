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
# ------------------------------- LOAD OPTIONS --------------------------------

my $gff = undef;
my $opt_help= 0;
my $force=undef;
my $outfile=undef;

# OPTION MANAGEMENT: split shared vs script options
my ($shared_argv, $script_argv) = split_argv_shared_vs_script(\@ARGV);

# Parse script-specific options
my $script_parser = Getopt::Long::Parser->new;
$script_parser->configure('bundling','no_auto_abbrev');
if ( ! $script_parser->getoptionsfromarray(
  $script_argv,
  'h|help'                  => \$opt_help,
  'gff|f=s'                 => \$gff,
  'force!'                  => \$force,
  'output|outfile|out|o=s'  => \$outfile,
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

# Prepare output
my $gffout = prepare_gffout( $outfile);

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>     MAIN     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

######################
### Parse GFF input #
my ($hash_omniscient) = slurp_gff3_file_JD({ input => $gff });

my $nbNameAdded=0;

foreach my $tag (keys %{$hash_omniscient->{'level1'}}){
  foreach my $id (keys %{$hash_omniscient->{'level1'}{$tag}}){

    my $feature=$hash_omniscient->{'level1'}{$tag}{$id};

    if($feature->has_tag('name')){
      my $name=$feature->_tag_value('name');
      create_or_replace_tag($feature,'Name', $name);
      $feature->remove_tag('name');
    }

    #Name already contained in the gene attribute.
    if($feature->has_tag('gene')){
      # we get the Name
      my $name=$feature->_tag_value('gene');

      # If no attribute Name or if we have to replace it
      if(! $feature->has_tag('Name') or ($force)){
        create_or_replace_tag($feature,'Name', $name);
        $nbNameAdded++;
      }
      elsif($feature->has_tag('Name') and ( ! $force)){
        dual_print1 "Feature contains already an attribute Name. You can force it replacement by using the option --force\n";
      }
      dual_print1 "Name found in gene attribute = $name\n";
    }# Name not found in gene attribute. So we try to get the name included in the inference attribute.
    elsif($feature->has_tag('inference')){
      my @inferenceAtt=$feature->get_tag_values('inference');
      if ($#inferenceAtt > 0){

        my @tab = split /\|/,$inferenceAtt[$#inferenceAtt]; # split the last value by the character |
        my $name = $tab[$#tab];

        # SKIP case
        if($name =~ /protein motif:Pfam:/i){
          next;
        }
        if($name =~ /protein motif:CLUSTERS:/i){
          next;
        }
        if($name =~ /similar to AA sequence:UniProtKB:/i){
          next;
        }
        # ELSE name contains the Uniprot header of the protein coming from "--proteins" option ( Fasta file of trusted proteins to first annotate from ).


        if(! $feature->has_tag('Name') or ($force)){
          create_or_replace_tag($feature,'Name', $name);
          $nbNameAdded++;
        }
        elsif($feature->has_tag('Name') and ( ! $force)){
          dual_print1 "Feature contains already an attribute Name. You can force it replacement by using the option --force\n";
        }
        dual_print1 "My Name get in inference attribute = $name\n";
      }
    }
  }
}
dual_print1 "We added $nbNameAdded Name attributes\n";

print_omniscient( {omniscient => $hash_omniscient, output => $gffout} );

# --- final messages ---
end_script();

__END__

=head1 NAME

agat_sp_prokka_infer_name_from_attributes.pl

=head1 DESCRIPTION

The script aims to fill a Name attribute based on <gene> attribute in a prokka gff
annotation file. If no gene attribute is present it take if from the <inference>
attribute.

=head1 SYNOPSIS

    agat_sp_prokka_infer_name_from_attributes.pl -gff file.gff  [ -o outfile ]
    agat_sp_prokka_infer_name_from_attributes.pl --help

=head1 OPTIONS

=over 8

=item B<--gff> or B<-f>

Input GTF/GFF file.

=item B<--force>

If Name attribute already exists, they will be replaced if a new one is found

=item B<-o> , B<--output> , B<--out> or B<--outfile>

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
