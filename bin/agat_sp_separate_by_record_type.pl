#!/usr/bin/env perl

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use AGAT::AGAT;
use Sort::Naturally;

start_script();
my $header = get_agat_header();
# ------------------------------- LOAD OPTIONS --------------------------------
my $opt_gfffile;
my $opt_output;
my $opt_help = 0;

# OPTION MANAGEMENT: split shared vs script options
my ($shared_argv, $script_argv) = split_argv_shared_vs_script(\@ARGV);

# Parse script-specific options
my $script_parser = Getopt::Long::Parser->new;
$script_parser->configure('bundling','no_auto_abbrev');
if ( ! $script_parser->getoptionsfromarray(
		$script_argv,
		'g|gff=s'         => \$opt_gfffile,
		'o|output=s'      => \$opt_output,
		'h|help!'         => \$opt_help,
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

if (! defined($opt_gfffile) ){
    pod2usage( {
           -message => "$header\nAt least 1 parameter is mandatory:\nInput reference gff file (-g).\n\n".
           "Ouptut is optional. Look at the help documentation to know more.\n",
           -verbose => 0,
           -exitval => 1 } );
}

# --- Manage config ---
my $shared_opts = parse_shared_options($shared_argv);
initialize_agat({ config_file_in => ( $shared_opts->{config} ), input => $opt_gfffile, shared_opts => $shared_opts });

######################
# Manage output file #

if (! $opt_output) {
	dual_print1 "Default output name: split_result\n";
	$opt_output="split_result";
}

if (-d $opt_output){
	die "The output directory choosen already exists. Please give me another Name.\n";
}
mkdir $opt_output;

                #####################
                #     MAIN          #
                #####################

######################
### Parse GFF input #
my ($hash_omniscient) = slurp_gff3_file_JD({ input => $opt_gfffile });

my $topfeatures = get_feature_type_by_agat_value($hash_omniscient, 'level1', 'topfeature');
my $standalones = get_feature_type_by_agat_value($hash_omniscient, 'level1', 'standalone');
my %handlers;
my $gffout;
#################
# == LEVEL 1 == #
foreach my $tag_l1 ( sort { ncmp ($a, $b) } keys %{$hash_omniscient->{'level1'}}){ # primary_tag_key_level1 = gene or repeat etc...
	# deal with topfeatures and standalone feature type
	if ( exists_keys ($topfeatures, ($tag_l1) ) or  exists_keys ($standalones, ($tag_l1) ) ){

		$gffout = prepare_gffout( $opt_output."/".$tag_l1.".gff");

		foreach my $key_l1 (keys %{$hash_omniscient->{'level1'}{$tag_l1}}){
			$gffout->write_feature($hash_omniscient->{'level1'}{$tag_l1}{$key_l1});
		}
		$gffout->close();
	}
	# deal with everything that is not topfeatures or standalone feature
	else{
		foreach my $key_l1 (keys %{$hash_omniscient->{'level1'}{$tag_l1}}){

	    #################
	    # == LEVEL 2 == #
	    my $level1_printed=undef;
	    foreach my $tag_l2 ( sort {$a cmp $b} keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...

	      if ( exists_keys ($hash_omniscient, ('level2', $tag_l2, $key_l1) ) ){
	        foreach my $feature_level2 ( @{$hash_omniscient->{'level2'}{$tag_l2}{$key_l1}}) {
	          #manage handler
	          if(! exists_keys ( \%handlers, ($tag_l2) ) ) {

	          $gffout = prepare_gffout( $opt_output."/".$tag_l2.".gff");

	            $handlers{$tag_l2}=$gffout;
	          }
	          $gffout = $handlers{$tag_l2};

	          #################
	          # == LEVEL 1 == #
	          if(! $level1_printed){
	            $gffout->write_feature($hash_omniscient->{'level1'}{$tag_l1}{$key_l1}); # print feature
	            $level1_printed=1;
	          }

	          #################
	          # == LEVEL 2 == #
	          $gffout->write_feature($feature_level2);

	          #################
	          # == LEVEL 3 == #
	          my $level2_ID = lc($feature_level2->_tag_value('ID'));

	          ###########
	          # Before tss
	          if ( exists_keys($hash_omniscient,('level3','tss',$level2_ID)) ){
	            foreach my $feature_level3 ( @{$hash_omniscient->{'level3'}{'tss'}{$level2_ID}}) {
	              $gffout->write_feature($feature_level3);
	            }
	          }

	          ######
	          # FIRST EXON
	          if ( exists_keys($hash_omniscient,('level3','exon',$level2_ID)) ){
	            foreach my $feature_level3 ( @{$hash_omniscient->{'level3'}{'exon'}{$level2_ID}}) {
	              $gffout->write_feature($feature_level3);
	            }
	          }
	          ###########
	          # SECOND CDS
	          if ( exists_keys($hash_omniscient,('level3','cds',$level2_ID)) ){
	            foreach my $feature_level3 ( @{$hash_omniscient->{'level3'}{'cds'}{$level2_ID}}) {
	              $gffout->write_feature($feature_level3);
	            }
	          }

	          ###########
	          # Last tts
	          if ( exists_keys($hash_omniscient,('level3','tts',$level2_ID)) ){
	            foreach my $feature_level3 ( @{$hash_omniscient->{'level3'}{'tts'}{$level2_ID}}) {
	              $gffout->write_feature($feature_level3);
	            }
	          }

	          ###########
	          # The rest
	          foreach my $primary_tag_key_level3 (keys %{$hash_omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...
	            if( ($primary_tag_key_level3 ne 'cds') and ($primary_tag_key_level3 ne 'exon') and ($primary_tag_key_level3 ne 'tss') and ($primary_tag_key_level3 ne 'tts')){
	              if ( exists_keys ($hash_omniscient, ('level3', $primary_tag_key_level3, $level2_ID) ) ){
	                foreach my $feature_level3 ( @{$hash_omniscient->{'level3'}{$primary_tag_key_level3}{$level2_ID}}) {
	                  $gffout->write_feature($feature_level3);
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

#Close all FH opened
foreach my $key (keys %handlers){
	$handlers{$key}->close();
}

# --- final messages ---
end_script();

__END__

=head1 NAME

agat_sp_separate_by_record_type.pl

=head1 DESCRIPTION

The script will separate the features from the gff input file into different files according to
the record type. A record represent all features linked collectively by Parent/ID relationships.
(e.g gene + mrna + exon + cds + utr of a locus).

a) When the record contains Level2 feature, the record type is the Level2 feature type (e.g tRNA,mRNA,ncRNA etc...)
b) Some features do not have children (top and standalone level1 features) e.g. location,region,chromosome.
In such case the record type is the level1 feature type.

=head1 SYNOPSIS

    agat_sp_separate_by_record_type.pl -g infile.gff [ -o outfolder ]
    agat_sp_separate_by_record_type.pl --help

=head1 OPTIONS

=over 8

=item B<-g>, B<--gff> or B<-ref>

Input GTF/GFF file.

=item B<-o> or B<--output>

Output folder.  If no output folder provided, the default name will be <split_result>.

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
