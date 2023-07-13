#!/usr/bin/perl -w

package AGAT::AGAT;

use strict;
use warnings;
use Exporter;

use AGAT::OmniscientI;
use AGAT::OmniscientO;
use AGAT::OmniscientTool;
use AGAT::Config;
use AGAT::Levels;
use AGAT::OmniscientStat;
use AGAT::Utilities;
use AGAT::PlotR;
use Bio::Tools::GFF;

our $VERSION     = "v1.1.0";
our @ISA         = qw(Exporter);
our @EXPORT      = qw(get_agat_header print_agat_version get_agat_config handle_levels);
sub import {
  AGAT::AGAT->export_to_level(1, @_); # to be able to load the EXPORT functions when direct call; (normal case)
}

=head1 SYNOPSIS

  Meta package for conveniency. It allows to call all packages needed in once to deal with Omniscient data structure:
  $omniscient{'other'){'header'}[value, value]
  $omniscient->{"level1"}{$primary_tag}{$id}=$feature;
  $omniscient->{"level2"}{$primary_tag}{$parent} = [$feature,$feature];
  $omniscient->{"level3"}{$primary_tag}{$parent} = [$feature,$feature,$feature];

	It also contains function to deal with general features e.g. configuration file,
	AGAT version, AGAT header...

=head1 DESCRIPTION

    Omniscient packages are non-OO packages use to handle any kind of gtf/gff data.

=head1 AUTHOR

    Jacques Dainat - jacques.dainat@nbis.se

=cut

# ==============================================================================
#                          			=== MAIN ===

# Provide version
sub print_agat_version{
	print $VERSION."\n";
}

# Provide header
sub get_agat_header{

  my $header = qq{
 ------------------------------------------------------------------------------
|   Another GFF Analysis Toolkit (AGAT) - Version: $VERSION                      |
|   https://github.com/NBISweden/AGAT                                          |
|   National Bioinformatics Infrastructure Sweden (NBIS) - www.nbis.se         |
 ------------------------------------------------------------------------------
};

	return $header;
}

# Provide AGAT information
sub print_agat_info{
	my $header = get_agat_header();
	print <<MESSAGE
	$header 
AGAT checks, fixes, pads missing information (features/attributes) of any
kind of GTF/GFF (GXF) files and create complete, sorted and standardised 
GFF/GTF formated files. Over the years it has been enriched by many many
tools to perform just about any tasks that is possible related to GTF/GFF
format files (sanitizing, conversions, merging, modifying, filtering, FASTA
sequence extraction, adding information, etc). Comparing to other methods
AGAT is robust to even the most despicable GTF/GFF files.

By default AGAT automatically selects the appropriate parser and generates
a GFF3 output. This can be tuned via the config file.

AGAT contains 2 types of scripts: 
=================================

1) _sp_ prefix (slurp)
---------------
Data is loaded into memory via the AGAT parser that removes duplicate features, 
fixes duplicated IDs, adds missing ID and/or Parent attributes, deflates factorized 
attributes (attributes with several parents are duplicated with uniq ID), add missing 
features when possible (e.g. add exon if only CDS described, add UTR if CDS and exon 
described), fix feature locations (e.g. check exon is embedded in the parent features 
mRNA, gene), etc.
The AGAT parser defines relationship between features using 3 levels.
(e.g Level1=gene; Level2=mRNA,tRNA; Level3=exon,cds,utr)
The feature type information is stored within the 3rd column of a GXF file.
The parser needs to know to which level a feature type is part of. This information
is stored by default in a yaml file provided with the tool. We have implemented the
most common feature types met in gff/gtf files. If a feature type is not yet handle
by the parser it will throw a warning. You can easily inform the parser how
to handle it (level1, level2 or level3) by modifying the feature_levels.yaml file.

To access the AGAT feature_levels file: agat levels --expose

The  yaml file will appear in the working folder. By default, AGAT uses the
feature_levels.yaml file from the working directory when any.

AGAT parser phylosophy:
 a) Parse by Parent/child relationship
          or gene_id/transcript_id
   b) ELSE Parse by a comon tag (an attribute value shared by feature that must be grouped together.
           By default we are using locus_tag and gene_id as locus tag, but you can specify the one of your choice vi the config.
     c) ELSE Parse sequentially (features are grouped in a bucket, and the bucket change at each level2 feature met, and bucket(s) are linked to the first l1 top feature met)


2) _sq_ prefix (sequential):
---------------------
The gff file is read and processed from its top to the end line by line via the bioperl parser.
This is memory efficient, but no sanity check will be performed by the AGAT parser.

Configuration
=============

AGAT has a configuration file: agat_config.yaml

To access the AGAT config file: agat config --expose

The config yaml file will appear in the working folder. By default, AGAT 
uses the config file from the working directory when any.
The configuration can be used to change output format, to merge loci,
to activate tabix output, etc. (For _sq_ scripts only input/output format 
configuration parameters are used).
MESSAGE
}

# load configuration file from local file if any either the one shipped with AGAT
# return a hash containing the configuration.
sub get_agat_config{
	my ($args)=@_;

	my ($config_file_provided);
	if( defined($args->{config_file_in}) ) { $config_file_provided = $args->{config_file_in};}

	my $config_file_checked = get_config({type => "local", config_file_in => $config_file_provided}); #try local first, if none will take the original
	my $config = load_config({ config_file => $config_file_checked});
	check_config({ config => $config});
	return $config;
}

# ==============================================================================
#										=== fonction for agat caller ====

# $general is a hash reference to the overall application
# $config  is a hash reference with options
# $args    is an array reference with "residual" cmd line arguments
sub handle_main {
		my ($general, $config, $args) = @_;

		my $version = $general->{configs}[-1]{version};
		my $tools = $general->{configs}[-1]{tools};
		my $help = $general->{configs}[-1]{help};
		my $info = $general->{configs}[-1]{info};

		if($version){
			print_agat_version();
		}

		if($info){
			print_agat_info();
		}

		if($tools){
			my ($package, $filename, $line) = caller;
			my $agat_bin = $ENV{'AGAT_BIN'};
			opendir my $dir, $agat_bin or die "Cannot open directory: $!";
			my @files = readdir $dir;
			closedir $dir;
			foreach my $file (sort @files){
				# only file starting by agat_
				if ( $file =~ /^agat_/){
					print $file."\n";
				}
			}
		}

		# if help was called (or not arg provided) we let AppEaser continue to print help
		my $nb_args = keys %{$general->{configs}[-1]};
		if(! $help and $nb_args != 0){ exit 0;}
}

# Function to manipulate levels from the agat caller
sub handle_levels {
		my ($general, $config, $args) = @_;

		my $expose = $general->{configs}[-1]{expose};
		my $help = $general->{configs}[-1]{help};

		# Deal with Expose feature OPTION
		if($expose){
			expose_levels();
			print "Feature_levels YAML file copied in your working directory\n";
		}

		# if help was called (or not arg provided) we let AppEaser continue to print help
		my $nb_args = keys %{$general->{configs}[-1]};
		if(! $help and $nb_args != 0){ exit 0;}

}

# Function to manipulate config from the agat caller
sub handle_config {
		my ($general, $config, $args) = @_;

		my $expose = $general->{configs}[-1]{expose};
		my $help = $general->{configs}[-1]{help};

		# config file option
		my $verbose = $general->{configs}[-1]{verbose};
		my $progress_bar = $general->{configs}[-1]{progress_bar};
		my $config_new_name = $general->{configs}[-1]{output};
		my $log = $general->{configs}[-1]{log};
		my $debug = $general->{configs}[-1]{debug};
		my $tabix = $general->{configs}[-1]{tabix};
		my $merge_loci = $general->{configs}[-1]{merge_loci};
		my $throw_fasta = $general->{configs}[-1]{throw_fasta};
		my $force_gff_input_version = $general->{configs}[-1]{force_gff_input_version};
		my $output_format = $general->{configs}[-1]{output_format};
		my $gff_output_version = $general->{configs}[-1]{gff_output_version};
		my $gtf_output_version = $general->{configs}[-1]{gtf_output_version};
		my $create_l3_for_l2_orphan = $general->{configs}[-1]{create_l3_for_l2_orphan};
		my $locus_tag = $general->{configs}[-1]{locus_tag};
		my $check_sequential = $general->{configs}[-1]{check_sequential};
		my $check_l2_linked_to_l3 = $general->{configs}[-1]{check_l2_linked_to_l3};
		my $check_l1_linked_to_l2 = $general->{configs}[-1]{check_l1_linked_to_l2};
		my $remove_orphan_l1 = $general->{configs}[-1]{remove_orphan_l1};
		my $check_all_level3_locations = $general->{configs}[-1]{check_all_level3_locations};
		my $check_cds = $general->{configs}[-1]{check_cds};
		my $check_exons = $general->{configs}[-1]{check_exons};
		my $check_utrs = $general->{configs}[-1]{check_utrs};
		my $check_all_level2_locations = $general->{configs}[-1]{check_all_level2_locations};
		my $check_all_level1_locations = $general->{configs}[-1]{check_all_level1_locations};
		my $check_identical_isoforms = $general->{configs}[-1]{check_identical_isoforms};
		my $prefix_new_id = $general->{configs}[-1]{prefix_new_id};

		# Deal with Expose feature OPTION
		if($expose){
			my $config_file = get_config({type => "original"});
			my $config = load_config({ config_file => $config_file});
			print "Config loaded\n";

			# set config params on the fly
			my $modified_on_the_fly = undef;

			# integer 0-4
			if( defined($verbose) ){
				$config->{ verbose } = $verbose;
				$modified_on_the_fly = 1;
			}
			# bolean
			if( defined($progress_bar) ){
				$config->{ progress_bar } = _make_bolean($progress_bar);
				$modified_on_the_fly = 1;
			}
			# bolean
			if( defined($log) ){
				$config->{ log } = _make_bolean($log);
				$modified_on_the_fly = 1;
			}
			# bolean
			if( defined($debug) ){
				$config->{ debug } = _make_bolean($debug);
				$modified_on_the_fly = 1;
			}
			# bolean
			if( defined($tabix) ){
				$config->{ tabix } = _make_bolean($tabix);
				$modified_on_the_fly = 1;
			}
			# bolean
			if( defined($merge_loci) ){
				$config->{ merge_loci } = _make_bolean($merge_loci);
				$modified_on_the_fly = 1;
			}
			# bolean
			if( defined($throw_fasta) ){
				$config->{ throw_fasta } = _make_bolean($throw_fasta);
				$modified_on_the_fly = 1;
			}
			# integer
			if( defined($force_gff_input_version) ){
				$config->{ force_gff_input_version } = $force_gff_input_version;
				$modified_on_the_fly = 1;
			}
			# string
			if( defined($output_format) ){
				$config->{ output_format } = lc($output_format);
				$modified_on_the_fly = 1;
			}
			# 
			# integer
			if( defined($gff_output_version) ){
				$config->{ gff_output_version } = $gff_output_version;
				$modified_on_the_fly = 1;
			}
			# string
			if( defined($gtf_output_version) ){
				$config->{ gtf_output_version } = lc($gtf_output_version);
				$modified_on_the_fly = 1;
			}
			# bolean
			if( defined($create_l3_for_l2_orphan) ){
				$config->{ create_l3_for_l2_orphan } = _make_bolean($create_l3_for_l2_orphan);
				$modified_on_the_fly = 1;
			}
			# string
			if( defined($locus_tag) ){
				my @list = split(/,/, $locus_tag);
				$config->{ locus_tag } = \@list;
				$modified_on_the_fly = 1;
			}
			# bolean
			if( defined($check_sequential) ){
				$config->{ check_sequential } = _make_bolean($check_sequential);
				$modified_on_the_fly = 1;
			}
			# bolean
			if( defined($check_l2_linked_to_l3) ){
				$config->{ check_l2_linked_to_l3 } = _make_bolean($check_l2_linked_to_l3);
				$modified_on_the_fly = 1;
			}
			# bolean
			if( defined($check_l1_linked_to_l2) ){
				$config->{ check_l1_linked_to_l2 } = _make_bolean($check_l1_linked_to_l2);
				$modified_on_the_fly = 1;
			}
			# bolean
			if( defined($remove_orphan_l1) ){
				$config->{ remove_orphan_l1 } = _make_bolean($remove_orphan_l1);
				$modified_on_the_fly = 1;
			}
			# bolean
			if( defined($check_all_level3_locations) ){
				$config->{ check_all_level3_locations } = _make_bolean($check_all_level3_locations);
				$modified_on_the_fly = 1;
			}
			# bolean
			if( defined($check_cds) ){
				$config->{ check_cds } = _make_bolean($check_cds);
				$modified_on_the_fly = 1;
			}
			# bolean
			if( defined($check_exons) ){
				$config->{ check_exons } = _make_bolean($check_exons);
				$modified_on_the_fly = 1;
			}
			# bolean
			if( defined($check_utrs) ){
				$config->{ check_utrs } = _make_bolean($check_utrs);
				$modified_on_the_fly = 1;
			}
			# bolean
			if( defined($check_all_level2_locations) ){
				$config->{ check_all_level2_locations } = _make_bolean($check_all_level2_locations);
				$modified_on_the_fly = 1;
			}
			# bolean
			if( defined($check_all_level1_locations) ){
				$config->{ check_all_level1_locations } = _make_bolean($check_all_level1_locations);
				$modified_on_the_fly = 1;
			}
			# bolean
			if( defined($check_identical_isoforms) ){
				$config->{ check_identical_isoforms } = _make_bolean($check_identical_isoforms);
				$modified_on_the_fly = 1;
			}
			# string
			if( defined($prefix_new_id) ){
				$config->{ prefix_new_id } = $prefix_new_id;
				$modified_on_the_fly = 1;
			}

			if ($modified_on_the_fly) {
					print "Config modified\n";
			}

			# check config
			check_config({ config => $config});
			print "Config checked\n";

			 
			
			if ($modified_on_the_fly) {
				expose_config_hash({ config_in => $config, config_file_out => $config_new_name})
			} else {
				expose_config_file({config_file_in => $config_file, config_file_out => $config_new_name});
			}

			# inform user
			my $config_file_used;
			if($config_new_name){
				$config_file_used = $config_new_name;
			} else { $config_file_used = "agat_config.yaml"; }
			print "Config file written in your working directory ($config_file_used)\n";
		}

		# if help was called (or not arg provided) we let AppEaser continue to print help
		my $nb_args = keys %{$general->{configs}[-1]};
		if(! $help and $nb_args != 0){ exit 0;}

}

# transform 0 into false and 1 into true
sub _make_bolean{
	my ($value) = @_;

	my $result="false";
	if($value and $value ne "false"){
		$result="true";
	}
	return $result;
}

1;
