#!/usr/bin/perl -w

package AGAT::Config;

use strict;
use warnings;
use YAML qw(DumpFile LoadFile);
use File::Copy;
use File::ShareDir ':ALL';
use AGAT::Utilities;
use Cwd qw(cwd);
use Exporter;


our @ISA = qw(Exporter);
our @EXPORT = qw( load_config expose_config_file check_config get_config expose_config_hash );
sub import {
  AGAT::Config->export_to_level(1, @_); # to be able to load the EXPORT functions when direct call; (normal case)
  AGAT::Config->export_to_level(2, @_); # to be able to load the EXPORT functions when called from one level up;
}

=head1 SYNOPSIS

This is the code to handle AGAT's config file. Accessing file is slow, the
values from the yaml files will be stored within OMNISCIENT when accessed the
first time.

=head1 DESCRIPTION

 This lib contains functions to parse YAML config files and store the values within
 OMNISCIENT and functions to assess the values.

=head1 AUTHOR

  Jacques Dainat - jacques.dainat@nbis.se

=cut

#	-----------------------------------CONSTANT-----------------------------------

my $config_file= ("agat_config.yaml");

#	------------------------------------GENERAL------------------------------------	 


# @Purpose: Load yaml file, check all is set, shift false to 0, return the config
# @input: 4 =>	verbose, config_file (path), log, debug
# @output: 1 => hash
# @Remark: none
sub load_config{
	my ($args) = @_;

	# -------------- INPUT --------------
	# -- Declare all variables and fill them --
	my ($verbose, $log, $debug, $path);

	if( ! defined($args->{verbose}) ) { $verbose = undef;} else{ $verbose = $args->{verbose}; }
	if( ! defined($args->{log}) ) { $log = undef;} else{ $log = $args->{log}; }
	if( ! defined($args->{debug}) ) { $debug = undef;} else{ $debug = $args->{debug}; }
	if( ! defined($args->{config_file}) ) { warn "No config file provided!"; exit;} else{ $path = $args->{config_file}; }

	my $config = LoadFile($path);

	# check if all needed parameters are present
	check_config({config => $config});

	# false does not exists in perl, must be replaced by 0
	foreach my $key (keys %{$config}){
		if ( lc($config->{$key}) eq "false" ){
			$config->{$key}=0;
		}
	}

	return $config;
}

# @Purpose: Select which config file to use (local or original shiped with AGAT)
# If type=original we only take the original config
# If type=local we try first to take the local one. If none we take the original one.
sub get_config{
	my ($args) = @_;

	# -------------- INPUT --------------
	# -- Declare all variables and fill them --
	my ( $verbose, $log, $debug, $type ) ;
	if( ! defined($args->{verbose}) ) { $verbose = undef;} else{ $verbose = $args->{verbose}; }
	if( ! defined($args->{log}) ) { $log = undef;} else{ $log = $args->{log}; }
	if( ! defined($args->{debug}) ) { $debug = undef;} else{ $debug = $args->{debug}; }
	if( ! defined($args->{type}) ) { $type = "local";} else{ $type = $args->{type};}
	if( !$type eq "local" and !$type eq "original" ){
		warn "type must be local or original. $type unknown!";exit;
	}

	#set run directory
	my $run_dir = cwd;

	my $path=undef;

	# get local config if any
	if ($type eq "local") {
		$path = $run_dir."/".$config_file;
		if (-e $path){
			dual_print($log, "Using $config_file file found in your working directory.\n", $verbose );
		} else {
			$path = undef;
		}
	}
	if (! $path) { #otherwise use the standard location ones
		$path = dist_file('AGAT', $config_file);
		dual_print($log, "Using standard $path file\n", $verbose );
	}
	return $path;
}

sub expose_config_hash{
	my ($args)=@_;

	my ($config);
	if( ! defined($args->{config}) ) { $config = undef;} else{ $config = $args->{config};}

	DumpFile("agat_config.yaml", $config);
}

# @Purpose: Write the config hash in a yaml file in the current directory 
sub expose_config_file{
	my ($args)=@_;

	my ($path);
	if( ! defined($args->{config_file}) ) { $path = undef;} else{ $path = $args->{config_file};}

	#set run directory
	if(! $path){
		$path = dist_file('AGAT', $config_file);
		print "Path where $config_file is standing according to dist_file: $path\n";
	}
	# copy the file locally
	my $run_dir = cwd;
	copy($path, $run_dir) or die print "Copy failed: $!";
}

# @Purpose: Check config value to be sure everything is set as expected
sub check_config{
	my ($args)=@_;

	my ($config);
	if( ! defined($args->{config}) ) { warn "Config file missing!\n";} else{ $config = $args->{config};}

	# check config parameters
	my $error;
	if( exists_keys($config,("verbose") ) ){
		if( $config->{ verbose } < 0 or $config->{ verbose } > 4){
			print "Verbose parameter must be between 0 and 4.\n";
			$error = 1;
		}
	}
	else{
		print "verbose parameter missing in the configuration file.\n";
		$error = 1;
	}
	if( !  exists_keys($config,("progress_bar") ) ){
		print "progress_bar parameter missing in the configuration file.\n";
		$error = 1;
	}
	if( !  exists_keys($config,("log") ) ){
		print "progress_bar parameter missing in the configuration file.\n";
		$error = 1;
	}
	if( !  exists_keys($config, ("debug") ) ){
		print "debug parameter missing in the configuration file.\n";
		$error = 1;
	}
	if( !  exists_keys($config, ("tabix") ) ){
		print "tabix parameter missing in the configuration file.\n";
		$error = 1;
	}
	if( !  exists_keys($config, ("merge_loci") ) ) {
		print "merge_loci parameter missing in the configuration file.\n";
		$error = 1;
	}
	if( !  exists_keys($config, ("throw_fasta") ) ) {
		print "throw_fasta parameter missing in the configuration file.\n";
		$error = 1;
	}
	if( exists_keys($config, ("force_gff_input_version") ) ) {
		my %values = (0 => 1, 1 => 1, 2 => 1, 2.5 => 1, 3 => 1);
		if (! exists_keys(\%values,($config->{ force_gff_input_version }) ) ) {
			print "force_gff_input_version parameter must be 0, 1, 2, 2.5 or 3.\n";
			$error = 1;
		}
	}
	else{
		print "force_gff_input_version parameter missing in the configuration file.\n";
		$error = 1;
	}

	if( exists_keys($config, ("output_format") ) ) {
		my %values = ( gff => 1, gtf => 1 ); # case was set to lowercase when loaded into hash
		$config->{ output_format } = lc($config->{ output_format }); # set it to lowercase
		if (! exists_keys(\%values,( $config->{ output_format } ) ) ) {
			print "output_format parameter must be GFF or GTF.\n";
			$error = 1;
		}
	}
	else{
		print "output_format parameter missing in the configuration file.\n";
		$error = 1;
	}

	if( exists_keys($config, ("gff_output_version") ) ) {
		my %values = ( 1 => 1, 2 => 1, 2.5 => 1, 3 => 1 );
		if (! exists_keys(\%values,($config->{ gff_output_version }) ) ) {
			print "gff_output_version parameter must be 1, 2, 2.5 or 3.\n";
			$error = 1;
		}
	}
	else{
		print "gff_output_version parameter missing in the configuration file.\n";
		$error = 1;
	}
	if( exists_keys($config, ("gtf_output_version") ) ) {
		my %values = ( 1 => 1, 2 => 1, 2.1 => 1, 2.2 => 1, 2.5 => 1, 3 => 1, relax => 1 );
		if (! exists_keys(\%values,($config->{ gtf_output_version }) ) ) {
			print "gtf_output_version parameter must be 1, 2, 2.1, 2.2, 2.5, 3 or relax.\n";
			$error = 1;
		}
	}
	else{
		print "gtf_output_version parameter missing in the configuration file.\n";
		$error = 1;
	}
	if( ! exists_keys($config, ("create_l3_for_l2_orphan") ) ) {
		print "create_level3_for_level2_orphan parameter missing in the configuration file.\n";
		$error = 1;
	}
	if( ! exists_keys($config, ("locus_tag") ) ) {
		print "locus_tag parameter missing in the configuration file.\n";
		$error = 1;
	}
	if( !  exists_keys($config, ("prefix_new_id") ) ) {
		print "prefix_new_id parameter missing in the configuration file.\n";
		$error = 1;
	}
	if( ! exists_keys($config, ("check_sequential") ) ) {
		print "check_sequential parameter missing in the configuration file.\n";
		$error = 1;
	}
	if( ! exists_keys($config, ("check_l2_linked_to_l3") ) ) {
		print "check_l2_linked_to_l3 parameter missing in the configuration file.\n";
		$error = 1;
	}
	if( ! exists_keys($config, ("check_l1_linked_to_l2") ) ) {
		print "check_l1_linked_to_l2 parameter missing in the configuration file.\n";
		$error = 1;
	}
	if( ! exists_keys($config, ("remove_orphan_l1") ) ) {
		print "remove_orphan_l1 parameter missing in the configuration file.\n";
		$error = 1;
	}
	if( ! exists_keys($config, ("check_all_level3_locations") ) ) {
		print "check_all_level3_locations parameter missing in the configuration file.\n";
		$error = 1;
	}
	if( ! exists_keys($config, ("check_cds") ) ) {
		print "check_cds parameter missing in the configuration file.\n";
		$error = 1;
	}
	if( ! exists_keys($config, ("check_exons") ) ) {
		print "check_exons parameter missing in the configuration file.\n";
		$error = 1;
	}
	if( ! exists_keys($config, ("check_utrs") ) ) {
		print "check_utrs parameter missing in the configuration file.\n";
		$error = 1;
	}
	if( ! exists_keys($config, ("check_all_level2_locations") ) ) {
		print "check_all_level2_locations parameter missing in the configuration file.\n";
		$error = 1;
	}
	if( ! exists_keys($config, ("check_all_level1_locations") ) ) {
		print "check_all_level1_locations parameter missing in the configuration file.\n";
		$error = 1;
	}
	if( ! exists_keys($config, ("check_identical_isoforms") ) ) {
		print "check_identical_isoforms parameter missing in the configuration file.\n";
		$error = 1;
	}

	# Now exit if one confiuguration parameter is missing
	if ($error){exit 1;}
}

1;
