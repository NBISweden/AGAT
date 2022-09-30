#!/usr/bin/perl -w

package AGAT::OmniscientYaml;

use strict;
use warnings;
use Bio::Tools::GFF;
use Bio::Seq;
use YAML qw(DumpFile LoadFile);
use File::Copy;
use Try::Tiny;
use File::ShareDir ':ALL';
use AGAT::Utilities;
use Cwd qw(cwd);
use Exporter;


our @ISA = qw(Exporter);
our @EXPORT = qw(load_config expose_config_file check_config get_config expose_config_hash
	               get_feature_type_by_agat_value get_levels_info load_levels load_json expose_levels);
sub import {
  AGAT::OmniscientYaml->export_to_level(1, @_); # to be able to load the EXPORT functions when direct call; (normal case)
  AGAT::OmniscientYaml->export_to_level(2, @_); # to be able to load the EXPORT functions when called from one level up;
}

=head1 SYNOPSIS

This is the code to handle AGAT's yaml value. Accessing file is slow, the
values from the yaml files must be stored within OMNISCIENT when accessed the
first time.

=head1 DESCRIPTION

 This lib contains functions to parse YAML files and store the values within
 OMNISCIENT and functions to access value from OMNISCIENT.

=head1 AUTHOR

  Jacques Dainat - jacques.dainat@nbis.se

=cut

#	------------------------------------GENERAL------------------------------------	 

#	------------------------------------CONFIG------------------------------------	 


my $config_file= ('config.yaml');

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

	DumpFile('config.yaml', $config);
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
	# copy the json files locally
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
			print "force_gff_input_version parameter must be 0,1,2,2.5 or 3.\n";
			$error = 1;
		}
	}
	else{
		print "force_gff_input_version parameter missing in the configuration file.\n";
		$error = 1;
	}
	if( exists_keys($config, ("gff_output_version") ) ) {
		my %values = ( 1 => 1, 2 => 1, 2.5 => 1, 3 => 1);
		if (! exists_keys(\%values,($config->{ gff_output_version }) ) ) {
			print "gff_output_version parameter must be 1,2,2.5 or 3.\n";
			$error = 1;
		}
	}
	else{
		print "gff_output_version parameter missing in the configuration file.\n";
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
# +----------------------------- CONFIG  END ----------------------------------+

#	--------------------------------FEATURE LEVELS--------------------------------	 

my $feature_levels_file = ('feature_levels.yaml');

# @Purpose: get value provided by the yaml file related to the feature
# @input: 2 => hash(omniscient hash), level, value
# @output: 1 hash => {tag}=value
sub get_feature_type_by_agat_value {
	my ($omniscient, $level, $value) = @_;

	my $list_features = {};
	my $hash = undef;
	# If info not in omniscient we append omniscient to include all info
	if (! exists_keys ($omniscient, ('other', 'level') ) ){
		$hash = get_levels_info({verbose => 0}) if (! $hash); # get from the file
		$omniscient->{'other'}{'level'} = $hash;
	}

	# Fill hash with the features and their values
	foreach my $tag ( keys %{$omniscient->{'other'}{'level'}{$level}} ){
		if($omniscient->{'other'}{'level'}{$level}{ lc($tag) } eq lc($value) ){
			$list_features->{ lc($tag) } = lc ($value);
		}
	}
	return $list_features;
}

# @Purpose: We load levels from agat's yaml file and save it in a hash similar as we save it in omniscient
# @input: 2 =>	hash, integer
# @output: 3 => hash
# @Remark: none
sub get_levels_info{
	my ($args) = @_;

	my ($hash, $verbose);
	# if the hash exist we will append it otherwise it will be a new one
	if( ! defined($args->{omniscient})) { $hash = {} } else{ $hash = $args->{omniscient}; }
	#size line
	if( ! defined($args->{verbose}) ) { $verbose = 0;} else{ $verbose = $args->{verbose}; }

	load_levels({ omniscient => $hash, verbose => $verbose});

	return $hash;
}

# @Purpose: set path to look at the json feature level files (If present locally we take them otherwise look at standard path).
# If expose option is activated, we copy the json files localy and exit
# We save the Levels in the LEVEL variable accessible here, and in the hash
# @input: 3 =>	hash, string (path), integer
# @output: 0 => none
# @Remark: none
sub load_levels{
	my ($args) = @_;

	# -------------- INPUT --------------
	# Check we receive a hash as ref
	if(ref($args) ne 'HASH'){ warn "Hash Arguments expected for load_levels. Please check the call.\n";exit;}
	# -- Declare all variables and fill them --
	my ($hash_omniscient, $verbose, $log, $debug);
	# string to print
	if( defined($args->{omniscient})) {$hash_omniscient = $args->{omniscient};} else{ warn "Omniscient input is mandatory for load_levels\n"; exit;}
	# character to fill the line with
	if( ! defined($args->{verbose}) ) { $verbose = undef;} else{ $verbose = $args->{verbose}; }
	# log
	if( ! defined($args->{log}) ) { $log = undef;} else{ $log = $args->{log}; }
	# log
	if( ! defined($args->{debug}) ) { $debug = undef;} else{ $debug = $args->{debug}; }

	#check first if exist locally
	my $run_dir = cwd;
	my $path = $run_dir."/".$feature_levels_file;
	my $message = "Using local";

	# If nothing local take the one shipped with AGAT
	if (! -e $path) {
		$path = dist_file('AGAT', $feature_levels_file);
		$message = "Using standard";
		dual_print ($log, "Path where $feature_levels_file is standing according to dist_file: $path\n", $verbose) if ($debug);
	
	}
	dual_print($log, "$message $path file\n", $verbose );

	# Load the yaml files as hash
	my $feature_levels_hash = LoadFile($path);

	# Save the data within omniscient
	$hash_omniscient->{'other'}{'level'} = $feature_levels_hash;
}


# @Purpose: copy the yaml level feature file in the current directory 
sub expose_levels{
	
	my	$path = dist_file('AGAT', $feature_levels_file);
	print "Path where $feature_levels_file is standing according to dist_file: $path\n";
	
	# copy the json files locally
	my $run_dir = cwd;
	copy($path, $run_dir) or die print "Copy failed: $!";
}

1;
