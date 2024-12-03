#!/usr/bin/perl -w

package AGAT::Levels;

use strict;
use warnings;
use YAML qw(DumpFile LoadFile);
use File::Copy;
use File::ShareDir ':ALL';
use AGAT::Utilities;
use Cwd qw(cwd);
use Exporter;


our @ISA = qw(Exporter);
our @EXPORT = qw( get_levels_info load_levels expose_levels get_feature_type_by_agat_value );

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

#	-----------------------------------CONSTANT-----------------------------------	 

my $feature_levels_file = ('feature_levels.yaml');

#	------------------------------------GENERAL------------------------------------	 

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

# @Purpose: set path to look at the feature level files (If present locally we take them otherwise look at standard path).
# If expose option is activated, we copy the file localy and exit
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
