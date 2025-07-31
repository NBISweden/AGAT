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
our @EXPORT = qw( load_levels expose_levels get_feature_type_by_agat_value );

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

	# Fill hash with the features and their values
	foreach my $tag ( keys %{$LEVELS->{$level}} ){
		if($LEVELS->{$level}{ lc($tag) } eq lc($value) ){
			$list_features->{ lc($tag) } = lc ($value);
		}
	}
	return $list_features;
}

# @Purpose: set path to look at the feature level files (If present locally we take them otherwise look at standard path).
# If expose option is activated, we copy the file localy and exit
# We save the Levels in the LEVEL variable accessible here, and in the hash
# @input: 3 =>	hash, string (path), integer
# @output: 0 => none
# @Remark: none
sub load_levels{

	#check first if exist locally
	my $run_dir = cwd;
	my $path = $run_dir."/".$feature_levels_file;
	my $message = "Using local";

	# If nothing local take the one shipped with AGAT
	if (! -e $path) {
		$path = dist_file('AGAT', $feature_levels_file);
		$message = "Using standard";
		dual_print ({ string => "Path where $feature_levels_file is standing according to dist_file: $path\n", debug_only => 1 });
	}
	dual_print ({ string => "$message $path file\n"});

	# Load the yaml files as hash
	my $feature_levels_hash = LoadFile($path);

	# retrun hash
	return $feature_levels_hash;
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
