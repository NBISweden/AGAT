#!/usr/bin/perl -w

package AGAT::OmniscientJson;

use strict;
use warnings;
use Bio::Tools::GFF;
use Bio::Seq;
use JSON;
use Try::Tiny;
use File::ShareDir ':ALL';
use AGAT::Utilities;
use Cwd qw(cwd);
use Exporter;



our @ISA = qw(Exporter);
our @EXPORT = qw(get_feature_type_by_agat_value get_levels_info load_levels load_json);
sub import {
  AGAT::OmniscientJson->export_to_level(1, @_); # to be able to load the EXPORT functions when direct call; (normal case)
  AGAT::OmniscientJson->export_to_level(2, @_); # to be able to load the EXPORT functions when called from one level up;
}

=head1 SYNOPSIS

This is the code to handle AGAT's json value. Accessing file is slow, the
values from the json files must be stored within OMNISCIENT when accessed the
first time.

=head1 DESCRIPTION

 This lib contains functions to parse JSON files and store the values within
 OMNISCIENT and functions to access value from OMNISCIENT.

=head1 AUTHOR

  Jacques Dainat - jacques.dainat@nbis.se

=cut

#				   +------------------------------------------------------+
#				   |+----------------------------------------------------+|
#				   || 			HANDLE OTHER part of OMNISCIENT		 ||
#				   |+----------------------------------------------------+|
#				   +------------------------------------------------------+

# @Purpose: get value provided by the json file related to the feature
# @input: 2 => hash(omniscient hash), level, value
# @output: 1 hash => {tag}=value
sub get_feature_type_by_agat_value {
	my ($omniscient, $level, $value) = @_;

	my $list_features = {};
	my $hash = undef;
	# If info not in omniscient we append omniscient to include all info
	if (! exists_keys ($omniscient, ('other', 'level', 'level1') ) ){
		$hash = get_levels_info() if (! $hash); # get from the file
		$omniscient->{'other'}{'level'}{'level1'} = $hash->{'other'}{'level'}{'level1'};
	}
	elsif(! exists_keys ($omniscient, ('other', 'level', 'level2') ) ){
		$hash = get_levels_info() if (! $hash); # get from the file
		$omniscient->{'other'}{'level'}{'level2'} = $hash->{'other'}{'level'}{'level2'};
	}
	elsif(! exists_keys ($omniscient, ('other', 'level', 'level3') ) ){
		$hash = get_levels_info() if (! $hash); # get from the file
		$omniscient->{'other'}{'level'}{'level3'} = $hash->{'other'}{'level'}{'level3'};
	}
	elsif(! exists_keys ($omniscient, ('other', 'level', 'spreadfeature') ) ){
		$hash = get_levels_info() if (! $hash); # get from the file
		$omniscient->{'other'}{'level'}{'spreadfeature'} = $hash->{'other'}{'level'}{'spreadfeature'};
	}

	foreach my $tag ( keys %{$omniscient->{'other'}{'level'}{$level}} ){
		if($omniscient->{'other'}{'level'}{$level}{ lc($tag) } eq lc($value) ){
			$list_features->{ lc($tag) } = lc ($value);
		}
	}
	return $list_features;
}

# @Purpose: We load levels from agat's json and save it in a hash similar as we save it in omniscient
# @input: 2 =>	hash, integer
# @output: 3 => hash
# @Remark: none
sub get_levels_info{
		my ($hash, $verbose) = @_ ;

		$hash = {} if (! $hash); # if the hash exist we will append it otherwise it will be a new one
		$verbose = 0 if(! defined ($verbose));
		load_levels($hash,undef,$verbose);
		return $hash;
}

# @Purpose: set path to look at the json feature level files (If present locally we take them otherwise look at standard path).
# If expose option is activated, we copy the json files localy and exit
# We save the Levels in the LEVEL variable accessible here, and in the hash
# @input: 3 =>	hash, string (path), integer
# @output: 0 => none
# @Remark: none
sub load_levels{
	my ($hash_omniscient, $expose_feature_levels, $verbose) = @_ ;

	$verbose = 0 if(! $verbose );

	print "	 Accessing the feature level files:\n" if($verbose > 0);
	#set original path to json files, order matter
	my @files = ('features_level1.json', 'features_level2.json', 'features_level3.json', 'features_spread.json');
	my @paths;
	foreach my $file ( @files ){
		my $path = dist_file('AGAT', $file);
		print "Path where $file is standing according to dist_file: $path\n" if ($verbose > 2);
		push @paths, $path;
	}

	#set run directory
	my $run_dir = cwd;
	# Check if it is asked to copy the json files locally
	if ($expose_feature_levels){
		foreach my $path (@paths) {
				copy($path, $run_dir) or die "Copy failed: $!";
		}
		print "			All json feature level files copied in your working directory\n" if ($verbose);
		exit;
	}
	# Load the json files
	else{
		my $cpt=1;
		foreach my $file (@files) {
			#check first if exist locally
			my $path = $run_dir."/".$file;
			if (-e $path) {

				print "			Using local $file file\n" if($verbose > 0);

				if ($cpt == 1){
					$hash_omniscient->{'other'}{'level'}{'level1'} = load_json($path);
				}
				elsif ($cpt == 2){
					$hash_omniscient->{'other'}{'level'}{'level2'} = load_json($path);
				}
				elsif ($cpt == 3){
					$hash_omniscient->{'other'}{'level'}{'level3'} = load_json($path);;
				}
				else {
					$hash_omniscient->{'other'}{'level'}{'spreadfeature'} = load_json($path);
				}
			}
			else{ #otherwise use the standard location ones

				print "			Using standard ".$paths[$cpt-1]." file\n" if($verbose > 0);

				if ($cpt == 1){
					$hash_omniscient->{'other'}{'level'}{'level1'} = load_json($paths[0]);
				}
				elsif ($cpt == 2){
					$hash_omniscient->{'other'}{'level'}{'level2'} =  load_json($paths[1]);
				}
				elsif ($cpt == 3){
					$hash_omniscient->{'other'}{'level'}{'level3'} = load_json($paths[2]);
				}
				else {
					$hash_omniscient->{'other'}{'level'}{'spreadfeature'} = load_json($paths[3]);
				}
			}
			$cpt++;
		}
	}
}

# @Purpose: load json data into variable
# @input: 3 =>	String path to the json file
# @output: 1 => hash reference with data
# @Remark: none
sub load_json{

	my ($file_path) = @_;

	my $result = undef;
	my $json_text = do {
		open(my $json_fh, "<:encoding(UTF-8)", $file_path)
		or die("load_json: Can't open $file_path: $!\n");
		local $/;
		<$json_fh>
	};

	my $json = JSON->new;
	try{
		$result = $json->decode($json_text);
	}
	catch{
		print "error while parsing $file_path. Please verify the sanity of your json file.\n";
	};

	return $result;
}


#				   +-------------------------  END -----------------------------+
1;
