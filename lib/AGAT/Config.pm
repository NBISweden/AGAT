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

# This is the default one, but can be changed by the user
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
	my ($path);
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
	my ( $type, $config_file_in ) ;

	if( ! defined($args->{type}) ) { $type = "local";} else{ $type = $args->{type};}
	if( !$type eq "local" and !$type eq "original" ){
		die "type must be local or original. $type unknown!";
	}
	if( ! defined($args->{config_file_in}) ) { $config_file_in = undef;} else{ $config_file_in = $args->{config_file_in}; }
	
	my $path=undef;
	my $log_message = "";
	# original approach trying to get the local and/or original config file
	if (! $config_file_in){
		#set run directory
		my $run_dir = cwd;

		# get local config if any
		if ($type eq "local") {
			$path = $run_dir."/".$config_file;
			if (-e $path){
				$log_message = "=> Using $config_file config file found in your working directory.\n";
			} else {
				$path = undef;
			}
		}
		#otherwise use the standard location ones
		if (! $path) { 
			$path = dist_file('AGAT', $config_file);
			$log_message = "=> Using standard $path config file\n";
		}
	}
	# Config file provided we must load this one !
	else{
		if (-e $config_file_in){
			$path = $config_file_in;
			$log_message = "=> Using provided config file $path.\n";
		} else{
			die "=> Config file provided $config_file_in does not exist! Please check the path!";
		}
	}
	return $path, $log_message;
}

sub expose_config_hash{
	my ($args)=@_;

	my ($config_in, $config_file_out);
	if( ! defined($args->{config_in}) ) { $config_in = undef;} else{ $config_in = $args->{config_in};}
	if( ! defined($args->{config_file_out}) ) { $config_file_out = $config_file;} else{ $config_file_out = $args->{config_file_out};}

	DumpFile($config_file_out, $config_in);
}

# @Purpose: Write the config hash in a yaml file in the current directory 
sub expose_config_file{
	my ($args)=@_;

	my ($path_in, $config_file_out);
	if( ! defined($args->{config_file_in}) ) { $path_in = undef;} else{ $path_in = $args->{config_file_in};}
	if( ! defined($args->{config_file_out}) ) { $config_file_out = $config_file;} else{ $config_file_out = $args->{config_file_out};}

	#set run directory
	if(! $path_in){
		$path_in = dist_file('AGAT', $config_file);
		print "Path where $config_file is standing according to dist_file: $path_in\n";
	}
	# copy the file locally
	my $run_dir = cwd;
	copy($path_in, $run_dir."/".$config_file_out) or die print "Copy failed: $!";
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
			warn "Verbose parameter must be between 0 and 4.\n";
			$error = 1;
		}
	}
	else{
		warn "verbose parameter missing in the configuration file.\n";
		$error = 1;
	}
	if( !  exists_keys($config,("progress_bar") ) ){
		warn "progress_bar parameter missing in the configuration file.\n";
		$error = 1;
	}
	if( !  exists_keys($config,("cpu") ) ){
		warn "cpu parameter missing in the configuration file.\n";
		$error = 1;
	}
	if( !  exists_keys($config,("minimum_chunk_size") ) ){
		warn "minimum_chunk_size parameter missing in the configuration file.\n";
		$error = 1;
	}
	if( !  exists_keys($config,("log") ) ){
		warn "log parameter missing in the configuration file.\n";
		$error = 1;
	}
	if( !  exists_keys($config, ("debug") ) ){
		warn "debug parameter missing in the configuration file.\n";
		$error = 1;
	}
	if( !  exists_keys($config, ("tabix") ) ){
		warn "tabix parameter missing in the configuration file.\n";
		$error = 1;
	}
	if( !  exists_keys($config, ("merge_loci") ) ) {
		warn "merge_loci parameter missing in the configuration file.\n";
		$error = 1;
	}
	if( !  exists_keys($config, ("throw_fasta") ) ) {
		warn "throw_fasta parameter missing in the configuration file.\n";
		$error = 1;
	}
	if( exists_keys($config, ("force_gff_input_version") ) ) {
		my %values = (0 => 1, 1 => 1, 2 => 1, 2.5 => 1, 3 => 1);
		if (! exists_keys(\%values,($config->{ force_gff_input_version }) ) ) {
			warn "force_gff_input_version parameter must be 0, 1, 2, 2.5 or 3.\n";
			$error = 1;
		}
	}
	else{
		warn "force_gff_input_version parameter missing in the configuration file.\n";
		$error = 1;
	}

	if( exists_keys($config, ("output_format") ) ) {
		my %values = ( gff => 1, gtf => 1 ); # case was set to lowercase when loaded into hash
		$config->{ output_format } = lc($config->{ output_format }); # set it to lowercase
		if (! exists_keys(\%values,( $config->{ output_format } ) ) ) {
			warn "output_format parameter must be GFF or GTF.\n";
			$error = 1;
		}
	}
	else{
		warn "output_format parameter missing in the configuration file.\n";
		$error = 1;
	}

	if( exists_keys($config, ("gff_output_version") ) ) {
		my %values = ( 1 => 1, 2 => 1, 2.5 => 1, 3 => 1 );
		if (! exists_keys(\%values,($config->{ gff_output_version }) ) ) {
			warn "gff_output_version parameter must be 1, 2, 2.5 or 3.\n";
			$error = 1;
		}
	}
	else{
		warn "gff_output_version parameter missing in the configuration file.\n";
		$error = 1;
	}
	if( exists_keys($config, ("gtf_output_version") ) ) {
		my %values = ( 1 => 1, 2 => 1, 2.1 => 1, 2.2 => 1, 2.5 => 1, 3 => 1, relax => 1 );
		if (! exists_keys(\%values,($config->{ gtf_output_version }) ) ) {
			warn "gtf_output_version parameter must be 1, 2, 2.1, 2.2, 2.5, 3 or relax.\n";
			$error = 1;
		}
	}
	else{
		warn "gtf_output_version parameter missing in the configuration file.\n";
		$error = 1;
	}
	if( ! exists_keys($config, ("deflate_attribute") ) ) {
		print "deflate_attribute parameter missing in the configuration file.\n";
		$error = 1;
	}
	if( ! exists_keys($config, ("create_l3_for_l2_orphan") ) ) {
		warn "create_level3_for_level2_orphan parameter missing in the configuration file.\n";
		$error = 1;
	}
	if( ! exists_keys($config, ("locus_tag") ) ) {
		warn "locus_tag parameter missing in the configuration file.\n";
		$error = 1;
	}
	if( !  exists_keys($config, ("prefix_new_id") ) ) {
		warn "prefix_new_id parameter missing in the configuration file.\n";
		$error = 1;
	}
	if( ! exists_keys($config, ("check_sequential") ) ) {
		warn "check_sequential parameter missing in the configuration file.\n";
		$error = 1;
	}
	if( ! exists_keys($config, ("check_l2_linked_to_l3") ) ) {
		warn "check_l2_linked_to_l3 parameter missing in the configuration file.\n";
		$error = 1;
	}
	if( ! exists_keys($config, ("check_l1_linked_to_l2") ) ) {
		warn "check_l1_linked_to_l2 parameter missing in the configuration file.\n";
		$error = 1;
	}
	if( ! exists_keys($config, ("remove_orphan_l1") ) ) {
		warn "remove_orphan_l1 parameter missing in the configuration file.\n";
		$error = 1;
	}
	if( ! exists_keys($config, ("check_all_level3_locations") ) ) {
		warn "check_all_level3_locations parameter missing in the configuration file.\n";
		$error = 1;
	}
	if( ! exists_keys($config, ("check_cds") ) ) {
		warn "check_cds parameter missing in the configuration file.\n";
		$error = 1;
	}
	if( ! exists_keys($config, ("check_exons") ) ) {
		warn "check_exons parameter missing in the configuration file.\n";
		$error = 1;
	}
	if( ! exists_keys($config, ("check_utrs") ) ) {
		warn "check_utrs parameter missing in the configuration file.\n";
		$error = 1;
	}
	if( ! exists_keys($config, ("check_all_level2_locations") ) ) {
		warn "check_all_level2_locations parameter missing in the configuration file.\n";
		$error = 1;
	}
	if( ! exists_keys($config, ("check_all_level1_locations") ) ) {
		warn "check_all_level1_locations parameter missing in the configuration file.\n";
		$error = 1;
	}
	if( ! exists_keys($config, ("check_identical_isoforms") ) ) {
		warn "check_identical_isoforms parameter missing in the configuration file.\n";
		$error = 1;
	}
	if( ! exists_keys($config, ("clean_attributes_from_template") ) ) {
		warn "clean_attributes_from_template parameter missing in the configuration file.\n";
		$error = 1;
	}

	# Now exit if one configuration parameter is missing
	if ($error){exit 1;}
}

1;
