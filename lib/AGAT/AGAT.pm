#!/usr/bin/perl -w

package AGAT::AGAT;

use strict;
use warnings;
use Exporter;
use Pod::Usage;
use POSIX qw(strftime);
use File::Basename;
use File::Path qw(remove_tree);
use AGAT::OmniscientI;
use AGAT::OmniscientO;
use AGAT::OmniscientTool;
use AGAT::Config;
use AGAT::Levels;
use AGAT::OmniscientStat;
use AGAT::Utilities;
use AGAT::PlotR;
use Bio::Tools::GFF;

our $VERSION     = "v1.7.0";
our @ISA         = qw(Exporter);
our @EXPORT      = qw(get_agat_header print_agat_version initialize_agat handle_levels create_log_file parse_shared_options apply_shared_options split_argv_shared_vs_script);
sub import {
    AGAT::AGAT->export_to_level(1, @_); # to be able to load the EXPORT functions when direct call; (normal case)
    AGAT::OmniscientI->export_to_level(1, @_);
    AGAT::OmniscientO->export_to_level(1, @_);
    AGAT::OmniscientTool->export_to_level(1, @_);
    AGAT::Config->export_to_level(1, @_);
    AGAT::Levels->export_to_level(1, @_);
    AGAT::OmniscientStat->export_to_level(1, @_);
    AGAT::Utilities->export_to_level(1, @_);
    AGAT::PlotR->export_to_level(1, @_);
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

# initialize_agat:
# set logging and save it in global $LOGGING
# load configuration file from local file if any either the one shipped with AGAT and save it in global $CONFIG 
# load level file from local file if any either the one shipped with AGAT and save it in global $LEVELS

sub initialize_agat{
	my ($args)=@_;

	# Needed if log activated
	my ($input, $config_file_provided, $shared_opts);
	if ( defined($args->{input}) ) { $input = $args->{input}; } 
	if ( defined($args->{config_file_in}) ) { $config_file_provided = $args->{config_file_in};}
	if ( defined($args->{shared_opts}) ) { $shared_opts = $args->{shared_opts}; }

	# Get the config file
	my ($config_file_checked, $log_message) = get_config({type => "local", config_file_in => $config_file_provided}); #try local first, if none will take the original
	# Load the config and check it
	$CONFIG = load_config({ config_file => $config_file_checked});	
	# Apply shared options overrides into global $CONFIG if provided
	apply_shared_options($shared_opts) if (defined $shared_opts);
	
	# --- logging --- LOGGING DEFINED INTO UTILITIES
	$LOGGING = {'verbose' => $CONFIG->{'verbose'}, 'debug_mode' => $CONFIG->{'debug'} };

	# Print the header. Put here because get_agat_config it the first function call for _sp_ and _sq_ screen
	print AGAT::AGAT::get_agat_header() if ( $CONFIG->{verbose} );
	
	# Create log file if needed
	if ( $CONFIG->{log} ){
		if (! $input){
			my ($package, $filename, $line, $subroutine) = caller(0);
			die "No input file provided for log naming set up! Called from subroutine: $subroutine at $filename line $line\n";
		}
		my $log = create_log_file({input => $input});
		$LOGGING->{'log'} = $log ;
		# +----------------- Print header ------------------+
		dual_print ({ string => AGAT::AGAT::get_agat_header(), log_only => 1 });
		dual_print ({ string => $log_message, log_only => 1 });
	}

	# --- set LEVELS variable ---
	$LEVELS = load_levels();
}

# Parse a common set of AGAT shared options from an argv arrayref
# Returns (\%shared_opts, \@residual_args)
sub parse_shared_options {
	my ($argv_ref) = @_;
	my @argv = @{$argv_ref // []};

	require Getopt::Long;
	# Use a local parser instance so configuration is scoped here only
	my $parser = Getopt::Long::Parser->new();
	# Avoid auto-abbrev so that --gff does not match --gff_output_version
	# No pass_through here: caller partitions argv for shared options only
	$parser->configure(qw(bundling no_auto_abbrev));

	my %opts;
	if ( !$parser->getoptionsfromarray(
		\@argv,
		'check_all_level1_locations!' => \$opts{check_all_level1_locations},
		'check_all_level2_locations!' => \$opts{check_all_level2_locations},
		'check_all_level3_locations!' => \$opts{check_all_level3_locations},
		'check_cds!'                  => \$opts{check_cds},
		'check_exons!'                => \$opts{check_exons},
		'check_identical_isoforms!'   => \$opts{check_identical_isoforms},
		'check_l1_linked_to_l2!'      => \$opts{check_l1_linked_to_l2},
		'check_l2_linked_to_l3!'      => \$opts{check_l2_linked_to_l3},
		'check_sequential!'           => \$opts{check_sequential},
		'check_utrs!'                 => \$opts{check_utrs},
		'clean_attributes_from_template!' => \$opts{clean_attributes_from_template},
		'config=s'                     => \$opts{config},
		'cpu|thread|core|job=i'        => \$opts{cpu},
		'create_l3_for_l2_orphan!'    => \$opts{create_l3_for_l2_orphan},
		'debug!'                      => \$opts{debug},
		'deflate_attribute!'          => \$opts{deflate_attribute},
		'force!'                      => \$opts{force},
		'force_gff_input_version=i'    => \$opts{force_gff_input_version},
		'gff_output_version=i'         => \$opts{gff_output_version},
		'gtf_output_version=s'         => \$opts{gtf_output_version},
		'locus_tag=s'                  => \$opts{locus_tag},
		'log!'                        => \$opts{log},
		'merge_loci!'                 => \$opts{merge_loci},
		'minimum_chunk_size=i'         => \$opts{minimum_chunk_size},
		'output_format=s'              => \$opts{output_format},
		'prefix_new_id=s'              => \$opts{prefix_new_id},
		'progress_bar!'               => \$opts{progress_bar},
		'remove_orphan_l1!'           => \$opts{remove_orphan_l1},
		'tabix!'                      => \$opts{tabix},
		'throw_fasta!'                => \$opts{throw_fasta},
		'url_encode_out!'             => \$opts{url_encode_out},
		'v|verbose=i'                  => \$opts{verbose},
	) ) {
    pod2usage( { -message => 'Failed to parse command line',
                 -verbose => 1,
                 -exitval => 1 } );
}
	return (\%opts);
}

# Partition argv into shared vs script-specific option lists
# Returns (\@shared_args, \@script_args)
sub split_argv_shared_vs_script {
	my ($argv_ref) = @_;
	my @argv = @{$argv_ref // []};

	my %shared = map { $_ => 1 } qw(
		config cpu thread core job minimum_chunk_size
		v verbose progress_bar log debug tabix merge_loci throw_fasta force
		force_gff_input_version output_format gff_output_version gtf_output_version
		deflate_attribute create_l3_for_l2_orphan clean_attributes_from_template
		locus_tag check_sequential check_l2_linked_to_l3 check_l1_linked_to_l2
		remove_orphan_l1 check_all_level3_locations check_cds check_exons check_utrs
		check_all_level2_locations check_all_level1_locations check_identical_isoforms
		prefix_new_id url_encode_out
	);
	
	# Boolean options that don't take arguments
	my %boolean_opts = map { $_ => 1 } qw(
		progress_bar log debug tabix merge_loci throw_fasta force
		deflate_attribute create_l3_for_l2_orphan clean_attributes_from_template
		check_sequential check_l2_linked_to_l3 check_l1_linked_to_l2
		remove_orphan_l1 check_all_level3_locations check_cds check_exons check_utrs
		check_all_level2_locations check_all_level1_locations check_identical_isoforms
		url_encode_out
	);
	
	my %script = map { $_ => 1 } qw(g gxf gtf gff o output h help);

	my (@shared_args, @script_args);
	for (my $i = 0; $i < @argv; $i++) {
		my $t = $argv[$i];
		if ($t =~ /^--?([^=]+)(?:=.*)?$/) {
			my $name = $1;
			
			# Remove no- prefix to get base name for lookup
			my $base_name = $name;
			$base_name =~ s/^no-//;

			if ($shared{$name} || $shared{$base_name}) {
				push @shared_args, $t;
				# Only capture next arg if NOT a boolean option and NOT using = syntax
				if ($t !~ /=/ && !$boolean_opts{$base_name} && ($i+1) < @argv && $argv[$i+1] !~ /^-/) {
					push @shared_args, $argv[++$i];
				}
			} elsif ($script{$name} || $script{$base_name}) {
				push @script_args, $t;
				if ($t !~ /=/ && ($i+1) < @argv && $argv[$i+1] !~ /^-/) {
					push @script_args, $argv[++$i];
				}
			} else {
				# Unknown option: default to script args to preserve behavior
				push @script_args, $t;
				if ($t !~ /=/ && ($i+1) < @argv && $argv[$i+1] !~ /^-/) {
					push @script_args, $argv[++$i];
				}
			}
		} else {
			# Non-option token: send to script args (future positional)
			push @script_args, $t;
		}
	}
	return (\@shared_args, \@script_args);
}

# Apply shared options into global $CONFIG (expects initialize_agat already run)
sub apply_shared_options {
	my ($opts_hr) = @_;
	my %opts = %{ $opts_hr // {} };

	# Numeric
	$CONFIG->{ cpu } = $opts{cpu} if defined $opts{cpu};
	$CONFIG->{ minimum_chunk_size } = $opts{minimum_chunk_size} if defined $opts{minimum_chunk_size};
	$CONFIG->{ verbose } = $opts{verbose} if defined $opts{verbose};
	$CONFIG->{ force_gff_input_version } = $opts{force_gff_input_version} if defined $opts{force_gff_input_version};
	$CONFIG->{ gff_output_version } = $opts{gff_output_version} if defined $opts{gff_output_version};

	# Strings / enums
	$CONFIG->{ output_format } = lc($opts{output_format}) if defined $opts{output_format};
	$CONFIG->{ gtf_output_version } = lc($opts{gtf_output_version}) if defined $opts{gtf_output_version};
	$CONFIG->{ prefix_new_id } = $opts{prefix_new_id} if defined $opts{prefix_new_id};

	# Booleans (stored as 'true'/'false')
	for my $k (qw(check_all_level1_locations check_all_level2_locations check_all_level3_locations check_cds check_exons check_identical_isoforms check_l1_linked_to_l2 check_l2_linked_to_l3 check_sequential check_utrs clean_attributes_from_template create_l3_for_l2_orphan debug deflate_attribute force log merge_loci progress_bar remove_orphan_l1 tabix throw_fasta url_encode_out)) {
		if (defined $opts{$k}) {
			# With ! specification: 0 or 1
			$CONFIG->{$k} = $opts{$k} ? 1 : 0;
		}
	}

	# List
	if (defined $opts{locus_tag}) {
		my @list = split(/,/, $opts{locus_tag});
		$CONFIG->{ locus_tag } = \@list;
	}

	# Re-check config validity
	check_config({ config => $CONFIG });

	# false does not exists in perl, must be replaced by 0
	normalize_false_to_zero($CONFIG);
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
		my $h = $general->{configs}[-1]{h};

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
		if(! $help and ! $h and $nb_args != 0){ exit 0;}
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

		# config file option
		my $check_all_level1_locations = $general->{configs}[-1]{check_all_level1_locations};
		my $check_all_level2_locations = $general->{configs}[-1]{check_all_level2_locations};
		my $check_all_level3_locations = $general->{configs}[-1]{check_all_level3_locations};
		my $check_cds = $general->{configs}[-1]{check_cds};
		my $check_exons = $general->{configs}[-1]{check_exons};
		my $check_identical_isoforms = $general->{configs}[-1]{check_identical_isoforms};
		my $check_l1_linked_to_l2 = $general->{configs}[-1]{check_l1_linked_to_l2};
		my $check_l2_linked_to_l3 = $general->{configs}[-1]{check_l2_linked_to_l3};
		my $check_sequential = $general->{configs}[-1]{check_sequential};
		my $check_utrs = $general->{configs}[-1]{check_utrs};
		my $clean_attributes_from_template = $general->{configs}[-1]{clean_attributes_from_template};
		my $config_new_name = $general->{configs}[-1]{output};
		my $cpu = $general->{configs}[-1]{cpu};
		my $create_l3_for_l2_orphan = $general->{configs}[-1]{create_l3_for_l2_orphan};
		my $debug = $general->{configs}[-1]{debug};
		my $deflate_attribute = $general->{configs}[-1]{deflate_attribute};
		my $expose = $general->{configs}[-1]{expose};
		my $force = $general->{configs}[-1]{force};
		my $force_gff_input_version = $general->{configs}[-1]{force_gff_input_version};
		my $gff_output_version = $general->{configs}[-1]{gff_output_version};
		my $gtf_output_version = $general->{configs}[-1]{gtf_output_version};
		my $help = $general->{configs}[-1]{help};
		my $locus_tag = $general->{configs}[-1]{locus_tag};
		my $log = $general->{configs}[-1]{log};
		my $merge_loci = $general->{configs}[-1]{merge_loci};
		my $minimum_chunk_size = $general->{configs}[-1]{minimum_chunk_size};
		my $output_format = $general->{configs}[-1]{output_format};
		my $prefix_new_id = $general->{configs}[-1]{prefix_new_id};
		my $progress_bar = $general->{configs}[-1]{progress_bar};
		my $remove_orphan_l1 = $general->{configs}[-1]{remove_orphan_l1};
		my $tabix = $general->{configs}[-1]{tabix};
		my $throw_fasta = $general->{configs}[-1]{throw_fasta};
		my $url_encode_out = $general->{configs}[-1]{url_encode_out};
		my $verbose = $general->{configs}[-1]{verbose};

		# Deal with Expose feature OPTION
		if($expose){
			my ($config_file, $log_info) = get_config({type => "original"});
			my $config = load_config({ config_file => $config_file});

			# set config params on the fly
			my $modified_on_the_fly = undef;

			# bolean
			if( defined($check_all_level1_locations) ){
				$config->{ check_all_level1_locations } = _make_bolean($check_all_level1_locations);
				$modified_on_the_fly = 1;
			}
			# bolean
			if( defined($check_all_level2_locations) ){
				$config->{ check_all_level2_locations } = _make_bolean($check_all_level2_locations);
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
			if( defined($check_identical_isoforms) ){
				$config->{ check_identical_isoforms } = _make_bolean($check_identical_isoforms);
				$modified_on_the_fly = 1;
			}
			# bolean
			if( defined($check_l1_linked_to_l2) ){
				$config->{ check_l1_linked_to_l2 } = _make_bolean($check_l1_linked_to_l2);
				$modified_on_the_fly = 1;
			}
			# bolean
			if( defined($check_l2_linked_to_l3) ){
				$config->{ check_l2_linked_to_l3 } = _make_bolean($check_l2_linked_to_l3);
				$modified_on_the_fly = 1;
			}
			# bolean
			if( defined($check_sequential) ){
				$config->{ check_sequential } = _make_bolean($check_sequential);
				$modified_on_the_fly = 1;
			}
			# bolean
			if( defined($check_utrs) ){
				$config->{ check_utrs } = _make_bolean($check_utrs);
				$modified_on_the_fly = 1;
			}
			# bolean
			if( defined($clean_attributes_from_template) ){
				$config->{ clean_attributes_from_template } = _make_bolean($clean_attributes_from_template);
				$modified_on_the_fly = 1;
			}
			# Integer
			if( defined($cpu) ){
				$config->{ cpu } = $cpu;
				$modified_on_the_fly = 1;
			}
			# bolean
			if( defined($create_l3_for_l2_orphan) ){
				$config->{ create_l3_for_l2_orphan } = _make_bolean($create_l3_for_l2_orphan);
				$modified_on_the_fly = 1;
			}
			# bolean
			if( defined($debug) ){
				$config->{ debug } = _make_bolean($debug);
				$modified_on_the_fly = 1;
			}
			# bolean
			if( defined($deflate_attribute) ){
				$config->{ deflate_attribute } = _make_bolean($deflate_attribute);
				$modified_on_the_fly = 1;
			}
			# bolean
			if( defined($force) ){
				$config->{ force } = _make_bolean($force);
				$modified_on_the_fly = 1;
			}
			# integer
			if( defined($force_gff_input_version) ){
				$config->{ force_gff_input_version } = $force_gff_input_version;
				$modified_on_the_fly = 1;
			}
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
			# string
			if( defined($locus_tag) ){
				my @list = split(/,/, $locus_tag);
				$config->{ locus_tag } = \@list;
				$modified_on_the_fly = 1;
			}
			# bolean
			if( defined($log) ){
				$config->{ log } = _make_bolean($log);
				$modified_on_the_fly = 1;
			}
			# bolean
			if( defined($merge_loci) ){
				$config->{ merge_loci } = _make_bolean($merge_loci);
				$modified_on_the_fly = 1;
			}
			# Integer
			if( defined($minimum_chunk_size) ){
				$config->{ minimum_chunk_size } = $minimum_chunk_size;
				$modified_on_the_fly = 1;
			}
			# string
			if( defined($output_format) ){
				$config->{ output_format } = lc($output_format);
				$modified_on_the_fly = 1;
			}
			# string
			if( defined($prefix_new_id) ){
				$config->{ prefix_new_id } = $prefix_new_id;
				$modified_on_the_fly = 1;
			}
			# bolean
			if( defined($progress_bar) ){
				$config->{ progress_bar } = _make_bolean($progress_bar);
				$modified_on_the_fly = 1;
			}
			# bolean
			if( defined($remove_orphan_l1) ){
				$config->{ remove_orphan_l1 } = _make_bolean($remove_orphan_l1);
				$modified_on_the_fly = 1;
			}
			# bolean
			if( defined($tabix) ){
				$config->{ tabix } = _make_bolean($tabix);
				$modified_on_the_fly = 1;
			}
			# bolean
			if( defined($throw_fasta) ){
				$config->{ throw_fasta } = _make_bolean($throw_fasta);
				$modified_on_the_fly = 1;
			}
			# bolean
			if( defined($url_encode_out) ){
				$config->{ url_encode_out } = _make_bolean($url_encode_out);
				$modified_on_the_fly = 1;
			}
			# integer 0-4
			if( defined($verbose) ){
				$config->{ verbose } = $verbose;
				$modified_on_the_fly = 1;
			}

			# Now that CLI overrides are applied, honor verbosity for prints
			if ($config->{verbose} && $config->{verbose} > 0) { print "Config loaded\n"; }

			if ($modified_on_the_fly) {
				if ($config->{verbose} && $config->{verbose} > 0) { print "Config modified\n"; }
			}

			# check config
			check_config({ config => $config});
			if ($config->{verbose} && $config->{verbose} > 0) { print "Config checked\n"; }
			
			if ($modified_on_the_fly) {
				expose_config_hash({ config_in => $config, config_file_out => $config_new_name})
			} else {
				expose_config_file({config_file_in => $config_file, config_file_out => $config_new_name, verbose => $config->{verbose}});
			}

			# inform user
			my $config_file_used;
			if($config_new_name){
				$config_file_used = $config_new_name;
			} else { $config_file_used = "agat_config.yaml"; }
			if ($config->{verbose} && $config->{verbose} > 0) { print "Config file written in your working directory ($config_file_used)\n"; }
		}

		# if help was called (or not arg provided) we let AppEaser continue to print help
		my $nb_args = keys %{$general->{configs}[-1]};
		if(! $help and $nb_args != 0){ exit 0;}

}

# transform 0, false, no into false and everything else into true for displaying into config file when exposed!
sub _make_bolean{
	my ($value) = @_;

	# Treat empty, undef, 0, "0", "false", "no" as false
	if (!defined($value) || $value eq "" || $value eq "0" || lc($value) eq "false" || lc($value) eq "no"){
		return "false";
	}
	return "true";
}

# +----------------- create a log file  ------------------+
sub create_log_file{
	my ($args)=@_;

	my ($input, $log);
	if( defined($args->{input}) ) { $input = $args->{input};}

	if( -f $input){
			my ($filename,$path,$ext) = fileparse($input,qr/\.[^.]*/);
			$AGAT_LOG = $AGAT_LOG."_".$filename;
			
			# create folder if not exist
			if (-d $AGAT_LOG) {
				remove_tree($AGAT_LOG) or die "Failed to delete $AGAT_LOG: $!";
			}
			# create a tmp directory
			mkdir $AGAT_LOG or die "Cannot create directory '$AGAT_LOG': $!";

			# create a log file
			open($log, '>', "$AGAT_LOG/main.log"  ) or
						warn "Can not open $AGAT_LOG/main.log for printing log: $!" && die;
			print $log file_text_line({ string => (strftime "%m/%d/%Y at %Hh%Mm%Ss", localtime),
																  char => " ",
																  extra => "\n"
																  });
		}
		else{
			die "File $input provided as input does not exits! Please verify your path and file existence!";
		}
	return $log;
}

1;
