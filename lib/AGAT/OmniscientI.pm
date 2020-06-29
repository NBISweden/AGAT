#!/usr/bin/perl -w

package AGAT::OmniscientI;

use strict;
use warnings;
use Try::Tiny;
use Bio::Tools::GFF;
use File::Basename;
use File::ShareDir ':ALL';
use POSIX qw(strftime);
use Sort::Naturally;
use LWP::UserAgent;
use Bio::OntologyIO::obo;
use Bio::Ontology::OntologyEngineI;
use Clone 'clone';
use Exporter;
use AGAT::Omniscient;
use AGAT::OmniscientTool;
use AGAT::OmniscientJson;
use AGAT::Utilities;

our @ISA = qw(Exporter);
our @EXPORT = qw(get_level select_gff_format
							modelate_utr_and_cds_features_from_exon_features_and_cds_start_stop
							slurp_gff3_file_JD);
sub import {
	AGAT::OmniscientI->export_to_level(1, @_); # to be able to load the EXPORT functions when direct call; (normal case)
	AGAT::OmniscientI->export_to_level(2, @_); # to be able to load the EXPORT functions when called from one level up;
}

=head1 SYNOPSIS

		This is the code to fill Omniscient data structure by parsing any kind of gtf/gff files.
		(It could also fill Omniscient data structure by parsing hash and list if data is provided in a proper way)

=head1 DESCRIPTION

		A library to convert handle any kind of gff file and save it in memory as GFF3 "full" and standard format.
		Full format means, we expand exon having several parents, we add ID everywhere (even if level3 ID is not mandatory), and Parent everywhere.
		Omniscient is a hash to save gff3 data in a specific 3 levels way: eg Level1: gene, Level2: mRNA, Level3:exon,cds,utr.
		Parser phylosophy: Parse by Parent/child relationship
												ELSE Parse by a common tag	(an attribute value shared by feature that must be grouped together. By default we are using locus_tag but can be set by parameter)
												 ELSE Parse by sequential (mean group features in a bucket, and the bucket change at each level2 feature, and bucket are join in a comon tag at each new L1 feature)

		/!\ Case with only level3 features (i.e rast or some prokka files, sequential will not work as expected. Indeed all features will be the child of only one newly created Parent.
			 To create a parent per feature or group of feature, a comon_tag must be used to regroup them correclty.)

		To resume by priority of way to parse: Parent/child relationship > locus_tag > sequential

		Omniscient data strucure: Ii is a hash to store all the gff feature in 3 levels structures:
			example at level1: $omniscient->{"level1"}{$primary_tag}{$id}=$feature;
			example at level2 and level3: $omniscient->{"levelX"}{$primary_tag}{$parent} = [$feature];
		It is easy then to pare the data structure by record (set of data linked to each other representative of a complex entity, e.g: gene + mRNA + exon + CDS + UTR features) from a top feature to its subfeature.
		When creating a Omniscientwe also return another hash less complex (mRNAGeneLink) allowing to parse the data from level2 to level1 (avoiding to go through the whole omniscient to retrieve this information).
		From Level3 to level 1 it's already possible.
			 example: $mRNAGeneLink->{lc($id)}=$parent;

=head1 AUTHOR

	 Jacques Dainat - jacques.dainat@nbis.se

=cut

#===== TO	DO =====
# When creating a parent check its type from the value in the constant hash

##########################
#		 DEFINE CONSTANT		#
use constant PREFIX_NEW_ID => "nbis"; # used when creating a new ID # old nbis_NEW
use constant PREFIX_ID_L1_NEW => "nbisL1"; # used when creating a new ID for a new Level1 feature # old nbis_noL1id
use constant PREFIX_ID_L2_NEW => "nbisL2"; # used when creating a new ID for a new Level2 feature # old nbis_noL2id

# Comon_tag is used in old gff format and in gtf (with gene_id) to group features together. Priority to comonTag compare to sequential read
# COMONTAG is accessible from the whole file. if a tag has been specified by a user, it is added to this list when slurp_gff3_file_JD is called
my @COMONTAG = ('locus_tag','gene_id');
my $createL3forL2orphan = 1; #case 28

# ====== PURPOSE =======:
# Save in omniscient hash (sorted in a specific way (3 levels)) a whole gff3 file
# Parser phylosophy: Parse by Parent/child ELSE
#						Parse by comon_tag	ELSE
#							Parse by sequential (mean group features in a bucket, and the bucket change at each level2 feature, and bucket are join in a comon tag at each new L1 feature)
#							So if only level3 feature (i.e rast or some prokka files, sequential will not work. A comon_tag must be provided)
# Priority Parent > locus_tag > sequential
# ====== INPUT =======:
# $file => string (file) / list / hash
# $locus_tag => tag to consider for gathering features (in top of the default one)
# $gff_version => Int (if is used, force the parser to use this gff parser instead of guessing)
# $verbose => define the deepth of verbosity
# $no_check is to deactivate sanity check. We can tune the deactivation of check steps using no_check_skip.
# $no_check_skip [list] is to avoid the deactivation of all check with the $no_check parameter. Check steps listed here will not be deactivated.
sub slurp_gff3_file_JD {

	my $start_run = time();
	my $previous_time = undef;
	my %omniscient; #Hast where all the feature will be saved

#	+-----------------------------------------+
#	|							HANDLE ARGUMENTS						|
#	+-----------------------------------------+
	my ($args) = @_	;

	# Check we receive a hash as ref
	if(ref($args) ne 'HASH'){ print "Hash Arguments expected for slurp_gff3_file_JD. Please check the call.\n"; exit;	}

	# Declare all variables and fill them
	my ($file, $gff_version, $locus_tag, $verbose, $no_check, $merge_loci, $no_check_skip, $expose_feature_levels, $log, $debug);
	#first define verbosity
	if( defined($args->{verbose}) ) {$verbose = $args->{verbose};}
		else{ $verbose = 1; } # verbose 0 is quite mode.
	if( defined($args->{log})){
			my $log_name = $args->{log};
			open($log, '>', $log_name  ) or
						dual_print($log, "Can not open $log_name for printing: $!", 1) && die;
			print $log AGAT::Omniscient::get_agat_header(); # print AGAT header
			print $log file_text_line({ string => (strftime "%m/%d/%Y at %Hh%Mm%Ss", localtime),
																  char => " ",
																  extra => "\n"});
	}
	if( defined($args->{debug})) {$debug = $args->{debug};}
	dual_print ($log, surround_text("- Start parsing -",80,"*"), $verbose);
	dual_print($log, file_text_line({ string => "parse options and metadata", char => "-" }), $verbose);
	#Secondly check if expose_feature_levels option
	if( defined($args->{expose_feature_levels})){ $expose_feature_levels = $args->{expose_feature_levels};
										dual_print ($log, "=> Expose feature level json files\n", $verbose );
	} # list of check to skip
	dual_print ($log, "=> Accessing the feature level json files\n", $verbose );
	load_levels( {omniscient => \%omniscient, expose => $expose_feature_levels, verbose => $verbose, log => $log, debug => $debug}); # 	HANDLE feature level
	if( defined($args->{input})) {$file = $args->{input};} 					 else{ dual_print($log, "Input data --input is mandatory when using slurp_gff3_file_JD!"); exit;}
	if( defined($args->{gff_version})) { $gff_version = $args->{gff_version}; } # force using gff parser version
	if( defined($args->{locus_tag})) {
		if( ref($args->{locus_tag}) ne 'ARRAY') {
			@COMONTAG = ($args->{locus_tag});
		}
		else{
			@COMONTAG = @{$args->{locus_tag}};
		}
	} #add a new comon tag to the list if provided.}
		dual_print($log, "=> Attribute used to group features when no Parent/ID relationship exists:\n", $verbose);
		foreach my $comTag (@COMONTAG){
			dual_print($log, "	* $comTag\n", $verbose);
		}
	if( defined($args->{no_check})) { $no_check = $args->{no_check}; dual_print($log, "=> no_check option activated\n", $verbose); } # skip checks
	if( defined($args->{no_check_skip})) {
			$no_check_skip = $args->{no_check_skip}	;

			if( ref($no_check_skip) ne 'ARRAY') {
				$no_check_skip = [$no_check_skip];
			}
			dual_print($log, "=> Check forced:\n", $verbose);
			foreach my $no_check (@$no_check_skip){
				dual_print($log, " * $no_check\n", $verbose);
			}
	}
	else {$no_check_skip = [];} 	 # arrayref of check to skip

	if( defined($args->{merge_loci})) { $merge_loci = $args->{merge_loci}; dual_print($log, "=> merge_loci option activated\n", $verbose); } # activat merge locus option
		else{ $merge_loci = undef;	dual_print($log, "=> merge_loci option deactivated\n", $verbose); }

#	+-----------------------------------------+
#	|	HANDLE GFF HEADER						|
#	+-----------------------------------------+
	my $gff3headerInfo = _check_header($file, $log);

#	+-----------------------------------------+
#	|	HANDLE SOFA (feature-ontology)			|
#	+-----------------------------------------+
	my $ontology = {};
	my $ontology_obj = _handle_ontology($gff3headerInfo, $verbose, $log);
	if($ontology_obj){
		$ontology = create_term_and_id_hash($ontology_obj, $verbose, $log, $debug);
	}

#	+-----------------------------------------+
#	|			HANDLE WARNING 					|
#	+-----------------------------------------+
	my %WARNS;
	my %globalWARNS;
	my $nbWarnLimit = $debug ? undef : 10; # limit number of warning if not debug mode
		local $SIG{__WARN__} = sub {
		my $message = shift;
		my @thematic=split /@/,$message ;

		if($thematic[0] eq "GLOBAL"){ #extract global warning (about feature type in ontology and if AGAT deal with it)
			push @{$globalWARNS{$thematic[1]}}, $thematic[2];
		}
		else{
			# Print  a limited amount of warning
			$WARNS{$thematic[0]}++;
			if($nbWarnLimit){
				if ($WARNS{$thematic[0]} <= $nbWarnLimit){
					dual_print($log, $message, $verbose);
				}
				if($WARNS{$thematic[0]} == $nbWarnLimit){
					dual_print($log, "$thematic[0] ************** Too much WARNING message we skip the next **************\n", $verbose);
				}
			}
			# Print all warning
			else{
				dual_print($log, $message, $verbose);
			}
		}
	};

#	+-------------------------------------------------------------------------+
#	|				HANDLE FEATUTRES PARSING ACCORDING TO TYPE OF INPUTS			|
#	+-------------------------------------------------------------------------+
	my %mRNAGeneLink; #Hast that keep track about link between l2 and l1
	my %duplicate;# Hash to store duplicated feature info
	my %miscCount;# Hash to store any counter. Will be use to create a new uniq ID
	my %uniqID;# Hash to follow up with an uniq identifier every feature
	my %uniqIDtoType; # Hash to follow up with an uniq identifier every feature type
	my %locusTAG; # Hash to follow up the locus tag met
	my %infoSequential;# Hash to store sequential features
	my $locusTAGvalue=undef;
	my $last_l1_f=undef;
	my $last_l2_f=undef;
	my $last_l3_f=undef;
	my $last_f=undef;# last feature handled
	my $lastL1_new =undef; # Bolean to check if last l1 feature is a newly created one. Important to deal with strict sequential

	dual_print($log, file_text_line({ string => "parse features", char => "-" }), $verbose);

	# ============================> ARRAY CASE <============================
	if(ref($file) eq 'ARRAY'){
		 foreach my $feature (@{$file}) {
			($locusTAGvalue, $last_l1_f, $last_l2_f, $last_l3_f, $last_f, $lastL1_new) =
					manage_one_feature($ontology, $feature, \%omniscient, \%mRNAGeneLink, \%duplicate, \%miscCount, \%uniqID, \%uniqIDtoType, \%locusTAG, \%infoSequential, $locusTAGvalue, $last_l1_f, $last_l2_f, $last_l3_f, $last_f, $lastL1_new, $verbose, $log, $debug);
	 		}
	}
	# ============================> HASH CASE <============================
	elsif(ref($file) eq 'HASH'){

		foreach my $level (keys %{$file}){
			# save header if any
			if ($level eq 'other'){
				foreach my $thing( keys %{$file->{'other'} } ){
					$omniscient{'other'}{$thing} = $file->{'other'}{$thing};
				}
				next;
			}
			if ( ref($file->{$level}) eq 'HASH'){ #Header,level1,level2,#level3
				foreach my $tag (keys %{$file->{$level}}){
					foreach my $id (keys %{$file->{$level}{$tag}}){
						if ( ref($file->{$level}{$tag}{$id}) eq 'ARRAY'){ #level2,#level3
							foreach my $feature ( @{$file->{$level}{$tag}{$id} }){
								($locusTAGvalue, $last_l1_f, $last_l2_f, $last_l3_f, $last_f, $lastL1_new) =
										manage_one_feature($ontology, $feature, \%omniscient, \%mRNAGeneLink, \%duplicate, \%miscCount, \%uniqID, \%uniqIDtoType, \%locusTAG, \%infoSequential, $locusTAGvalue, $last_l1_f, $last_l2_f, $last_l3_f, $last_f, $lastL1_new, $verbose, $log, $debug);
							}
						}
						else{ #level1
							my $feature = $file->{$level}{$tag}{$id};
							($locusTAGvalue, $last_l1_f, $last_l2_f, $last_l3_f, $last_f, $lastL1_new) =
									manage_one_feature($ontology, $feature, \%omniscient, \%mRNAGeneLink, \%duplicate, \%miscCount, \%uniqID, \%uniqIDtoType, \%locusTAG, \%infoSequential, $locusTAGvalue, $last_l1_f, $last_l2_f, $last_l3_f, $last_f, $lastL1_new, $verbose, $log, $debug);
						}
					}
				}
			}
			else{ #extra list of feature
				foreach my $feature ( @{$file->{$level}} ) {
					($locusTAGvalue, $last_l1_f, $last_l2_f, $last_l3_f, $last_f, $lastL1_new) =
							manage_one_feature($ontology, $feature, \%omniscient, \%mRNAGeneLink, \%duplicate, \%miscCount, \%uniqID, \%uniqIDtoType, \%locusTAG, \%infoSequential, $locusTAGvalue, $last_l1_f, $last_l2_f, $last_l3_f, $last_f, $lastL1_new, $verbose, $log, $debug);
	 				}
			}
		}
	}
	# ============================> FILE CASE <============================
	else{
		# take care of headers
		my $header = get_header_lines($file, $verbose, $log);
		$omniscient{'other'}{'header'}=$header if $header;

		#GFF format used for parser
		my $format;
		if($gff_version){$format = $gff_version;}
		else{ $format = select_gff_format($file, $verbose, $log);}
		dual_print( $log, "=> GFF version parser used: $format\n", $verbose );
		my $gffio = Bio::Tools::GFF->new(-file => $file, -gff_version => $format);

		#read every lines
		while( my $feature = $gffio->next_feature()) {
			if($format eq "1"){_gff1_corrector($feature, $verbose);} # case where gff1 has been used to parse.... we have to do some attribute manipulations
			($locusTAGvalue, $last_l1_f, $last_l2_f, $last_l3_f, $last_f, $lastL1_new) =
					manage_one_feature($ontology, $feature, \%omniscient, \%mRNAGeneLink, \%duplicate, \%miscCount, \%uniqID, \%uniqIDtoType, \%locusTAG, \%infoSequential, $locusTAGvalue, $last_l1_f, $last_l2_f, $last_l3_f, $last_f, $lastL1_new, $verbose, $log, $debug);
			}

			#close the file
			$gffio->close();
	}

	#------- Inform user about warnings encountered during parsing ---------------
	foreach my $thematic (keys %WARNS){
		my $nbW = $WARNS{$thematic};
		dual_print($log, "$nbW warning messages: $thematic\n", $verbose);
	}

	# Parsing time
	dual_print ($log, surround_text("- End parsing -\ndone in ".(time() - $start_run)." seconds",80,"*","\n"), $verbose );
	$previous_time = time();
	my $check_time = $previous_time;

	#	+-----------------------------------------+
	#	|					 	 CHECK OMNISCIENT				 			|
	#	+-----------------------------------------+

	# -------------------- Mandatory checks --------------------
	my $check_cpt = 1;
	dual_print ($log, surround_text("- Start checks -",80,"*"), $verbose);

	dual_print ($log, file_text_line({ string => "Check$check_cpt: feature types", char => "-" }), $verbose );
	_handle_globalWARNS({ warning => \%globalWARNS, ontology => $ontology, log => $log, type => "ontology", verbose => $verbose });
	_handle_globalWARNS({ warning => \%globalWARNS, ontology => $ontology, log => $log, type => "agat", verbose => $verbose });
	delete $globalWARNS{$_} for (keys %globalWARNS); # re-initialize the hash
	delete $WARNS{$_} for (keys %WARNS); # re-initialize the hash
	dual_print ($log, file_text_line({ string => "	 done in ".(time() - $previous_time)." seconds", char => "-", extra => "\n" }), $verbose );
	$check_cpt++; $previous_time = time();

	#report detected duplicates
	dual_print ($log, file_text_line({ string => "Check$check_cpt: duplicates", char => "-" }), $verbose );
	_check_duplicates($log, \%duplicate, \%omniscient, $verbose);
	dual_print ($log, file_text_line({ string => "	 done in ".(time() - $previous_time)." seconds", char => "-" }), $verbose );
	$check_cpt++; $previous_time = time();

	# -------------------- Extra checks --------------------

	if(! $no_check or	grep( /check_sequential/, @$no_check_skip ) ) {
		#Check sequential if we can fix cases. Hash to be done first, else is risky that we remove orphan L1 feature ... that are not yet linked to a sequential bucket
		dual_print ($log, file_text_line({ string => "Check$check_cpt: sequential bucket", char => "-", prefix => "\n"}), $verbose );
		if( keys %infoSequential ){ #hash is not empty
				_check_sequential($debug, $log, \%infoSequential, \%omniscient, \%miscCount, \%uniqID, \%uniqIDtoType, \%locusTAG, \%mRNAGeneLink, $verbose);
				undef %infoSequential;
		}
		else{
				dual_print ($log, "Nothing to check as sequential !\n", $verbose);
		}
		dual_print ($log, file_text_line({ string => "	 done in ".(time() - $previous_time)." seconds", char => "-" }), $verbose );
		$check_cpt++; $previous_time = time();
	}

	if(! $no_check or	grep( /check_l2_linked_to_l3/, @$no_check_skip ) ) {
			#Check relationship between l3 and l2
			dual_print ($log, file_text_line({ string => "Check$check_cpt: l2 linked to l3", char => "-", prefix => "\n" }), $verbose );
			_check_l2_linked_to_l3($log, \%omniscient, \%mRNAGeneLink, \%miscCount, \%uniqID, \%uniqIDtoType, $verbose, $debug); # When creating L2 missing we create as well L1 if missing too
			dual_print ($log, file_text_line({ string => "	 done in ".(time() - $previous_time)." seconds", char => "-" }), $verbose );
			$check_cpt++; $previous_time = time();
	}

	if(! $no_check or	grep( /check_l1_linked_to_l2/, @$no_check_skip ) ) {
			#Check relationship between mRNA and gene.	/ gene position are checked! If No Level1 we create it !
			dual_print ($log, file_text_line({ string => "Check$check_cpt: l1 linked to l2", char => "-", prefix => "\n" }), $verbose );
			_check_l1_linked_to_l2($log, \%omniscient, \%miscCount, \%uniqID, \%uniqIDtoType, $verbose, $debug);
			dual_print ($log, file_text_line({ string => "	 done in ".(time() - $previous_time)." seconds", char => "-" }), $verbose );
			$check_cpt++; $previous_time = time();
	}

	if(! $no_check or	grep( /remove_orphan_l1$/, @$no_check_skip ) ) {
			#check level1 has subfeature else we remove it
			dual_print ($log, file_text_line({ string => "Check$check_cpt: remove orphan l1", char => "-", prefix => "\n" }), $verbose );
			dual_print ($log, "We remove only those not supposed to be orphan\n");
			_remove_orphan_l1(\%omniscient, \%miscCount, \%uniqID, \%uniqIDtoType, \%mRNAGeneLink, $verbose, $log, $debug); #or fix if level2 is missing (refseq case)
			dual_print ($log, file_text_line({ string => "	 done in ".(time() - $previous_time)." seconds", char => "-" }), $verbose );
			$check_cpt++; $previous_time = time();
	}

	if(! $no_check or	grep( /check_exons/, $no_check_skip ) ) {
			#Check relationship L3 feature, exons have to be defined... / mRNA position are checked!
			dual_print ($log, file_text_line({ string => "Check$check_cpt: check exons", char => "-", prefix => "\n" }), $verbose );
			_check_exons($debug, $log, \%omniscient, \%mRNAGeneLink, \%miscCount, \%uniqID,	\%uniqIDtoType, $verbose);
			dual_print ($log, file_text_line({ string => "	 done in ".(time() - $previous_time)." seconds", char => "-"}), $verbose );
			$check_cpt++; $previous_time = time();
	}

	if(! $no_check or	grep( /check_utrs/, @$no_check_skip ) ) {
			#Check relationship L3 feature, exons have to be defined... / mRNA position are checked!
			dual_print ($log, file_text_line({ string => "Check$check_cpt: check utrs", char => "-", prefix => "\n" }), $verbose );
			_check_utrs($debug, $log, \%omniscient, \%mRNAGeneLink, \%miscCount, \%uniqID,	\%uniqIDtoType, $verbose);
			dual_print ($log, file_text_line({ string => "	 done in ".(time() - $previous_time)." seconds", char => "-" }), $verbose );
			$check_cpt++; $previous_time = time();
	}

	if(! $no_check or	grep( /check_all_level2_positions/, @$no_check_skip ) ) {
		# Check rna positions compared to its l2 features
		dual_print ($log, file_text_line({ string => "Check$check_cpt: all level2 locations", char => "-", prefix => "\n" }), $verbose );
		check_all_level2_positions( { omniscient => \%omniscient, verbose => $verbose, log => $log } );
		dual_print ($log, file_text_line({ string => "	 done in ".(time() - $previous_time)." seconds", char => "-" }), $verbose );
		$check_cpt++; $previous_time = time();
	}

	if(! $no_check or	grep( /check_all_level1_positions/, @$no_check_skip ) ) {
		# Check gene positions compared to its l2 features
		dual_print ($log, file_text_line({ string => "Check$check_cpt: all level1 locations", char => "-", prefix => "\n" }), $verbose );
		check_all_level1_positions( { omniscient => \%omniscient, verbose => $verbose, log => $log } );
		dual_print ($log, file_text_line({ string => "	 done in ".(time() - $previous_time)." seconds", char => "-"}), $verbose );
		$check_cpt++; $previous_time = time();
	}

	#check loci names (when overlap should be the same if type is the same)
	if ( $merge_loci ){
		# Better probably to keep it before check 10 anyway
		dual_print ($log, file_text_line({ string => "Check$check_cpt: merge overlaping features into same locus", char => "-", prefix => "\n" }), $verbose );
		_merge_overlap_features($log, \%omniscient, \%mRNAGeneLink, $verbose);
		dual_print ($log, file_text_line({ string => "	 done in ".(time() - $previous_time)." seconds", char => "-" }), $verbose );
		$check_cpt++; $previous_time = time();
	}

	if(! $no_check or	grep( /check_identical_isoforms/, @$no_check_skip ) ) {
		#check identical isoforms
		dual_print ($log, file_text_line({ string => "Check$check_cpt: remove identical isoforms", char => "-", prefix => "\n"}), $verbose );
		_check_identical_isoforms($log, \%omniscient, \%mRNAGeneLink, $verbose);
		dual_print ($log, file_text_line({ string => "	 done in ".(time() - $previous_time)." seconds", char => "-" }), $verbose );
		$check_cpt++; $previous_time = time();
	}

	#------- Inform user about warnings encountered during checking ---------------
	foreach my $thematic (keys %WARNS){
		my $nbW = $WARNS{$thematic};
		if($nbW > $nbWarnLimit){
			dual_print($log, "$nbW warning messages: $thematic\n", $verbose);
		}
	}

	if(! $no_check ){
		dual_print ($log, surround_text("- End checks -\ndone in ".(time() - $check_time)." seconds",80,"*","\n"), $verbose );
	}

 dual_print ($log, "=> OmniscientI total time: ".(time() - $start_run)." seconds\n", $verbose );

	#return
	return \%omniscient, \%mRNAGeneLink	;
}

##==============================================================================
##==============================================================================

# ====== PURPOSE =======:
# The method read a gff3 feature, Check for the sanity according to what will has been read before, and the whole features that will be read.
# Designed to be used within a loop going through a big amout of feature
# ====== INPUT =======:
# $feature => gff feature object
# 		example: scaffold625	maker	exon	341518	341628	.	+	.	ID=CLUHART00000008717:exon:1406;Parent=CLUHART00000008717
# $omniscient => hash to store all the gff feature in 3 levels structures
# 		example at level1: $omniscient->{"level1"}{$primary_tag}{$id}=$feature;
# 		example at other levels: $omniscient->{"levelX"}{$primary_tag}{$parent}= [$feature];
# $mRNAGeneLink => hash to keep track link from L2 to l1 (avoid to go through the whole omniscient to retrieve this information)
# 		example: $mRNAGeneLink->{lc($id)}=$parent;
# $duplicate => hash to keep track of duplicates found
# 		example: duplicate->{$level}{$primary_tag}{$id} = [$feature];
# $miscCount => hash that contains a counter for each feature type. It is used to create uniq ID
# 		example: $miscCount->{$primary_tag}++;
# $uniqID => hash of uniqID (UniqID link to the original ID)
# 		example: $uniqID->{$uID}=$id;
# $uniqIDtoType => hash to keep track about the feature type linked to the uniqID (useful for handling SPREADFEATURES)
# 		example: $uniqIDtoType->{$uID}=$primary_tag;
# $locusTAG_uniq => hash of comon tag found when reading the grouped features sequentialy
# 		example: $locusTAG_uniq->{'level1'}{$id}=$id;
# $infoSequential => hash that contains features grouped together in a sequential order. (Useful as example when no Parent tag and/or locus tag missing)
# 		structure at level1: $infoSequential->{$id}{'level1'}=$id;
# 		structure at other levels: $infoSequential->{$locusTAGvalue}{$l2_id}{'level3'}}, [$feature1,$feature2] ;
# $last_locusTAGvalue => String: Last locus tag that has been met when parsing the file (If no locus tag found it will be the last feature L1 ID)
# 		example: CLUHARG00000008717
# $last_l1_f => String: Last l1 feature that has been met when parsing the file
# 		example: scaffold625	maker	gene	341518	341628	.	+	.	ID=CLUHARG00000008717
# $last_l2_f => String: Last L2 feature that has been met when parsing the file
# 		example: scaffold625	maker	mRNA	341518	341628	.	+	.	ID=CLUHART00000008717;Parent=CLUHARG00000008717
# $last_l3_f => String: Last L3 feature that has been met when parsing the file
# 		example: scaffold625	maker	exon	341518	341628	.	+	.	ID=CLUHART00000008717:exon:1406;Parent=CLUHART00000008717
# $verbose => INT: Verbose level. Bigger the value is, deeper the information sould be. 0 = quiet
# 		example: 2
# ====== OUTPUT======= : Omniscient Hash
sub manage_one_feature{

		my ($ontology, $feature, $omniscient, $mRNAGeneLink, $duplicate, $miscCount,
		$uniqID, $uniqIDtoType, $locusTAG_uniq, $infoSequential, $last_locusTAGvalue,
		$last_l1_f, $last_l2_f, $last_l3_f, $last_f, $lastL1_new, $verbose, $log, $debug)=@_;

		my $seq_id = $feature->seq_id;					#col1
		my $source_tag = lc($feature->source_tag);		#col2
		my $primary_tag = lc($feature->primary_tag);	#col3
		my $start = $feature->start;					#col4
		my $end = $feature->end;						#col5
		my $score = $feature->score;					#col6
		my $strand = $feature->strand;					#col7
		my $frame = $feature->frame;					#col8
		# Attribute => tag=value tuples								#col9
		my $id= undef;
		my $parent= undef;
		my $locusTAGvalue=undef;

#		+-------------------------------+
#		|	 CKECK SEQUENCE ONTOLOGY			|
#		+-------------------------------+
		# We have ontology values (otherwise we inform the user in _handle_globalWARNS)
		if( keys %{$ontology} ){
			# The tag is not part of the ontology let's save this info
			if(! exists_keys($ontology, ($primary_tag) ) ) {
				warn "GLOBAL@"."ontology1@".$primary_tag."@";
			}
		}

#	+----------------------------------------------------------------------------+
#							MANAGE LEVEL1 => feature WITHOUT parent
#	+----------------------------------------------------------------------------+
		if( get_level($omniscient, $feature) eq 'level1' ) {

				##########
				# Deal with standalone and topfeature that do not expect children
				if ($omniscient->{'other'}{'level'}{'level1'}{$primary_tag} eq 'standalone' or
							$omniscient->{'other'}{'level'}{'level1'}{$primary_tag} eq 'topfeature'){

					$id = lc(_check_uniq_id($omniscient, $miscCount, $uniqID, $uniqIDtoType, $feature));
					if(! _it_is_duplication($duplicate, $omniscient, $uniqID, $feature)){
						$omniscient->{"level1"}{$primary_tag}{$id}=$feature;
						dual_print ($log, "::::::::::0Push-L1-omniscient level1 || $primary_tag || $id = ".$feature->gff_string()."\n", $verbose) if ($debug);
					}
					return $last_locusTAGvalue, $last_l1_f, $last_l2_f, $last_l3_f, $last_l1_f, $lastL1_new;
				}

				##########
				# get ID #
				$id = lc(_check_uniq_id($omniscient, $miscCount, $uniqID, $uniqIDtoType, $feature));
				_save_common_tag_value_top_feature($feature, $locusTAG_uniq, 'level1');

				#####################
				# Ckeck duplication #
				if(! _it_is_duplication($duplicate, $omniscient, $uniqID, $feature)){

						################
						# Save feature #
						$last_l1_f = $feature;
						dual_print ($log, "::::::::::Push-L1-omniscient level1 || $primary_tag || $id = ".$feature->gff_string()."\n", $verbose) if ($debug);
						$omniscient->{"level1"}{$primary_tag}{$id}=$feature;
 				 		$locusTAG_uniq->{'level1'}{$id}=$id;

 				 		#############
 				 		# COMON TAG #
						$locusTAGvalue =_get_comon_tag_value($feature, $locusTAG_uniq, 'level1');

						if($locusTAGvalue){
								dual_print ($log, "::::::::::Push-L1-sequential $locusTAGvalue || level1 == $id\n", $verbose) if ($debug);
								$locusTAG_uniq->{'level1'}{$locusTAGvalue}=$id;
								$infoSequential->{$id}{'level1'}=$id;
						}
						else{
								$locusTAGvalue=$id;
						}

						#################
						#reinitialization
						$last_l2_f=undef; #opposite to I have a comon tag
						$last_l3_f=undef;
				}
				return $id, $last_l1_f, $last_l2_f, $last_l3_f, $last_l1_f, $lastL1_new;
		 }

# +----------------------------------------------------------------------------+
#					MANAGE LEVEL2 => feature WITH child and WITH parent
# +----------------------------------------------------------------------------+
		elsif ( get_level($omniscient, $feature) eq 'level2' ) {

				#reinitialization
				$last_l3_f=undef;

				##########
				# get ID #
				$id = _check_uniq_id($omniscient, $miscCount, $uniqID, $uniqIDtoType, $feature);

				##############
				# get Parent #
				#GFF case
				if($feature->has_tag('Parent')){
							$parent = lc($feature->_tag_value('Parent'));
							$locusTAGvalue=$parent;
							_save_common_tag_value_top_feature($feature, $locusTAG_uniq, 'level2');
				}

				#GTF case
				elsif($feature->has_tag('gene_id') ){
						$parent = lc($feature->_tag_value('gene_id'));
						create_or_replace_tag($feature,'Parent',$feature->_tag_value('gene_id')); #modify Parent To keep only one
						$locusTAGvalue=$parent;
						_save_common_tag_value_top_feature($feature, $locusTAG_uniq, 'level2');
				}
				else{
						warn "WARNING level2: No Parent attribute found @ for the feature: ".$feature->gff_string()."\n";

						#################
		 				# COMON TAG PART1
						$locusTAGvalue =_get_comon_tag_value( $feature, $locusTAG_uniq, 'level2');


						######################
						# NEED THE LEVEL1 ID #
						my $l1_ID="";
						# If I don't have a last_l1_f I create one. The Id can be used as comonTag.
						# The feature can also be used later if comon_tag was existing, but mising for one of the feature.
						# If we have comon tag, we check we are changing from the previous one before to create a new level1 feature.
						# It's to deal with potential level2 (like mRNA isoforms).
						if(! $last_l1_f or ($locusTAGvalue and ($locusTAGvalue ne $last_locusTAGvalue) ) ){
								dual_print ($log, "create L1 feature\n", $verbose) if ($debug);
								$l1_ID = _create_ID($miscCount, $uniqID, $uniqIDtoType, $primary_tag, $id, PREFIX_ID_L1_NEW);
								$last_l1_f = clone($feature);
								create_or_replace_tag($last_l1_f,'ID',$l1_ID); #modify Parent To keep only one
								$last_l1_f->primary_tag('gene');
								$lastL1_new = 1;
						}
						else{ # case where previous level1 exists
								# Stricly sequential at level2 feature. We create a new L1 at every L2 met except if two L2 are in a row
								if ( ( $lastL1_new and not exists_keys( $omniscient, ('other', 'level', 'level2', $last_f->primary_tag ) ) ) and
										(!$locusTAGvalue or ($locusTAGvalue ne $last_locusTAGvalue) ) ){ # if previous L1 newly created and last feature is not f2 (if several f2 in a row we attach them to the same newly created l1 feature)

										dual_print ($log, "create L1 feature stritcly\n", $verbose ) if ($debug);
										$l1_ID = _create_ID($miscCount, $uniqID, $uniqIDtoType, $primary_tag, $id, PREFIX_ID_L1_NEW);
										$last_l1_f = clone($feature);
										create_or_replace_tag($last_l1_f,'ID',$l1_ID); #modify Parent To keep only one
										$last_l1_f->primary_tag('gene');
										$lastL1_new = 1;
										$lastL1_new = 1;
								}
								else{
										dual_print ($log, "take last L1 feature\n", $verbose ) if ($debug);
										$l1_ID=$last_l1_f->_tag_value('ID');
										$lastL1_new = undef;
								}
						}
						create_or_replace_tag($feature,'Parent',$l1_ID); #modify Parent To keep only one

						#################
		 				# COMON TAG PART2
						if($locusTAGvalue){ #Previous Level up feature had a comon tag
								dual_print ($log, "::::::::::Push-L2-Sequential-1 $locusTAGvalue || ".lc($id)." || level2 == ".$feature->gff_string."\n", $verbose ) if ($debug);
								$infoSequential->{$locusTAGvalue}{$id}{'level2'} = $feature ;

								return $locusTAGvalue, $last_l1_f, $feature, $last_l3_f, $feature, $lastL1_new;								#### STOP HERE AND RETURN
						}
						else{

								dual_print ($log, "::::::::::Push-L2-omniscient-2: level2 || ".$primary_tag." || ".lc($l1_ID)." == ".$feature->gff_string."\n", $verbose ) if ($debug);
						 		push (@{$omniscient->{"level2"}{$primary_tag}{lc($l1_ID)}}, $feature);

								# keep track of link between level2->leve1 #
								if (! exists ($mRNAGeneLink->{lc($id)})){
										$mRNAGeneLink->{lc($id)}=$l1_ID;
				 				}

								return lc($l1_ID) , $last_l1_f, $feature, $last_l3_f, $feature, $lastL1_new;								#### STOP HERE AND RETURN
						}
				}

				#####################
				# Ckeck duplication #
				if(! _it_is_duplication($duplicate, $omniscient, $uniqID, $feature)){

						############################################
						# keep track of link between level2->leve1 #
						if (! exists ($mRNAGeneLink->{lc($id)})){
								$mRNAGeneLink->{lc($id)}=$parent;
				 		}

				 		####################
				 		# SAVE THE FEATURE #
				 		dual_print ($log, "::::::::::Push-L2-omniscient-3 level2 || $primary_tag || $parent == ".$feature->gff_string."\n", $verbose) if ($debug);
						push (@{$omniscient->{"level2"}{$primary_tag}{lc($parent)}}, $feature);
				}
				return $last_locusTAGvalue, $last_l1_f, $feature, $last_l3_f, $feature, $lastL1_new;
		}

# +----------------------------------------------------------------------------+
#									MANAGE LEVEL3 => feature WITHOUT child
# +----------------------------------------------------------------------------+
		elsif ( get_level($omniscient, $feature) eq 'level3' ){

				# get ID #
				$id = _check_uniq_id($omniscient, $miscCount, $uniqID, $uniqIDtoType, $feature);

				#		+-------------------------------+
				#		|					GET PARENT L3					 |
				#		+-------------------------------+
				my @parentList;

				###################
				#	 GFF case			#
				if($feature->has_tag('Parent')){
						_fix_parent_attribute_when_id_l1_l2_identical($feature, $last_l1_f, $last_l2_f, $uniqID, 'Parent', $verbose); # must be before @parentList
						@parentList = $feature->get_tag_values('Parent');
						$locusTAGvalue=$last_locusTAGvalue;
						_save_common_tag_value_top_feature($feature, $locusTAG_uniq, 'level3');
				}

				###################
				#		GTF case		#
				elsif($feature->has_tag('transcript_id') ){
						_fix_parent_attribute_when_id_l1_l2_identical($feature, $last_l1_f, $last_l2_f, $uniqID, 'transcript_id', $verbose); # must be before @parentList
						@parentList = $feature->get_tag_values('transcript_id');
						create_or_replace_tag($feature,'Parent',$feature->_tag_value('transcript_id')); #modify Parent To keep only one
						$locusTAGvalue=$last_locusTAGvalue;
						_save_common_tag_value_top_feature($feature, $locusTAG_uniq, 'level3');
				}

				################### In that case we create a uniq parentID to create a proper omniscient structure.
				# No parent case	# But the feature itself stay intact without parentID.
				else{
						my $play_this_game=1;
						warn "WARNING level3: No Parent attribute found @ for the feature: ".$feature->gff_string()."\n";

						#################
						# COMON TAG PART1
						$locusTAGvalue = _get_comon_tag_value( $feature, $locusTAG_uniq, 'level3');

						######################
						# NEED THE LEVEL2 ID #
						my $l2_id="";

						#To keep track of locus tag that has been spread over the file, and a piece is found later
						my $skip_last_l2=undef;
						if($last_l2_f and $locusTAGvalue){

								if(exists_keys ($locusTAG_uniq, ('linkl2l1', lc($last_l2_f->_tag_value('ID') ) ) ) ){
										if (lc($locusTAG_uniq->{'linkl2l1'}{lc( $last_l2_f->_tag_value('ID') )}) ne lc($locusTAGvalue)){
												$skip_last_l2=1;
												dual_print ($log, "skip last l2\n", $verbose) if ( $debug );
										}
								}
						}

						# Just to avoid to have parent undef in case there is no parent feature define for the last_l2_f
						my $parent_of_last_l2 = "@@@@";
						if($last_l2_f and $last_l2_f->has_tag('Parent')){ $parent_of_last_l2 = lc($last_l2_f->_tag_value('Parent')); }

						# case where No level2 feature defined yet - I will need a bucketL2
						# OR common tag changed (= level1/level2 different) so we have to create a new level2 tag
						# but only if the last_comon tag is different as the parent of the last_l2_f
						# (In that case we can use the last L2 feature. It was missing the comon tag in it).
						if(! $last_l2_f or
							($locusTAGvalue and ($locusTAGvalue ne $last_locusTAGvalue)
							and $last_locusTAGvalue ne $parent_of_last_l2 or $skip_last_l2) ){
								dual_print ($log, "Come in the complex case L3!!!\n", $verbose) if ($debug);
								#######################
								# Change referentiel => based on the last L2 link to this locus
								#######################

								# case were locus already met before (feature are spread within the file), we link the L3 to the last l2 of this locus.
								if( exists_keys($locusTAG_uniq, ('level2', $locusTAGvalue) ) ){
										dual_print ($log, "Complex case L3 1 !!!\n", $verbose) if ($debug);
										$last_l2_f = @{$locusTAG_uniq->{'level2'}{$locusTAGvalue}}[$#{$locusTAG_uniq->{'level2'}{$locusTAGvalue}}];
										$l2_id = $last_l2_f->_tag_value('ID');
										foreach my $tag_l1 (keys %{$omniscient->{'level1'}}){
												if( exists_keys($omniscient,('level1', $tag_l1, $locusTAG_uniq->{'linkl2l1'}{lc($l2_id)}))){
														$last_l1_f = $omniscient->{'level1'}{$tag_l1}{$locusTAG_uniq->{'linkl2l1'}{lc($l2_id)}};
												}
										}
								}
								else{
										dual_print ($log, "Complex case L3 2 !!!\n", $verbose) if ($debug);
										# case we can catch parent from previous feature L3
										if($last_l3_f and $last_l3_f->has_tag('Parent')){ # Need to not be the first feature because we need a previous feature

												dual_print ($log, "Complex case 2.1 !!!\n", $verbose) if ($debug);
												my $previousL3 =_get_comon_tag_value( $last_l3_f, $locusTAG_uniq, 'level3');

												if ($locusTAGvalue and $previousL3 and ($locusTAGvalue eq $previousL3)){

														dual_print ($log, "Complex case 2.1.1 !!!\n", $verbose) if ($debug);
														$l2_id = $last_l3_f->_tag_value("Parent");
														create_or_replace_tag($feature,'Parent', $l2_id);
														push @parentList, $l2_id;
														$play_this_game=undef; # Only place where we skip this game
												}
										}
										if ($play_this_game){
												dual_print ($log, "Complex case 2.2 !!!\n", $verbose) if ($debug);
												$l2_id = _create_ID($miscCount, $uniqID, $uniqIDtoType, $primary_tag, $id, PREFIX_ID_L2_NEW);
												$last_l2_f = clone($feature);
												create_or_replace_tag($last_l2_f,'ID',$l2_id); #modify Parent To keep only one
												$last_l2_f->primary_tag('RNA');
										}
								}
						}
						# case where previous level2 exists
						else{
							dual_print ($log, "case where previous level2 exists\n", $verbose) if ($debug);
							$l2_id=$last_l2_f->_tag_value('ID');
						}

						# Let's play the no parent case. We will return from that part of the code
						# We don't play this game only if we decided finally to take the same parent
						# value as previous feature
						if ($play_this_game){

								create_or_replace_tag($feature,'Parent',$l2_id); #modify Parent To keep only one

								#############
				 				# COMON TAG	Part2
								if($locusTAGvalue){ #Previous Level up feature had a comon tag
										dual_print ($log, "::::::::::Push-L3-sequential-1 $locusTAGvalue || ".lc($l2_id)." || level3 == ".$feature->gff_string."\n", $verbose) if ($debug);
										#if($feature->_tag_value("Parent") eq ""){exit;}
										### TAKE LAST L2 of the locus tag iF exist !
										push( @{$infoSequential->{$locusTAGvalue}{$l2_id}{'level3'}}, $feature );
										return $locusTAGvalue, $last_l1_f, $last_l2_f, $feature, $feature, $lastL1_new;								#### STOP HERE AND RETURN
								}
								else{# No comon tag found
										######################
										# NEED THE LEVEL1 ID #
										if(!$last_l1_f and $last_l3_f){ #particular case : Two l3 that follow each other, but first one has locus_tag but not the second
												dual_print ($log, "::::::::::Push-L3-sequential-2 $last_locusTAGvalue || ".lc($l2_id)." || level3 == ".$feature->gff_string."\n", $verbose) if ($debug);
												push( @{$infoSequential->{$last_locusTAGvalue}{$l2_id}{'level3'}}, $feature );
												return $last_locusTAGvalue, $last_l1_f, $last_l2_f, $feature, $feature, $lastL1_new;
										}
										else{
												my $l1_id="";
												if($last_l1_f){ # case where previous level1 exists
														$l1_id=$last_l1_f->_tag_value('ID');
												}
												else{ # case where No level1 feature defined yet - I will need a bucketL1
														$l1_id = _create_ID($miscCount, $uniqID, $uniqIDtoType, $primary_tag, $id, PREFIX_ID_L1_NEW);
														$last_l1_f = clone($feature);
														create_or_replace_tag($last_l1_f,'ID',$l1_id); #modify Parent To keep only one
														$last_l1_f->primary_tag('gene');
												}

												#push( @{$infoSequential->{lc($l1_id)}{lc($l2_id)}{'level3'}}, $feature );
												dual_print ($log, "::::::::::Push-L3-omiscient-3: level3 ".$primary_tag." || ".lc($l2_id)." == ".$feature->gff_string."\n", $verbose) if ($debug);
												push (@{$omniscient->{"level3"}{$primary_tag}{lc($l2_id)}}, $feature);
												return $l2_id, $last_l1_f, $last_l2_f, $feature, $feature, $lastL1_new;								#### STOP HERE AND RETURN
										}
								}
						}
				}
			# END No parent case	#
			#######################

			#		+--------------------------------------+
			#		|					HANDLE PARENT(S) L3					 |
			#		+--------------------------------------+
			#	Save feature and check duplicates
			# (treat also cases where there is multiple parent. => In that case we expand to create a uniq feature for each)
			my $cptParent=0; # to check if it is a multiple parent case
					my $allParent = scalar @parentList;
					foreach my $parent (@parentList){ # first feature level3 with this primary_tag linked to the level2 feature
				$cptParent++;

				#Level3 key doesn't exist
				if(! exists_keys($omniscient,('level3',$primary_tag,lc($parent)))){

					# It is a multiple parent case
					if($allParent > 1){

						# Not the first parent, we have to clone the feature !!
						if($cptParent > 1){

							my $feature_clone=clone($feature);
							create_or_replace_tag($feature_clone,'Parent',$parent); #modify Parent To keep only one
							_check_uniq_id($omniscient, $miscCount, $uniqID, $uniqIDtoType, $feature_clone); #Will change the ID if needed

							dual_print ($log, "::::::::::Push-L3-omniscient-4 level3 || $primary_tag || ".lc($parent)." == ".$feature->gff_string."\n", $verbose) if ($debug);
							push (@{$omniscient->{"level3"}{$primary_tag}{lc($parent)}}, $feature_clone);
						}

						# It is the first parent. Do not clone the feature
						else{
							create_or_replace_tag($feature,'Parent',$parent); #modify Parent To keep only one
							dual_print ($log, "::::::::::Push-L3-omniscient-5 level3 || $primary_tag || ".lc($parent)." == ".$feature->gff_string."\n", $verbose) if ($debug);
							push (@{$omniscient->{"level3"}{$primary_tag}{lc($parent)}}, $feature);
						}
					}
					else{ #the simpliest case. One parent only
						dual_print ($log, "::::::::::Push-L3-omniscient-6 level3 || $primary_tag || ".lc($parent)." == ".$feature->gff_string."\n", $verbose) if ($debug);
						push (@{$omniscient->{"level3"}{$primary_tag}{lc($parent)}}, $feature);
					}
				}

				#Level3 key exists
				else{

					# It is a multiple parent case # Not the first parent, we have to clone the feature !!
					if($cptParent > 1){ #several parent, and we are not looking the first one

						my $feature_clone=clone($feature);
						create_or_replace_tag($feature_clone,'Parent',$parent); #modify Parent To keep only one
						_check_uniq_id($omniscient, $miscCount, $uniqID, $uniqIDtoType, $feature_clone); #Will change the ID if needed

						if( ! _it_is_duplication($duplicate, $omniscient, $uniqID, $feature_clone) ){
							dual_print ($log, "::::::::::Push-L3-omniscient-8 level3 || $primary_tag || ".lc($parent)." == ".$feature_clone->gff_string."\n", $verbose) if ($debug);
							push (@{$omniscient->{"level3"}{$primary_tag}{lc($parent)}}, $feature_clone);
						}

					}
					elsif($allParent > 1){ # It is a multiple parent case #several parent, but we are looking the first one

						# It is the first parent. Do not clone the feature
						create_or_replace_tag($feature,'Parent',$parent); #modify Parent To keep only one
						if( ! _it_is_duplication($duplicate, $omniscient, $uniqID, $feature) ){
							dual_print ($log, "::::::::::Push-L3-omniscient-9 level3 || $primary_tag || ".lc($parent)." == ".$feature->gff_string."\n", $verbose) if ($debug);
							push (@{$omniscient->{"level3"}{$primary_tag}{lc($parent)}}, $feature);
						}

					}
					#It is not a duplicated feature => save it in omniscient
					elsif( ! _it_is_duplication($duplicate, $omniscient, $uniqID, $feature) ){
 						#the simpliest case. One parent only
						dual_print ($log, "::::::::::Push-L3-omniscient-10 level3 || $primary_tag || ".lc($parent)." == ".$feature->gff_string."\n", $verbose) if ($debug);
						push (@{$omniscient->{"level3"}{$primary_tag}{lc($parent)}}, $feature);
					}
				}
					}
					return $last_locusTAGvalue, $last_l1_f, $last_l2_f, $feature, $feature, $lastL1_new;
				}

# +----------------------------------------------------------------------------+
# |MANAGE THE REST => feature UNKNOWN | # FEATURE NOT DEFINE IN ANY OF THE 3 LEVELS YET
# +----------------------------------------------------------------------------+
		else{
				warn "gff3 reader warning: primary_tag error @ ".$primary_tag." still not taken in account!".
				" Please modify the json files to define the feature in one of the levels.\n";
				warn "GLOBAL@"."parser1@".$primary_tag."@";
				return $locusTAGvalue, $last_l1_f, $last_l2_f, $last_l3_f, $feature, $lastL1_new;
		}

			print "Congratulation ! Read this line is not normal !! Please contact the developer !!!\n";exit;
}

##==============================================================================
##==============================================================================

# case where l1 and l2 feature have same ID, l2 ID has been modified,
# consequently parent of L3 has to be updated
sub _fix_parent_attribute_when_id_l1_l2_identical{
	my ($feature, $last_l1_f, $last_l2_f, $uniqID, $tag_parent, $verbose) = @_ ;

  if (! $last_l1_f or ! $last_l2_f){return;}
	my $last_l1_id = $last_l1_f->_tag_value('ID');
	my $last_l2_id = $last_l2_f->_tag_value('ID');
	if (! exists_keys ($uniqID, ( $last_l2_id ) ) ){ return;} # case the previous L2 was a duplicate
	my $previous_l2_original_id = lc( $uniqID->{ lc($last_l2_id) } );

	if ( lc($last_l1_id) eq lc($previous_l2_original_id) ){
		print "l1 and l2 had same ID, not normal but I will fix it.\n" if ($verbose > 2);
		print "previous_l2_original_id $previous_l2_original_id\n";
		print $last_l1_f->gff_string."\n";
		print $last_l2_f->gff_string."\n";
		print $feature->gff_string."\n";
		create_or_replace_tag($feature, $tag_parent, $last_l2_f->_tag_value('ID') );
	}
}

# /!\ $feature must have a parent if not level1
# Keep track to recover from sequential locus tag share whith feature saved in Omniscient
sub _save_common_tag_value_top_feature{
	my ($feature, $locusTAG_uniq, $level)=@_;

	my $locusName=undef;

	foreach my $tag (@COMONTAG){

			#check if we have the tag
		if($feature->has_tag($tag)){
			$locusName=lc($feature->_tag_value($tag)); #get the value

			if ( !( exists_keys ( $locusTAG_uniq, ('topfeature', $locusName, 'level1') ) ) and ($level eq 'level2') ) {
				$locusTAG_uniq->{'topfeature'}{$locusName}{'level1'}{'ID'} = lc($feature->_tag_value('ID'));
				last;
			}

			if ( !( exists_keys ( $locusTAG_uniq, ('topfeature', $locusName, 'level2') ) ) and ($level eq 'level2') ) {
					$locusTAG_uniq->{'topfeature'}{$locusName}{'level2'} = [lc($feature->_tag_value('ID')), lc($feature->_tag_value('Parent'))];
					last;
			}

			if ( !( exists_keys ( $locusTAG_uniq, ('topfeature', $locusName, 'level3') ) ) and ($level eq 'level3') ) {
				$locusTAG_uniq->{'topfeature'}{$locusName}{'level3'} = [lc($feature->_tag_value('ID')), lc($feature->_tag_value('Parent'))];
				last;
			}
		}
	}
}

#check if the comom tag is present among the attributes
# return lower case value
sub _get_comon_tag_value{
	my ( $feature, $locusTAG_uniq, $level)=@_;

	my $locusName=undef;
	foreach my $tag (@COMONTAG){

	 		#check if we have the tag
		if($feature->has_tag($tag)){
				$locusName=lc($feature->_tag_value($tag)); #get the value

				if(exists_keys ($locusTAG_uniq, ('level1',$locusName) ) ){
					$locusName = $locusTAG_uniq->{'level1'}{$locusName};
					last;
				}
				else{
					$locusTAG_uniq->{'level1'}{$locusName} = $locusName; #save it
					last;
			}
		}
	}

	if($level eq 'level2' and $locusName){
		if(! exists_keys ($locusTAG_uniq, ('level2',$locusName, lc($feature->_tag_value('ID'))) ) ){
			push @{$locusTAG_uniq->{'level2'}{$locusName}}, $feature;
			$locusTAG_uniq->{'linkl2l1'}{lc($feature->_tag_value('ID'))} =	$locusName;
		}
	}

	# In case where no parent, no comon tag, and no sequential, we cannot deal at all with it !!!!
	if(! $locusName and $level ne 'level1'){
		warn "WARNING gff3 reader: Hmmm, be aware that your feature doesn't contain any Parent and locus tag.".
		" No worries, we will handle it by considering it as striclty sequential.".
		" If you disagree, please provide an ID or a comon tag by locus.".
		" @ the feature is:\n".$feature->gff_string()."\n";
	}

	return $locusName;
}

#feature is not yet saved in omniscient !
sub _it_is_duplication{
	my ($duplicate, $omniscient, $uniqID, $feature)=@_;

	my $is_dupli=undef;
	my $potentialList=undef;

	my $level = get_level($omniscient, $feature);
	my $primary_tag = lc($feature->primary_tag);

	my $id = $uniqID->{ lc($feature->_tag_value('ID') ) }; # check the original ID

	if($level eq "level1"){
		if(! exists_keys($omniscient,($level, $primary_tag, lc($id) ))){
			return $is_dupli; #return is not a dupli
		}
		else{
			$potentialList=$omniscient->{$level}{$primary_tag}{lc($id)}; #push the feature L1 in potentialList
		}
	}
	else{ #feature l2 or l3

		my @parent = $feature->get_tag_values('Parent');
		foreach my $one_parent_ID (@parent){

			my $one_parent_uID = $one_parent_ID; # In case where the parent have not yet been processed, we cannot have his uID, we will check the current ID
			if ( exists_keys( $uniqID, ( lc($one_parent_ID) ) ) ){
				$one_parent_uID = lc($uniqID->{ lc($one_parent_ID) } ); # check the original ID
			}

			foreach my $primary_tag ( keys %{$omniscient->{$level}} ){
				if (exists_keys($omniscient,($level, $primary_tag, $one_parent_uID))){
					push @{$potentialList}, @{$omniscient->{$level}{$primary_tag}{$one_parent_uID}};
				}
			}
		}
		if(! $potentialList){ #potential list empty
				return $is_dupli; #return is not a dupli
		}
	}


	#check if list is a feature or a reference of a list of feature
	my @list_feature; # will be a list of feature.

	if(ref($potentialList) eq 'ARRAY'){
		@list_feature=@$potentialList; #it's a reference of a list of feature
	}
	else{push (@list_feature, $potentialList);} # it's a feature

	#### PREPARE THE SENTENCE TO CHECK
	my $current_string=_create_comparison_string($omniscient, $uniqID, $feature);

	#Check all the level2 list element
	foreach my $feature_in_omniscient ( @list_feature ){

		my $string=_create_comparison_string($omniscient, $uniqID, $feature_in_omniscient);
		if($current_string eq $string){
			$is_dupli=1;
			push (@{$duplicate->{$level}{$primary_tag}{$id}}, $feature);
			delete $uniqID->{ lc($feature->_tag_value('ID') ) }; # delete uniq ID that has been created for nothing
			last;
		}
	}

	if(! $is_dupli and $level eq "level1" and $omniscient->{"level1"}{$primary_tag}{$id}){
		warn "WARNING level1: This feature level1 is not a duplicate".
		" (different than others feature level1) but has an ID already used.".
		" We do not deal with that currently. @ the feature is: ".$feature->gff_string()."\n";
		# Indeed If we change the ID we do not know wich sub-feature parent attribute value to modify
		# (It could occur that we link sub feature not related)
	}
	return $is_dupli;
}

# find the level of the feature tested
sub get_level{
	my ($omniscient, $feature)=@_;

	my $source_tag = lc($feature->source_tag);
	my $primary_tag = lc($feature->primary_tag);

	my $level=undef;

	#########################################
	## PECULIARITIES FROM HAVANA / ENSEMBL ##
	if ($source_tag eq "ensembl" ){
		if ( $primary_tag eq "rna" ) {return 'level1';} #particularity ENSEMBL
	}
	if ( ($source_tag =~ "havana" or $source_tag =~ "ensembl") and ($primary_tag eq	"processed_transcript" or $primary_tag eq	"pseudogene" ) ){ #By default processed_transcript is l2 and pseudogene is l1
		if ($feature->has_tag('Parent')){return "level2" ;}
		else{return "level1" ;}
	}
	## PECULIARITIES FROM HAVANA / ENSEMBL ##
	#########################################

	if ( exists_keys($omniscient, ('other','level','level1', $primary_tag) ) ){
		return 'level1';
	}
	elsif( exists_keys($omniscient, ('other','level','level2', $primary_tag) ) ){
		return 'level2';
	}
	elsif( exists_keys($omniscient, ('other','level','level3', $primary_tag) ) ){
		return 'level3';
	}
}

# create string that should be uniq by feature
# will be used to detect duplicated features
sub _create_comparison_string{
	my ($omniscient, $uniqID, $feature)=@_;

	my $string=$feature->seq_id().$feature->primary_tag().$feature->start().$feature->end();
	my $primary_tag = lc($feature->primary_tag);

	# get the ID
	my $ID=undef;
	if(exists_keys($uniqID,( lc($feature->_tag_value('ID') ) ) ) ){
		$ID = $uniqID->{ lc($feature->_tag_value('ID') )};
	}
	else{
		$ID = $feature->_tag_value('ID')
	}

	# If we are checking a level1 feature no need to go further
	if ( exists_keys( $omniscient, ('other', 'level', 'level1', $primary_tag) ) ){
		$string .= $ID; # compare with original ID
		return $string;
	}

	# If we are not checking a level1 feature, let's take the parent info too
	my $Parent=undef;

	if(exists_keys($uniqID,( lc($feature->_tag_value('Parent') ) ) ) ) {
		$Parent = $uniqID->{ lc($feature->_tag_value('Parent') ) };
	}
	else{
		$Parent = $feature->_tag_value('Parent')
	}

	if ( exists_keys( $omniscient, ('other', 'level', 'level2', $primary_tag) ) ){
		$string .= $ID; # compare with original ID
		$string .= $Parent; # compare with original Parent
	}
	if ( exists_keys( $omniscient, ('other', 'level', 'level3', $primary_tag) ) ){
		$string .= $Parent; # compare with original Parent
	}

	return $string;
}

# create an ID uniq. Don't give multi-parent feature !
# If we have to create new ID for a SPREADFEATURES they will not have a shared ID.
sub _check_uniq_id{
	my	($omniscient, $miscCount, $uniqID, $uniqIDtoType, $feature)=@_;

	my $uID=undef;
	my $primary_tag = lc($feature->primary_tag);

	my $id=undef;
	if($feature->has_tag('ID')){ #has the tag
		$id = $feature->_tag_value('ID');
	}
	elsif($feature->has_tag($primary_tag."_id") ){
		$id = $feature->_tag_value($primary_tag."_id");
		create_or_replace_tag($feature, 'ID', $id);
	}
	elsif( get_level($omniscient, $feature) eq 'level1' and $feature->has_tag("gene_id") ){
		$id = $feature->_tag_value("gene_id");
		create_or_replace_tag($feature, 'ID', $id);
	}
	elsif( get_level($omniscient, $feature) eq 'level2' and $feature->has_tag("transcript_id") ){
		$id = $feature->_tag_value("transcript_id");
		create_or_replace_tag($feature, 'ID', $id);
	}

	# CHECK THE ID TO SEE IF IT's uniq, otherwise we have to create a new uniq ID
	if($id){
		# In case of non-spreadfeature (avoid CDS and UTR that can share identical IDs)
		if(! exists_keys($omniscient,('other', 'level', 'spreadfeature', $primary_tag) ) ){
			$uID = _create_ID($miscCount, $uniqID, $uniqIDtoType, $primary_tag, $id, PREFIX_NEW_ID); #method will push the uID
			if(	$id ne $uID ){ #push the new ID if there is one
				create_or_replace_tag($feature, 'ID', $uID);
			}
		}
		# In case of spreadfeature ( CDS and UTR that can share identical IDs)
		else{
			# First time we see this ID => No problem;
	 		if(! exists($uniqID->{ lc($id) })){
			 	#push the uID
			 	$uID = $id;
			 	$uniqID->{lc($uID)}=$id;
			 	$uniqIDtoType->{lc($id)}=$primary_tag;
			}
		  # NOT the first time we have this ID
			# check if it's the same type (To not mix a same ID between UTR and CDS);
			elsif( $uniqIDtoType->{lc($id)} eq $primary_tag ){ # Same type, so we can keep this ID, let's continue
			 	$uID = $id;
			}
			else{ # The spreadfeature type is different
				# Let's check if one of the same type is already in omniscient (THE ID could be linked to a non-spreadfeature), in that case we keep the ID already given.
				if( $feature->has_tag('Parent') ){
					if ( exists_keys( $omniscient, ('level3', $primary_tag, lc($feature->_tag_value('Parent')) ) ) ){
						$uID = 	@{ $omniscient->{'level3'}{$primary_tag}{lc($feature->_tag_value('Parent'))} }[0]->_tag_value('ID');
					}
				}
				if(! $uID){ #ID already taken by another feature type, and we do not have ID already existing of this feature type within omniscient, let's create a new ID
					$uID = _create_ID($miscCount, $uniqID, $uniqIDtoType, $primary_tag, $id, PREFIX_NEW_ID); #method will push the uID
				}
				if(	$id ne $uID ){ #push the new ID if there is one
				 		create_or_replace_tag($feature, 'ID', $uID);
				}
			}
		}
	}
	else{ #tag absent
		my $level = get_level($omniscient, $feature);
		if($level ne 'level3'){
			warn "gff3 reader error ".$level .": No ID attribute found".
			" @ for the feature: ".$feature->gff_string()."\n";
		}
		$miscCount->{$primary_tag}++;
		$id = $primary_tag."-".$miscCount->{$primary_tag}; # create an ID and then check if not already taken
		$uID = _create_ID($miscCount, $uniqID, $uniqIDtoType, $primary_tag, $id, PREFIX_NEW_ID); #method will push the uID
		create_or_replace_tag($feature, 'ID', $uID);
	}

	return $uID;
}

# create the ID and add it to the feature.
sub _create_ID{
	my	($miscCount, $uniqID, $uniqIDtoType, $primary_tag, $id, $prefix)=@_;

	my $key;

	if($prefix){
		$key=$prefix."-".$primary_tag;
	}
	else{
		$key=$primary_tag;
	}

	my $uID= $id ? $id : $key."-1";

	while( exists_keys($uniqID, (lc ($uID) ) )){	 #loop until we found an uniq tag
		$miscCount->{$key}++;
		$uID = $key."-".$miscCount->{$key};
	}

	#push the new ID
	$uniqID->{lc($uID)}=$id;
	$uniqIDtoType->{lc($uID)}=$primary_tag;

	return $uID;
}

# check if mRNA have a Parental feature existing. If not we create it.
sub _check_l1_linked_to_l2{
	my ($log, $hash_omniscient, $miscCount, $uniqID, $uniqIDtoType, $verbose, $debug)=@_;
	my $resume_case=undef;

	foreach my $primary_tag_l2 ( sort {$a cmp $b} keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
		foreach my $id_l1 ( sort {$a cmp $b} keys %{$hash_omniscient->{'level2'}{$primary_tag_l2}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
			my $l1_exist=undef;
			foreach my $primary_tag_l1 (keys %{$hash_omniscient->{'level1'}}){ # primary_tag_key_level1 = gene or repeat etc...
				if(exists_keys ($hash_omniscient, ('level1', $primary_tag_l1, $id_l1))){
					$l1_exist=1;
					last;
				}
			}
			# L1 is missing
			if(! $l1_exist){
				$resume_case++;
				my $gene_feature=clone($hash_omniscient->{'level2'}{$primary_tag_l2}{$id_l1}[0]);#create a copy of the first mRNA feature;

				#Deal case where we reconstruct other thing than a gene
				my $primary_tag_l1=undef;
				if(lc($gene_feature->primary_tag) =~ /match/){ $primary_tag_l1="match"; }
				else{ $primary_tag_l1="gene"; }

				$gene_feature->remove_tag('Parent'); # remove parent ID because, none.
				$gene_feature->primary_tag($primary_tag_l1); # change primary tag

				# ID and Parent ID identic, we have to update Parent ID
				my $ParentID = $hash_omniscient->{'level2'}{$primary_tag_l2}{$id_l1}[0]->_tag_value('Parent') ;
				my $l2_id = $hash_omniscient->{'level2'}{$primary_tag_l2}{$id_l1}[0]->_tag_value('ID');

				# Parent ID has same id as l2 feature !! We must modify it
				if ( lc( $l2_id ) eq  lc($ParentID) ){
					$ParentID = _create_ID($miscCount, $uniqID, $uniqIDtoType, $primary_tag_l1, $ParentID, PREFIX_ID_L1_NEW);
					dual_print($log, "Parent ID and ID are the same. Here is the new parent ID created $ParentID.\n") if ($debug);

					# Update new parent id to all feature l2 related
					foreach	my $l2_feature ( @{$hash_omniscient->{'level2'}{$primary_tag_l2}{$id_l1} } ){
						create_or_replace_tag($l2_feature,'Parent', $ParentID);
					}
					$hash_omniscient->{'level2'}{$primary_tag_l2}{lc($ParentID)} = delete $hash_omniscient->{'level2'}{$primary_tag_l2}{$id_l1};
				}

				# check start stop if isoforms exists
				check_level1_positions({ omniscient => $hash_omniscient, feature => $gene_feature } );

				# now save it in omniscient
				create_or_replace_tag($gene_feature, 'ID', $ParentID); #modify ID to replace by parent value
				$hash_omniscient->{"level1"}{$primary_tag_l1}{lc($ParentID)}=$gene_feature;
				dual_print($log, "No Parent feature found for ".$l2_id.". We create one: ".$gene_feature->gff_string()."\n", 0); # print only in log
			}
		}
	}

	if($resume_case){
		dual_print($log, "$resume_case cases fixed where L2 features have parent features missing\n", $verbose);
	}
	else{
		dual_print($log, "No problem found\n", $verbose);
	}
}

# @Purpose: Remove the level1 feature that havn't subfeature linked to it. Before to remove it check if L3 is linked to it. In that case it is a format error that we will fix !
# @input: 1 => hash(omniscient hash)
# @output: none
sub _remove_orphan_l1{
	my ($omniscient, $miscCount, $uniqID, $uniqIDtoType, $mRNAGeneLink, $verbose, $log, $debug)=@_;
	my $resume_case=undef;

 	foreach my $tag_l1 (keys %{$omniscient->{'level1'}}){
 			foreach my $id_l1 (keys %{$omniscient->{'level1'}{$tag_l1}} ){
				if ( exists_keys( $omniscient, ('other', 'level', 'level1', $tag_l1) ) ){
					if ($omniscient->{'other'}{'level'}{'level1'}{$tag_l1} eq 'standalone' or
								$omniscient->{'other'}{'level'}{'level1'}{$tag_l1} eq 'topfeature'){

						dual_print($log, "skip $tag_l1 because is suppose to be orphan\n") if $debug;
						next;
					}
				}

				my $neverfound="yes";
 				foreach my $tag_l2 (keys %{$omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
 						if ( exists_keys ( $omniscient,('level2',$tag_l2,$id_l1) ) ){
 							$neverfound=undef;last
 						}
 				}
 				if($neverfound){
 					$resume_case++;
 					print "removing ".$omniscient->{'level1'}{$tag_l1}{$id_l1}->gff_string."\n" if ($verbose >= 3);
					delete $omniscient->{'level1'}{$tag_l1}{$id_l1}; # delete level1 // In case of refseq the feature has been cloned and modified, it is why we nevertheless remove it
				}
 	 	}
 	}
	if($resume_case){
		dual_print($log, "$resume_case cases removed where L1 features do not have children (while they are suposed to have children).\n", $verbose);
	}
	else{
		dual_print($log, "None found\n", $verbose);
	}
}

# @Purpose: Check relationship betwwen L3 and L2. If L2 is missing we create it. When creating L2 missing we create as well L1 if missing too.
# @input: 4 => hash(omniscient hash), hash(mRNAGeneLink hash), hash(miscCount hash), hash(uniqID hash)
# @output: none
sub _check_l2_linked_to_l3{
	my ($log, $hash_omniscient, $mRNAGeneLink, $miscCount, $uniqID, $uniqIDtoType, $verbose, $debug)=@_;
	my $resume_case=undef;

 	foreach my $tag_l3 (sort {$a cmp $b} keys %{$hash_omniscient->{'level3'}}){

 			foreach my $id_l2 (sort {$a cmp $b} keys %{$hash_omniscient->{'level3'}{$tag_l3}}){

 				#check if L2 exits
 				if (! exists_keys($mRNAGeneLink, ( $id_l2 ) ) ) {
 					$resume_case++;

	 				#L3 linked directly to L1
	 				foreach my $tag_l1 (keys %{$hash_omniscient->{'level1'}}){ # primary_tag_key_level1 = gene or repeat etc...
						my $has_l1_feature = undef;
						my $id_l2_to_replace = undef;

						if(exists_keys ($hash_omniscient, ('level1', $tag_l1, $id_l2))){
							# case where it's linked by parent/ID attribute
							$has_l1_feature = $hash_omniscient->{"level1"}{$tag_l1}{$id_l2};
						}
						else{
							# Check if one as a common tag value == to L1 common tag value (then when creating l2 in check3 add parent for L2 of the L1 Id)
							foreach my $id_l1 (keys %{$hash_omniscient->{'level1'}{$tag_l1}}){
								my $l1_feature = $hash_omniscient->{"level1"}{$tag_l1}{$id_l1};
								foreach my $tag (@COMONTAG){
									#check if we have the tag
									if($l1_feature->has_tag($tag)){
										my $l1_ct_value=lc($l1_feature->_tag_value($tag)); #get the value
										foreach my $l3_feature (@{$hash_omniscient->{'level3'}{$tag_l3}{$id_l2}}){
											if($l3_feature->has_tag($tag) and	lc($l3_feature->_tag_value($tag)) eq	$l1_ct_value){
												$has_l1_feature = $l1_feature;
												$id_l2_to_replace = $l3_feature->_tag_value('Parent');
												# case where it's linked by comon_tag attribute
												last;
											}
										}
									}
									last if ($has_l1_feature);
								}
								last if ($has_l1_feature);
							}
				 		}

				 		if ($has_l1_feature){

							my $l1_ID = $has_l1_feature->_tag_value('ID');
							my $l2_feature = clone($has_l1_feature);#create a copy of the first mRNA feature;

							if (exists_keys($hash_omniscient,("level3",'cds', $id_l2) )	){
									$l2_feature->primary_tag('mRNA');
							}
							else{ #we cannot guess
								$l2_feature->primary_tag('RNA');
							}

	 						#Modify parent L2 (and L1 id if necessary)
	 						create_or_replace_tag($l2_feature,'Parent', $l1_ID); #modify ID to replace by parent value
	 						create_or_replace_tag($l2_feature,'ID', $id_l2_to_replace) if ($id_l2_to_replace); #modify ID to replace by parent value
							check_level2_positions($hash_omniscient, $l2_feature);

							if ( exists_keys ($uniqID,(lc($id_l2) ) ) ){ #the easiest is to modifiy the gene id

								my $new_l1id = _check_uniq_id($hash_omniscient, $miscCount, $uniqID, $uniqIDtoType, $has_l1_feature); # to check if ID was already in used by level1 feature
								create_or_replace_tag($l2_feature,'Parent', $new_l1id); #modify ID to replace by parent value
								my $primary_tag_l1 =$has_l1_feature->primary_tag();
								$hash_omniscient->{"level1"}{lc($primary_tag_l1)}{lc($new_l1id)} = delete $hash_omniscient->{"level1"}{lc($primary_tag_l1)}{lc($l1_ID)}; # now save it in omniscient
								#fill the $mRNAGeneLink hash
								$mRNAGeneLink->{ $id_l2 } = $new_l1id;
								push(@{$hash_omniscient->{"level2"}{lc($l2_feature->primary_tag)}{lc($new_l1id)}}, $l2_feature);
							}
	 						else{
	 							#fill the $mRNAGeneLink hash
	 							$mRNAGeneLink->{ $id_l2 } = $l1_ID; # Always need to keep track about l2->l1, else the method _check_l2_linked_to_l3 will recreate a l1 thinking this relationship is not fill
								push(@{$hash_omniscient->{"level2"}{lc($l2_feature->primary_tag)}{lc($l1_ID)}}, $l2_feature);
							}
							dual_print($log, "L3 was directly linked to L1. Corrected by creating the intermediate L2 feature from L1 feature:".$l2_feature->gff_string()."\n", 0);
							last
						}
					}

				if (! exists($mRNAGeneLink->{ $id_l2 }) ) { # it was not previous case (L3 linked directly to L1)

	 				#start fill L2
	 				my $l2_feature=clone($hash_omniscient->{'level3'}{$tag_l3}{$id_l2}[0]);#create a copy of the first mRNA feature;
					$l2_feature->frame(".") if ($l2_feature->frame ne "."); # If we clone a CDS there will be a frame information to remove.
					my $new_ID = $l2_feature->_tag_value('Parent');
					create_or_replace_tag($l2_feature,'ID', $new_ID); #modify ID to replace by parent value
					my $primary_tag_l2;
					if( exists_keys($hash_omniscient,('level3', 'cds', $id_l2)) ) {
						$primary_tag_l2="mRNA";
					}
					else{
						$primary_tag_l2="RNA";
					}
					$l2_feature->primary_tag($primary_tag_l2); # change primary tag
					check_level2_positions($hash_omniscient, $l2_feature);	# check start stop if isoforms exists

					#fill L1
					my $l1_feature=clone($hash_omniscient->{'level3'}{$tag_l3}{$id_l2}[0]);#create a copy of the first mRNA feature;
					$l1_feature->frame(".") if ($l1_feature->frame ne "."); # If we clone a CDS there will be a frame information to remove.
					$l1_feature->remove_tag('Parent'); # remove parent ID because, none.

					#Deal case where we reconstruct other thing than a gene
					my $primary_tag_l1=undef;
					if(lc($l1_feature->primary_tag) =~ /match/){ $primary_tag_l1="match"; }
					else{ $primary_tag_l1="gene"; }
					$l1_feature->primary_tag($primary_tag_l1); # change primary tag

					my $new_ID_l1 = _check_uniq_id($hash_omniscient, $miscCount, $uniqID, $uniqIDtoType, $l1_feature);
					create_or_replace_tag($l1_feature,'ID', $new_ID_l1); #modify ID to replace by parent value

				#finish fill Level2
					create_or_replace_tag($l2_feature, 'Parent', $new_ID_l1); # remove parent ID because, none.
					#save new feature L2
					push (@{$hash_omniscient->{"level2"}{lc($primary_tag_l2)}{lc($new_ID_l1)}}, $l2_feature);

				#finish fill Level1
					check_level1_positions({ omniscient => $hash_omniscient, feature => $l1_feature});	# check start stop if isoforms exists
					#save new feature L1
					$hash_omniscient->{"level1"}{lc($primary_tag_l1)}{lc($new_ID_l1)} = $l1_feature; # now save it in omniscient
					$mRNAGeneLink->{lc($id_l2)} = $new_ID_l1;

					dual_print($log, "L1 and L2 created:".$l1_feature->gff_string()."\n".$l2_feature->gff_string()."\n", 0);
	 				}
	 			}
 			}
 	}
	if($resume_case){
 		dual_print($log, "$resume_case cases fixed where L3 features have parent feature(s) missing\n", $verbose);
	}
	else{
		dual_print($log, "No problem found\n", $verbose);
	}
}

# @Purpose: Check L3 features. If exon are missing we create them. We go through all features of level3 and check them by type, if two should be merged, we do it (CDS 1-50 and 51-100, must be CDS 1-100).
# @input: 3 =>	hash(omniscient hash), hash(miscCount hash), hash(uniqID hash)
# @output: none
sub _check_exons{
	my ($debug, $log, $hash_omniscient, $mRNAGeneLink, $miscCount, $uniqID, $uniqIDtoType, $verbose)=@_;
	my $resume_case=undef; my $resume_case2=undef; my $resume_case3=undef; my $resume_case4=undef;

	my %checked;
	foreach my $tag_l3 ( sort {$a cmp $b} keys %{$hash_omniscient->{'level3'}}){
		if ($tag_l3 ne "exon"){
 				foreach my $id_l2 ( sort { (($a =~ /(\d+)$/)[0] || 0) <=> (($b =~ /(\d+)$/)[0] || 0) } keys %{$hash_omniscient->{'level3'}{$tag_l3}}){

 					if( ! exists_keys(\%checked,($id_l2)) ){ #l2 already checked
 						dual_print($log, "Check: $id_l2\n") if ($debug);
 						my $feature_example=undef; # will be used to create the exon features
	 					my $list_location_Exon=undef;
	 					my $list_location_NoExon=undef;
	 					my $list_location_NoExon_overlap=undef;
	 					my %createIT; # will be usefull to list the feature to create
#				 	+-----------------------------------------------------
#					| 			Go through l3 and save info needed		 |
#				 	+-----------------------------------------------------

	 					foreach my $tag_l3 ( sort {$a cmp $b} keys %{$hash_omniscient->{'level3'}}){

	 						# LIST NON-EXON LOCATIONS THAT NEED TO BE IN AN EXON LOCATION
	 						if ($tag_l3 ne "exon" and $hash_omniscient->{'other'}{'level'}{'level3'}{$tag_l3} eq "exon" ){

				 				if( exists_keys($hash_omniscient,('level3',$tag_l3, $id_l2)) ){

				 					my $list_location_l3=[];
				 					foreach my $l3_feature (@{$hash_omniscient->{'level3'}{$tag_l3}{$id_l2}}){

				 						if(! $feature_example){
				 							$feature_example=$l3_feature;
				 						}

				 						my $locationRefList=[[[$l3_feature->_tag_value('ID')] ,int($l3_feature->start), int($l3_feature->end)]];
				 						$list_location_l3 = _manage_location($locationRefList, $list_location_l3, 'adjacent', 0); # we check first in overlap mode to check if badly define features exists
				 					}

				 					#Rare case when a feature of the same type is badly defined
				 					if(@{$list_location_l3} < @{$hash_omniscient->{'level3'}{$tag_l3}{$id_l2}}){
				 						my $message = "Peculiar rare case, we found ".@{$list_location_l3}." ".$tag_l3." while ".@{$hash_omniscient->{'level3'}{$tag_l3}{$id_l2}}." expected.\n";
										$message .=	"Either some are supernumerary or some have been merged because they overlap or are adjacent while they are not suppose to.\n";
				 						$message .= "In case you were using gtf file as input (no parent/id attributes), check you provide the attribute (i.e comon_tag) used to group features together (e.g. locus_tag, gene_id, etc.).\n";
				 						$message .= "(In case your file contains only CDS features, and your organism is prokaryote (e.g rast file), using ID as comon_tag might be the solution.)\n";
				 						warn $message;
				 					}
				 					push @{$list_location_NoExon_overlap}, @{$list_location_l3}; #list of all feature that has been checked in overlap mode
				 				}
				 			}

				 			# LIST EXON LOCATIONS
				 			elsif($tag_l3 eq "exon"){

							if( exists_keys($hash_omniscient,('level3',$tag_l3, $id_l2)) ){

				 					foreach my $l3_feature (@{$hash_omniscient->{'level3'}{$tag_l3}{$id_l2}}){

				 						if(! $feature_example){
				 							$feature_example=$l3_feature;
				 						}
				 						my $locationRefList=[[[$l3_feature->_tag_value('ID')] ,int($l3_feature->start), int($l3_feature->end)]];
				 						$list_location_Exon = _manage_location($locationRefList, $list_location_Exon , 'adjacent', 0);
				 					}

				 					#Rare case when a features are badly defined
				 					# This approch works for exon because they have uniq ID
				 					if(@{$list_location_Exon} < @{$hash_omniscient->{'level3'}{$tag_l3}{$id_l2}}){
				 						dual_print($log, "Peculiar rare case, we have to remove existing exon which are supernumerary. Parent is $id_l2\n", 0); #print only in log
				 						#remove the non needed features (they where wrong, and are unecessary)
									my @id_list2=();
									foreach my $locations (@{$list_location_Exon}){
										# If several value in the ID list, we have to avoid the first value (which is the one to keep), and remove all the other ones which are the Id to remove
										if(@{$locations->[0]} > 1){
											my $correct_ID = shift @{$locations->[0]};
											push @id_list2, @{$locations->[0]};
											@{$locations->[0]} = $correct_ID;
											foreach my $l3_feature (@{$hash_omniscient->{'level3'}{'exon'}{$id_l2} } ){
			 										if($l3_feature->_tag_value('ID') eq	$correct_ID){
			 											$l3_feature->start($locations->[1]);
			 											$l3_feature->end($locations->[2]);
			 										}
			 									}
										}
									}
				 					my @tag_list = ('all');
				 					my @id_list=($id_l2);
				 					$resume_case3 += @id_list2;
				 					dual_print($log, "We remove the supernumerary @id_list2 exon(s)\n", 0); # print only in log
									remove_element_from_omniscient(\@id_list, \@id_list2, $hash_omniscient, 'level3', 'false', \@tag_list);
				 				}
				 			}
				 		}
				 	}

				 	## check all NOn-exon in adjacent mater to have the correct list (allows to merge UTR and CDS to recreate exon )
				 	foreach my $location (@{$list_location_NoExon_overlap}){
				 		$list_location_NoExon = _manage_location([$location], $list_location_NoExon, 'adjacent', 0);
				 	}

				 	#print "list_location_Exon: ".Dumper($list_location_Exon) if ($verbose >= 3);
				 	#print "list_location_NOEXON: ".Dumper($list_location_NoExon) if ($verbose >= 3);

#				 	+--------------------------------------------------------------------------------------------------------+
#					| 				COMPARE EXONS POSITION TO THOSE DESCRIBED BY NON-EXON FEATURES 						 |
#				 	+--------------------------------------------------------------------------------------------------------+

 						#Case where exon feature exists, we have to check them
	 					if( exists_keys($hash_omniscient,('level3','exon', $id_l2)) ){ #When thre are l3 features but no exon among them... we need to recreate them.

	 						if ($list_location_NoExon){ #We have features like UTR,CDS,etc allowing to check the exon locations.
		 						#create string to comapre the 2 lists.
		 						my $list_location_Exon_joined="";
		 						foreach my $location (sort {$a->[1] <=> $b->[1] } @{$list_location_Exon}){
		 							$list_location_Exon_joined .= $location->[1].$location->[2];
		 						}
		 						my $list_location_NoExon_joined="";
		 						foreach my $location (sort {$a->[1] <=> $b->[1] } @{$list_location_NoExon}){
		 							$list_location_NoExon_joined .= $location->[1].$location->[2];
		 						}
		 						#If two lists are different we have to check/fix the difference
		 						# If no overlap we create the exon:
		 						# If overlap:	Redefine internal exon ; Redfine external exon only if too short.
							if($list_location_Exon_joined ne $list_location_NoExon_joined ){
								dual_print($log, "Problem for $id_l2! coordinates of exons found not in\n".
								"agreement with those expected according to the other features\n".
								"(i.e. CDS and/or UTR, etc). Let's check that (We will create exon, or modify\n".
								"coordinates, depending of cases. If creation of UTR is needed it will be done in\n".
								"a next step) !!\n") if($debug);

								my $location_cpt=0;
			 						foreach my $location (sort {$a->[1] <=> $b->[1] } @{$list_location_NoExon}){
			 							$location_cpt++;

			 							my $create_exon=1;
			 							my $new_location;
			 							my $overlap;

			 							foreach my $exon_location (sort {$a->[1] <=> $b->[1] } @{$list_location_Exon}){

			 								($new_location, $overlap) = _manage_location_lowLevel_adjacent($location, $exon_location); #there is an overlap if $new_location != $exon_location. If it's the same, we should check $overlap to be sure

			 								if($new_location->[1] < $exon_location->[1] or $new_location->[2] > $exon_location->[2] ){ #The exon_location has been modified by location... We have to remodelate the exon (only if fit some conditions) location to take the modification into account
				 								$create_exon=undef; # We must avoid to create exon because there is an overlap.

			 									my $redefine_left=undef;
			 									my $redefine_right=undef;
			 									#first location => check left
				 								if($location_cpt == 1){
				 									if($new_location->[1] <	$exon_location->[1]){ $redefine_left = $new_location->[1];} # Modify only if it's more left
				 								}
				 								#=> check left and right
				 								if($location_cpt != 1 and $location_cpt != @$location){
				 									if($new_location->[1] <	$exon_location->[1]){ $redefine_left = $new_location->[1];}	# Modify only if it's more left
				 									if($new_location->[2] >	$exon_location->[2]){ $redefine_right = $new_location->[2];} # Modify only if it's more right
				 								}
				 								#last location => check right
				 								if($location_cpt == @$location){
				 									if($new_location->[2] >	$exon_location->[2]){ $redefine_right = $new_location->[2];} # Modify only if it's more right
				 								}

				 								foreach my $l3_feature (@{$hash_omniscient->{'level3'}{'exon'}{$id_l2} } ){
				 									if($l3_feature->_tag_value('ID') eq $exon_location->[0][0]){

				 										if($redefine_left){
															dual_print($log, "Modify the left location for ".$l3_feature->_tag_value('ID')." from ".$l3_feature->start()." to ".$new_location->[1]."\n", 0); #print only in log
				 											$l3_feature->start($new_location->[1]); $resume_case2++;
				 										}
														else{$redefine_left = $exon_location->[1];}

				 										if($redefine_right){
															dual_print($log, "Modify the right location for ".$l3_feature->_tag_value('ID')." from ".$l3_feature->end()." to ".$new_location->[2]."\n", 0); #print only in log
				 											$l3_feature->end($new_location->[2]); $resume_case2++;
				 										}
														else{$redefine_right = $exon_location->[2];}

				 										last;
				 									}
				 								}
				 							}
				 							elsif($overlap){ #location not modified but no overlap, so it means the exon is not defined ! <= ?? Not sure this comment is good 27th Nov 2018
				 								$create_exon=undef;
				 							}
				 						}

				 						if($create_exon){
				 							push @{$createIT{'exon'}}, $location;
				 						}
		 						}
		 					}
		 				}
		 				else{ dual_print($log, "No other feature to check the exon locations (e.g CDS, UTR, etc...). We can trust them then.\n") if ($debug)}
	 					}
	 					else{ $createIT{'exon'}=$list_location_NoExon;} # no exon exists, we have to create all of them

					# NOW CREATE EXON IF NECESSARY
					if(keys %createIT){
						foreach my $tag (keys %createIT){
					 		foreach my $location (@{$createIT{$tag}}){
					 			$resume_case++;
					 			my $feature_exon = clone($feature_example);#create a copy of a random feature l3;
								$feature_exon->start($location->[1]);
					 			$feature_exon->end($location->[2]);
					 			$feature_exon->frame(".");
					 			$feature_exon->primary_tag($tag);
					 			my $uID = _check_uniq_id($hash_omniscient, $miscCount, $uniqID, $uniqIDtoType, $feature_exon);
					 			create_or_replace_tag($feature_exon, 'ID', $uID); # remove parent ID because, none.
								#save new feature L2
								dual_print($log, "Create one Exon for $id_l2\n:".$feature_exon->gff_string."\n", 0); #print only in log
								push (@{$hash_omniscient->{"level3"}{$tag}{$id_l2}}, $feature_exon);
					 			}
					 	}
				 	}

				 	#Check extremities of exons (If exon is shorter we adapt it to the mRNA size, else we adapt the L2 size to the exon size)
	 					my $id_l1 = lc($mRNAGeneLink->{lc($id_l2)});
	 					my $getout=undef;
	 					foreach my $tag_l2 ( %{$hash_omniscient->{'level2'}} ){
	 						if( exists_keys($hash_omniscient,('level2', $tag_l2, $id_l1)) ){
	 							foreach my $l2_feature ( @{$hash_omniscient->{'level2'}{$tag_l2}{$id_l1}} ){
	 								if( lc($l2_feature->_tag_value('ID')) eq $id_l2 ){
	 								 	if( exists_keys ($hash_omniscient, ('level3', 'exon', $id_l2))) { # If no exon it could be a case whre no L3 feature need an exon like non_canonical_three_prime_splice_site (they are out of exon). So the list of exon does not exist.
		 								 	my $myLeftExtremity=$l2_feature->start();
		 								 	my $myRightExtremity=$l2_feature->end();

					 					 	my @list_exon = sort {$a->start <=> $b->start} @{$hash_omniscient->{'level3'}{'exon'}{$id_l2}};

					 					 	if( int($list_exon[0]->start) >	int($myLeftExtremity) ){
												dual_print($log, "Modify the exon ".$list_exon[0]->_tag_value('ID')." LEFT extremity for $id_l2! From ".$list_exon[0]->start." to ".$myLeftExtremity."\n", 0); #print only in log
					 					 		$list_exon[0]->start($myLeftExtremity);
												$resume_case2++;
					 					 	}
					 					 	if($list_exon[0]->start <	$myLeftExtremity){	#modify L2
												dual_print($log, "Modify the L2 $id_l2 LEFT extremity! From ".$l2_feature->start."to".$list_exon[0]->start."\n", 0); #print only in log
					 					 		$l2_feature->start($list_exon[0]->start);
												$resume_case4++;
					 					 	}

					 					 	if($list_exon[$#list_exon]->end <	$myRightExtremity){
					 					 		dual_print($log, "Modify the exon ".$list_exon[$#list_exon]->_tag_value('ID')." RIGHT extremity for $id_l2! From ".$list_exon[$#list_exon]->end." to ".$myRightExtremity."\n", 0); #print only in log
					 					 		$list_exon[$#list_exon]->end($myRightExtremity);
												$resume_case2++;
					 						}
					 						elsif($list_exon[$#list_exon]->end >	$myRightExtremity){ #modify L2
												dual_print($log, "Modify the L2 $id_l2 RIGHT extremity! From ".$l2_feature->end."to".$list_exon[$#list_exon]->end."\n", 0); #print only in log
					 							$l2_feature->end($list_exon[$#list_exon]->end);
												$resume_case4++;
					 						}

					 					 	$getout=1;
					 					 	last;
					 					}
					 				}
	 							}
	 							if($getout){
	 								last;
	 							}
	 						}
	 					}

				 	#keep track of l2 checked (as we loop over L3, we meet several time the same l2)
		 				$checked{$id_l2}++;
	 				}
 				}
 			}
 	}
	if( $resume_case ){
		dual_print($log, "$resume_case exons created that were missing\n", $verbose);
	}
	else{ dual_print($log, "No exons created\n", $verbose); }
	if( $resume_case2 ){
		dual_print($log, "$resume_case2 exons locations modified that were wrong\n", $verbose);
	}
	else{ dual_print($log, "No exons locations modified\n", $verbose); }
	if( $resume_case3 ){
		dual_print($log, "$resume_case3 exons removed that were supernumerary\n", $verbose);
	}
	else{ dual_print($log, "No supernumerary exons removed\n", $verbose);}
	if( $resume_case4 ){
		dual_print($log, "$resume_case4 level2 locations modified\n", $verbose);
	}
	else{ dual_print($log, "No level2 locations modified\n", $verbose);}
}

# @Purpose: Check L3 features. If UTRS are missing we create them.
# @input: 3 =>	hash(omniscient hash), hash(miscCount hash), hash(uniqID hash)
# @output: none
sub _check_utrs{
	my ($debug, $log, $hash_omniscient, $mRNAGeneLink, $miscCount, $uniqID, $uniqIDtoType, $verbose)=@_;
	my $resume_case=undef;my $resume_case2=undef; my $resume_case3=undef;

	my %checked;
	foreach my $tag_l3 (sort {$a cmp $b} keys %{$hash_omniscient->{'level3'}}){
		if ($tag_l3 ne "exon"){
 				foreach my $id_l2 (sort { (($a =~ /(\d+)$/)[0] || 0) <=> (($b =~ /(\d+)$/)[0] || 0) }	keys %{$hash_omniscient->{'level3'}{$tag_l3}}){

 					if( ! exists_keys(\%checked,($id_l2)) ){ #l2 already checked

 						my $feature_example=undef; # will be used to create the exon features
	 					my $list_location_Exon=undef;
	 					my $list_location_CDS=undef;
	 					my $list_location_UTR=undef;

#				 		+-----------------------------------------------------
#						| 			Go through l3 and save info needed					 |
#				 		+-----------------------------------------------------
	 					foreach my $tag_l3 (keys %{$hash_omniscient->{'level3'}}){

				 			# LIST CDS LOCATIONS
	 						if ($tag_l3 eq "cds"){
				 				if( exists_keys($hash_omniscient,('level3',$tag_l3, $id_l2)) ){

				 					foreach my $l3_feature (@{$hash_omniscient->{'level3'}{$tag_l3}{$id_l2}}){

				 						my $locationRefList=[[[$l3_feature->_tag_value('ID')] ,int($l3_feature->start), int($l3_feature->end)]];
				 						$list_location_CDS = _manage_location($locationRefList, $list_location_CDS, 'adjacent', $verbose);
				 					}
				 				}
				 			}

				 			# LIST UTR LOCATIONS
	 						if ($tag_l3 =~ "utr"){
				 				if( exists_keys($hash_omniscient,('level3',$tag_l3, $id_l2)) ){

				 					foreach my $l3_feature (@{$hash_omniscient->{'level3'}{$tag_l3}{$id_l2}}){

				 						my $locationRefList=[[[$l3_feature->_tag_value('ID')] ,int($l3_feature->start), int($l3_feature->end)]];
				 						$list_location_UTR = _manage_location($locationRefList, $list_location_UTR, 'adjacent', $verbose);
				 					}
				 				}
				 			}

				 			# LIST EXON LOCATIONS
				 			elsif($tag_l3 eq "exon"){
								if( exists_keys($hash_omniscient,('level3',$tag_l3, $id_l2)) ){

				 					foreach my $l3_feature (@{$hash_omniscient->{'level3'}{$tag_l3}{$id_l2}}){

				 						if(! $feature_example){
				 							$feature_example=$l3_feature;
				 						}
				 						#print "exonFeature= ".$l3_feature->gff_string."\n";
				 						push @{$list_location_Exon}, [ [$l3_feature->_tag_value('ID')], int($l3_feature->start), int($l3_feature->end)] ;
				 					}
				 				}
				 			}
				 	}

#				 	+-----------------------------------------------------
#					| 				HANDLE UTRs 						 										|
#				 	+-----------------------------------------------------
					if( exists_keys($hash_omniscient,('level3','cds', $id_l2)) ){ #Check UTR only if CDS exists

						# Create list of UTR expected:
						my $list_location_UTR_expected=undef;
						my $expected_utr=1;

						foreach my $exon_location (sort {$a->[1] <=> $b->[1] } @{$list_location_Exon}){

			 					my $new_location;
			 					my $overlap;
			 					my $never_overlap="yes";
							foreach my $location_cds (sort {$a->[1] <=> $b->[1] } @{$list_location_CDS}){

								if( $location_cds->[1] > $exon_location->[2]){last;}
								if( $location_cds->[2] < $exon_location->[1]){next;}

								($new_location, $overlap) = _manage_location_lowLevel_inversed($location_cds, $exon_location, $verbose);

								if($overlap eq "perfect"){ $never_overlap=undef; $expected_utr=undef;last;}

								if($new_location->[1] != $exon_location->[1] and $new_location->[2] != $exon_location->[2] ){ #two UTR expected				=========================	exon
									dual_print($log, "creation utr push1\n") if ($debug);
									push @{$list_location_UTR_expected}, [undef, $exon_location->[1], $location_cds->[1]-1];				#								=======			CDS
									push @{$list_location_UTR_expected}, [undef, $location_cds->[2]+1, $exon_location->[2]];
									$never_overlap=undef;
									last;
								}
								elsif($new_location->[1] != $exon_location->[1] or $new_location->[2] != $exon_location->[2] ){ #two UTR expected	{
									#print "creation utr push2\n".Dumper($new_location)."\n" if($verbose >= 3);
									push @{$list_location_UTR_expected}, $new_location;
									$never_overlap=undef;
									last;
								}
							}
							if($never_overlap){ #case of UTR that match++ fully the exon
								push @{$list_location_UTR_expected}, $exon_location;
							}
						}

						#print "list_location_UTR_expected: ".Dumper($list_location_UTR_expected) if ($verbose >= 3);
						#print "list_location_UTR: ".Dumper($list_location_UTR) if ($verbose >= 3);

		 					# Compare UTR Present and UTR expected
	 						my $list_utr_to_create=undef;

	 						if($list_location_UTR){ #List UTR not empty
								if($list_location_UTR_expected){ #List UTR not empty
				 					foreach my $UTRexp_location (sort {$a->[1] <=> $b->[1] } @{$list_location_UTR_expected} ){

			 							my $create_utr=1;
			 							my $new_location;
			 							my $overlap;
			 							foreach my $UTR_location (sort {$a->[1] <=> $b->[1] } @{$list_location_UTR}){

			 								($new_location, $overlap) = _manage_location_lowLevel_inversed($UTR_location, $UTRexp_location, $verbose); #just to check that it overlaps

			 								if($overlap and ( $UTR_location->[1] != $UTRexp_location->[1] or $UTR_location->[2] != $UTRexp_location->[2] ) ){ #It overlaps and at least one location is different. We have to re-modelate the utr location to take the modification into account
				 								$create_utr=undef;

				 								foreach my $tag_l3 (keys %{$hash_omniscient->{'level3'} } ){
				 									if($tag_l3 =~"utr"){
				 										if( exists_keys($hash_omniscient,('level3', $tag_l3, $id_l2)) ){
							 								foreach my $l3_feature (@{$hash_omniscient->{'level3'}{$tag_l3}{$id_l2} } ){
							 									if($l3_feature->_tag_value('ID') eq $UTR_location->[0][0] and $l3_feature->start eq $UTR_location->[1] and $l3_feature->end eq $UTR_location->[2]){ # we have to check position to be sure we modify the correct one, because UTR could share the same ID
																	$l3_feature->start($UTRexp_location->[1]);
							 										$l3_feature->end($UTRexp_location->[2]);
																	$resume_case2++;
																	dual_print($log, "Modify the UTR: ".$UTR_location->[0][0]." location for $id_l2! From ".$UTR_location->[1]." ".$UTR_location->[2]." to ".$UTRexp_location->[1]." ".$UTRexp_location->[2]."\n", 0); # print log only
							 										last;
							 									}
							 								}
							 							}
							 						}
							 					}
				 							}
				 							elsif($overlap and $overlap eq "perfect"){ #An UTR that match perfectly already exists !
				 								$create_utr=undef;
				 							}
				 						}

				 						if($create_utr){
				 							push @{$list_utr_to_create}, $new_location;
				 						}
			 					}
			 				}
			 				else{
								dual_print($log, "Error in the file, we have an UTR but none".
								" is expected according to the described exons. @ Level2 studied: $id_l2 \n") if($debug); # print in log only

								my $original_nb_utr=0;
								#let's remove UTR we will re-create new ones
								foreach my $tag_l3 (sort {$a cmp $b} keys %{$hash_omniscient->{'level3'}}){
									if ($tag_l3 =~ "utr"){
						 				if( exists_keys($hash_omniscient,('level3',$tag_l3, $id_l2)) ){
											$original_nb_utr = @{$hash_omniscient->{'level3'}{$tag_l3}{$id_l2}};
						 					delete $hash_omniscient->{'level3'}{$tag_l3}{$id_l2};
						 				}
									}
								}
								# list UTR to create if any has to be created
								my $nb_supernumary_utrs = $original_nb_utr;
								if($list_location_UTR_expected){
									$nb_supernumary_utrs = $original_nb_utr - @{$list_location_UTR_expected};
									$list_utr_to_create=$list_location_UTR_expected;# no UTR exists, we have to create all of them
								}
								$resume_case3 += $nb_supernumary_utrs;
								dual_print($log, "Remove $nb_supernumary_utrs supernumerary UTRs for $id_l2\n", 0)# print in log only
							}
		 				}
	 					else{
	 						if($list_location_UTR_expected){
	 							$list_utr_to_create=$list_location_UTR_expected;# no UTR exists, we have to create all of them
	 						}
 						}
 						#print "list_utr_to_create: ".Dumper($list_utr_to_create) if ($verbose >= 3);

						# NOW CREATE UTR IF NECESSARY
						my @cds_sorted = sort {$a->[1] <=> $b->[1]} @{$list_location_CDS};

						my $extremLeftCDS = $cds_sorted[0]->[1];
						my $extremRightCDS = $cds_sorted[$#cds_sorted]->[2];

						if($list_utr_to_create){

					 			foreach my $location (@{$list_utr_to_create}){
					 				$resume_case++;

									my $feature_utr = clone($feature_example);#create a copy of a random feature l3;
									$feature_utr->start($location->[1]);
									$feature_utr->end($location->[2]);
									$feature_utr->frame(".");

					 				#HANDLE primary tag
					 				my $primary_tag = "UTR";
					 				if($location->[2] < $extremLeftCDS){
					 					if($feature_utr->strand == 1){
					 						$primary_tag = "five_prime_UTR";
					 					}
					 					else{
					 						$primary_tag = "three_prime_UTR";
					 					}
					 				}
					 				elsif($location->[1] > $extremRightCDS){
					 					if($feature_utr->strand == 1){
					 						$primary_tag = "three_prime_UTR";
					 					}
					 					else{
					 						$primary_tag = "five_prime_UTR";
					 					}
					 				}

									$feature_utr->primary_tag($primary_tag);

									my $uID = _check_uniq_id($hash_omniscient, $miscCount, $uniqID, $uniqIDtoType, $feature_utr);
									create_or_replace_tag($feature_utr, 'ID', $uID); # remove parent ID because, none.
									#save new feature L2
									dual_print($log, "Create one UTR for $id_l2\n:".$feature_utr->gff_string."\n", 0); #print only in log
									push (@{$hash_omniscient->{"level3"}{lc($primary_tag)}{$id_l2}}, $feature_utr);
								}
						}
					}

					#keep track of l2 checked (as we loop over L3, we meet several time the same l2)
					$checked{$id_l2}++;
				}
			}
		}
 	}
	if( $resume_case ){
		dual_print($log, "$resume_case UTRs created that were missing\n", $verbose);
	}
	else{ dual_print($log, "No UTRs created\n", $verbose); }
	if( $resume_case2 ){
		dual_print($log, "$resume_case2 UTRs locations modified that were wrong\n", $verbose);
	}
	else{ dual_print($log, "No UTRs locations modified\n", $verbose); }
	if( $resume_case3 ){
		dual_print($log, "$resume_case3 UTRs removed that were supernumerary\n", $verbose);
	}
	else{ dual_print($log, "No supernumerary UTRs removed\n", $verbose);}
}

# @Purpose: Will merge a list of "location" (tuple of integer), and another list of location.
#					 If two location overlap or are adjacent, only one location will be kept that represent the most extrem values
# @input: 3 =>	list of 3 values([[S,X,Y][S,Z,W]] or [[[S],X,Y]]),	list of integer tuple, verbose option for debug
# @output: list of list
sub _manage_location{
	my ($locationRefList, $locationTargetList, $method, $verbose) = @_;

	my @new_location_list; #new location list that will be returned once filled

	if ($locationTargetList and @$locationTargetList >= 1){ #check number of location -> List not empty

		my @locations = (@{$locationRefList},@{$locationTargetList});
		my @locations_sorted = sort {$a->[1] <=> $b->[1]} @locations;

		my $location_modified = undef;
		my $overlap = undef;
		my $location1=undef;
		my $location2 = undef;

		foreach my $i (0 .. $#locations_sorted-1){

			if($location_modified){
				$location1 = $location_modified;
				$location_modified = undef;
			}
			else{
				$location1 = $locations_sorted[$i];
			}
			$location2 = $locations_sorted[$i+1];

			if ($location2->[1] > $location1->[2]+1 and $i+1 != $#locations_sorted){ #locations do not overlap and not before last of the list
				push @new_location_list, [@$location1];
			}
			else{
				my $location_back = undef;
				if($method eq 'adjacent'){
					($location_back, $overlap) = _manage_location_lowLevel_adjacent($location1, $location2);
					if($overlap){$location_modified=$location_back;}
				}
				else{
					($location_back, $overlap) = _manage_location_lowLevel_overlap($location1, $location2);
					if($overlap){$location_modified=$location_back;}
				}
			}
		}
		#if last round was not overlaping we should add the last value in the list
		if( $overlap ){
			push @new_location_list, [@$location_modified];
		}
		else{
			push @new_location_list, [@$location1];
			push @new_location_list, [@$location2];
		}
	}
	else{#check number of location -> none
		return \@{$locationRefList};
	}
	return \@new_location_list;
}

#	===================== location1
#		===================== location2
#	 ========================= <= New location2 returned
# @Purpose: Modify the location2 if it overlap the location1 by keeping the extrem values. Return the location2 intact if no overlap. /!\ The locations are merged if they are contigu
# @input: 2 =>	integer tuple [[ID],X,Y],	list of integer tuple
# @output: 2 => ref of a list of 2 element, boolean
sub _manage_location_lowLevel_adjacent{
	my ($location, $location2) = @_;

	my $new_location = [@{$location2}];
	my $overlap=undef;

	if ( ($location2->[1] <= $location->[2]+1) and ($location2->[2]+1 >= $location->[1]) ){ #it overlaps or are consecutive

		#Manage Id to avoid same IDs
		my %params = map { $_ => 1 } @{$new_location->[0]};
		foreach my $id ( @{$location->[0]}){
			if(! exists($params{$id})){
				push @{$new_location->[0]}, $id ; #append from the end the list of ID
			}
		}
		$overlap=1;

		if($location2->[1] > $location->[1]){
			$new_location->[1]=$location->[1];
		}
		if($location->[2] > $location2->[2]){
			$new_location->[2]=$location->[2];
		}
	}
	return $new_location, $overlap;
}

# ===================== location1
#		 ===================== location2
# ========================= <= New location2 returned
# @Purpose: Modify the location2 if it overlap the location1 by keeping the extrem values. Return the location2 intact if no overlap. /!\ We append the ID list by the end (as push) when there is an overlap
# @input: 2 =>	integer tuple [[ID],X,Y],	list of integer tuple
# @output: 2 => ref of a list of 2 element, boolean
sub _manage_location_lowLevel_overlap{
	my ($location, $location2) = @_;

	my $new_location = [@{$location2}];
	my $overlap=undef;

	if ( ($location2->[1] <= $location->[2]) and ($location2->[2] >= $location->[1]) ){ #it overlaps or are consecutive

		#Manage Id to avoid same IDs
		my %params = map { $_ => 1 } @{$new_location->[0]};
		foreach my $id ( @{$location->[0]}){
			if(! exists($params{$id})){
				push	@{$new_location->[0]}, $id ; #append from the end the list of ID
			}
		}

		$overlap=1;

		if($location2->[1] > $location->[1]){
			$new_location->[1]=$location->[1];
		}
		if($location->[2] > $location2->[2]){
			$new_location->[2]=$location->[2];
		}
	}

	return $new_location, $overlap;
}


#	================= 			location1 (cds)
#		===================== location2 (exon)
#										======== <= New location2 returned
sub _manage_location_lowLevel_inversed{
	my ($location, $location2, $verbose) = @_;

	my $new_location = [@{$location2}];
	my $overlap=undef;

	if ( ($location2->[1] == $location->[1]) and ($location2->[2] == $location->[2]) ){ #it overlaps perfectly
		return $new_location, "perfect";
	}

	if ( ($location2->[1] <= $location->[2]) and ($location2->[2] >= $location->[1]) ){ #it overlaps

		$overlap=1;

		if($location2->[1] < $location->[1]){
			$new_location->[2] = $location->[1]-1;
		}
		if($location->[2] < $location2->[2]){
			$new_location->[1] = $location->[2]+1;
		}
	}
	return $new_location, $overlap;
}

#===============================================================================
#Explanation: Case where part of the locus BBBBBB has been seen before to meet
#   its Parent feature (see below) = a parent feature ID has been created on
#   the fly during the parsing.
#   We now need to remove the wrong Parent ID and link them to the correct one.
#seq1	maker	CDS	561401	561519	.	+	2	ID=CLUHART00000006146:cds;locus_tag=AAAAA
#seq1	maker	UTR	337818	337914	.	+	.	ID=CLUHART00000008717:five_prime_utr;locus_tag=BBBBBB
#seq1	maker	UTR	343034	343277	.	+	.	ID=CLUHART00000008717:three_prime_utr;locus_tag=BBBBBB
#seq1	maker	CDS	564171	564235	.	+	0	ID=CLUHART00000006146:cds;locus_tag=AAAAA
#...
#seq1	maker	gene	337818	343277	.	+	.	ID=CLUHARG00000005458;locus_tag=BBBBBB
#
# HIS<=>Hash InfoSequential
sub _deinterleave_sequential{
	my ($infoSequential, $locusTAGuniq, $verbose, $log) = @_;
	my $resume_case=undef;

	foreach my $locusNameHIS (keys %{$infoSequential} ){

	 	if(exists_keys($locusTAGuniq,('level1', $locusNameHIS))){

	 		my $locusNameUniq = $locusTAGuniq->{'level1'}{$locusNameHIS};
	 		if($locusNameHIS ne $locusNameUniq ){
				dual_print($log, "Locus $locusNameHIS interleaved\n", 0); # print only in log
	 			$resume_case++;

	 			# The locusNameUniq already exists, we have to fill it with the part of
				# information missing that is contained in$infoSequential->{$locusNameHIS}
	 			if(exists_keys ($infoSequential,($locusNameUniq) ) ){
	 				foreach my $bucket (keys %{$infoSequential->{$locusNameHIS}} ){
	 					if ($bucket eq 'level1'){next;}

	 					my $prefix= lc(PREFIX_ID_L2_NEW); #when a l2 start with this prefix it means we created the l2 on the fly (the real l2 if exists, had not been met yet)
	 					if($bucket =~ /^$prefix/i){
	 						my $idok=undef;
	 						foreach my $feature ( @{$locusTAGuniq->{'level2'}{ $locusNameUniq}}){
	 							if(lc($feature->_tag_value('ID')) !~ /^$prefix/i){
	 								$idok = $feature->_tag_value('ID'); # @{$locusTAGuniq->{'level2'}{ $locusNameUniq }}[$cpt] is the first l2 feature that has been realy met
	 								last;
									# We make the assumption that the pieces of the locus that were lost before to describe its real l2 is part of the first real l2 met.
	 								# ====================================================================================================================================
	 							}
	 						}

	 				 		if(exists_keys ($infoSequential,($locusNameUniq, $idok) ) ){

	 				 			foreach my $level (keys %{$infoSequential->{$locusNameHIS}{$bucket}} ){
	 				 				push @{$infoSequential->{$locusNameUniq}{$idok}{$level}}, @{$infoSequential->{$locusNameHIS}{$bucket}{$level}};
	 				 			}
	 				 			delete $infoSequential->{$locusNameHIS}{$bucket};
	 				 			if(! %{$infoSequential->{$locusNameHIS}}){delete $infoSequential->{$locusNameHIS};} # remove because nothing linked to it anymore
	 				 		}
	 				 		else{
	 				 			$infoSequential->{$locusNameUniq}{$idok} = delete $infoSequential->{$locusNameHIS}{$bucket}; #delete the lod key but transfer the data to a new key
	 				 		}
	 				 	}
	 			 	}
	 			}
	 			else{ # The locusNameUniq didn't exists, we can directly shift the old locusNameHIS by the locusNameUniq
	 				$infoSequential->{$locusNameUniq} = delete $infoSequential->{$locusNameHIS}; # link to the first l2 #delete the old key but transfer the data to a new key
	 			}
			}
	 	}
	}
	dual_print($log, "$resume_case cases of interleaved locus. ( Features of the locus data are defined earlier in the file.\n", $verbose) if($resume_case);
}

#
#
# All Level1 feature are for sure in omniscient, we have only the ID in sequential,
# other level feature are in sequential only if no parent has been found
sub _check_sequential{ # Goes through from L3 to l1
 	my ($debug, $log, $infoSequential, $omniscient, $miscCount, $uniqID, $uniqIDtoType, $locusTAGuniq, $mRNAGeneLink, $verbose) = @_;
 	my $resume_case_l2 = undef;
	my %resume_case_l1;

 	_deinterleave_sequential($infoSequential, $locusTAGuniq, $verbose, $log); # PART OF LOCUS LOST BEFORE TO MEET ITS L2 or L1 ... we catch them and re-link everything as it should be

 	foreach my $locusNameHIS ( sort { ncmp($a,$b) } keys %{$infoSequential} ){ #comon tag was l1 id when no real comon tag present

 		foreach my $bucket (sort { ncmp($a,$b) } keys %{$infoSequential->{$locusNameHIS} } ){ #bucket = level1 or Id L2
 			dual_print($log, "\nlocusNameHIS $locusNameHIS bucket $bucket\n", $verbose) if($debug);

 			if ($bucket eq 'level1'){next;} #skip case level1 - structure of the hash different

 			my $must_create_l2=undef;
 			my $feature_l2 = undef;

			#Bucket is an uniq ID created during the reading process. So it can be used as uniq ID.
 			if(! exists_keys($infoSequential,($locusNameHIS, $bucket, 'level3') ) ){

	 				# Link the l2 to the L1 feature
					$feature_l2=$infoSequential->{$locusNameHIS}{$bucket}{'level2'};
					dual_print($log, "level2 in sequential does not have L3 feature associated in sequential".
						" - $locusNameHIS $bucket! ".$feature_l2->gff_string."\n", $verbose) if($debug);

					# We add it to omniscient and to mRNAGeneLink
					if(! exists($mRNAGeneLink->{lc($bucket)}) ){
						dual_print($log, "level2 does not exits in omniscient!".
															$feature_l2->gff_string."\n", $verbose) if($debug);

						push (@{$omniscient->{"level2"}{lc($feature_l2->primary_tag)}{lc($feature_l2->_tag_value('Parent'))} }, $feature_l2);
						$mRNAGeneLink->{lc($feature_l2->_tag_value('ID'))} = $feature_l2->_tag_value('Parent');
					}
					if( ! exists_keys($omniscient,('level3', "exon", lc($feature_l2->_tag_value("ID")))) ){ #check if an exon exist in the omniscient
						if ($createL3forL2orphan){
							# create the exon missing if option agreed
			 				dual_print($log, "create single level3 exon feature	!\n", $verbose);
			 				my $feature_l3 = clone($feature_l2);#create a copy of the l2 feature;
			 				$feature_l3->primary_tag('exon');
			 				create_or_replace_tag($feature_l3,'Parent', $feature_l3->_tag_value('ID')); # change parentID
							$feature_l3->remove_tag('ID');
							# create ID
							my $id =	_create_ID($miscCount, $uniqID, $uniqIDtoType, 'exon', undef, PREFIX_NEW_ID);
							create_or_replace_tag($feature_l3,'ID', $id); # change ID
							#my $id = _check_uniq_id($omniscient, $miscCount, $uniqID, $uniqIDtoType, $feature_l3);
							push (@{$omniscient->{"level3"}{lc($feature_l3->primary_tag)}{lc($feature_l3->_tag_value('Parent'))} }, $feature_l3);
						}
					}
	 				#warn "Not normal, we have feature L2	without L3 feature associated.\n";
					#We cannot guess the structure except if it is prokaryote or single exon in eucaryote.
 			}

			else{
 				foreach my $feature_L3 (@{$infoSequential->{$locusNameHIS}{$bucket}{'level3'}} ){
					if(! $feature_l2){
	 					if(! exists_keys($infoSequential,($locusNameHIS, $bucket,'level2'))	){
		 						dual_print($log, "_check_sequential level2 does not exits in sequential !\n", $verbose) if($debug);
								my $common_tag = _get_comon_tag_value($feature_L3, $locusTAGuniq, 'level1'); # check presence of common_tag, maybe we will play a different game

		 						#take L2 from omniscient if already exits
		 						if(exists($mRNAGeneLink->{lc($bucket)}) ){

		 							my $l1_id = $mRNAGeneLink->{lc($bucket)};
		 							foreach my $tag_l2 (keys %{$omniscient->{'level2'}}){
		 								if(exists_keys($omniscient, ('level2', $tag_l2, lc($l1_id) ) ) ){
				 							foreach my $featureL2 (@{$omniscient->{'level2'}{$tag_l2}{lc($l1_id)}}){
				 								if(lc($featureL2->_tag_value('ID')) eq lc($bucket) ){
				 									dual_print($log, "_check_sequential level2 exits in omniscient !\n", $verbose) if($debug);
				 									$feature_l2 = $featureL2;
				 									last;
				 								}
				 							}
				 							if($feature_l2){last;}
				 						}
			 						}
		 							$infoSequential->{$locusNameHIS}{$bucket}{'level2'} = $feature_l2;
								}
								#If locus_tag check from omniscient if feature has same locus tag
								elsif ( $common_tag and ( exists_keys($locusTAGuniq,('topfeature', $common_tag) ) ) and (! exists_keys($locusTAGuniq,('topfeature', $common_tag,'level1') ) ) ) {
										my $id_level2 = undef;
										if ( exists_keys($locusTAGuniq,('topfeature', $common_tag, 'level2') ) ){
											$id_level2 = $locusTAGuniq->{'topfeature'}{$common_tag}{'level2'}[0];
										}
										elsif ( exists_keys($locusTAGuniq,('topfeature', $common_tag, 'level3') ) ){
											$id_level2 = $locusTAGuniq->{'topfeature'}{$common_tag}{'level3'}[1];
										}
										dual_print($log, "FeatureA has the common tag value shared with a featureX from omniscient.".
										" We use same parent ID as featureX, and inject FeatureA in omniscient: $common_tag\n", $verbose) if($debug);

										create_or_replace_tag($feature_L3, 'Parent', $id_level2);
										push (@{$omniscient->{"level3"} {lc($feature_L3->primary_tag)} {$id_level2} }, $feature_L3);
										next;
								}
		 						else{#create l2
			 						$must_create_l2=1;
			 						$feature_l2 = clone($infoSequential->{$locusNameHIS}{$bucket}{'level3'}[0]);#create a copy of the first mRNA feature;

									#manage primary tag
									my $primary_tag_l2='RNA';
									foreach my $feature_L3 (@{$infoSequential->{$locusNameHIS}{$bucket}{'level3'}} ){

			 							if ( lc($feature_L3->primary_tag) eq 'cds'){
			 								$primary_tag_l2 ='mRNA';
			 								last;
			 							}
			 						}
			 						$feature_l2->primary_tag($primary_tag_l2);

		 							#Manage ID
									create_or_replace_tag($feature_l2,'ID', $bucket); #modify ID to replace by parent value
									dual_print($log, "_check_sequential ID created for level2: $bucket\n", $verbose) if($debug);

									#Manage Parent
									my $parentID = undef;
								 	if( exists_keys($infoSequential,($locusNameHIS,'level1'))	){ # parent ID exists in infoSequential
		 								$parentID = lc($infoSequential->{$locusNameHIS}{'level1'}); # Parent ID it correct case ???
		 								dual_print($log, "_check_sequential Parent ID for level2 taken from infoSequential: $parentID\n", $verbose) if($debug);
		 							}
									else{
										my $IDgoodCast = _id_exists_in_l1_omniscient($omniscient, $locusNameHIS);
										if($IDgoodCast){
												$parentID = $IDgoodCast;
												dual_print($log, "_check_sequential Parent ID for level2 taken from omniscient: $parentID\n", $verbose) if($debug);
										}

										if( ! $parentID ){ #In that case level1 feature doesn't exists in $infoSequential and in $omniscient. I will be created by the method check_gene_link_to_mrna
											#my	($miscCount, $uniqID, $primary_tag, $id, $prefix)=@_;
											$parentID =	_create_ID($miscCount, $uniqID, $uniqIDtoType, 'gene', undef, PREFIX_NEW_ID);
											dual_print($log, "_check_sequential Parent ID created for level2: $parentID\n", $verbose) if($debug);
											$infoSequential->{$locusNameHIS}{'level1'}=$parentID;
										}
									}
									create_or_replace_tag($feature_l2,'Parent', $parentID ); # change parentID
									push (@{$omniscient->{"level2"}{lc($primary_tag_l2)}{lc($parentID)}}, $feature_l2);
									$mRNAGeneLink->{lc($bucket)} = $parentID; # Always need to keep track about l2->l1, else the method check_l3_link_to_l2 will recreate a l1 thinking this relationship is not fill
									$infoSequential->{$locusNameHIS}{$bucket}{'level2'} = $feature_l2;
									dual_print($log, "feature level2 created: ".$feature_l2->gff_string."\n", $verbose);
									dual_print($log, "push-omniscient: level2 || ".lc($primary_tag_l2)." || ".lc($parentID)." == ".$feature_l2->gff_string."\n", $verbose) if($debug);
							 	}
	 					}
						else{
							#MUST push L2 in omniscient if absent !
							$feature_l2=$infoSequential->{$locusNameHIS}{$bucket}{'level2'};
							dual_print($log, "level2 exits in sequential - $locusNameHIS $bucket! ".$feature_l2->gff_string."\n", $verbose) if($debug);

							if(! exists($mRNAGeneLink->{lc($bucket)}) ){
								push (@{$omniscient->{"level2"}{lc($feature_l2->primary_tag)}{lc($feature_l2->_tag_value('Parent'))} }, $feature_l2);
								$mRNAGeneLink->{lc($feature_l2->_tag_value('ID'))} = $feature_l2->_tag_value('Parent');
								dual_print($log, "level2 does not exits in mRNAGeneLink(omniscient) !".$feature_l2->gff_string."\n", $verbose) if($debug);
							}
						}
					}

					my $primary_tag_L3 =	lc($feature_L3->primary_tag);
					create_or_replace_tag($feature_L3,'Parent', $feature_l2->_tag_value('ID')); #modify ID to replace by parent value
					push (@{$omniscient->{"level3"}{$primary_tag_L3}{lc($bucket)}}, $feature_L3);
					dual_print($log, "push-omniscient: level3 || ".$primary_tag_L3." || ".$bucket." == ".$feature_L3->gff_string."\n", $verbose) if($debug);
				}
			}

			$resume_case_l2++;
			#$resume_case_l1{lc( $feature_l2->_tag_value('Parent') )}++;

 			if($must_create_l2){
 				check_level2_positions($omniscient, $feature_l2);
			}
		}
		#LEVEL 1 IS taking care later
	}
	if(%resume_case_l1){
		my $nb_cases = scalar keys %resume_case_l1;
		dual_print($log, "We found $nb_cases sequential cases (holding $resume_case_l2 L2 features).\n", $verbose);
	}
	else{
		dual_print($log, "None found\n", $verbose);
	}
}

#print return ID if exists original cast
sub _id_exists_in_l1_omniscient{
	my ($omniscient, $id)=@_;

	my $id_good_cast=undef;
	foreach my $tag_l1 (keys %{$omniscient->{'level1'}} ){
		if(exists_keys($omniscient, ('level1',$tag_l1, $id))){
			$id_good_cast = $omniscient->{'level1'}{$tag_l1}{$id}->_tag_value('ID');
			return $id_good_cast;
		}
	}
	return $id_good_cast;
}

# L1: LocusID->level->typeFeature->ID->[ID,start,end]
# LocusID->level->typeFeature->Parent->[ID,start,end]
# @Purpose: When two feature overlap at level3, and are the same type level 2 they have to be merged under the same level 1 feature.
# @input: 2 =>	hash,	integer for verbosity
# @output: 0
sub _merge_overlap_features{
	my ($log, $omniscient, $mRNAGeneLink, $verbose) = @_;
	my $resume_case=undef;

	my $sortBySeq = _gather_and_sort_l1_by_seq_id_and_strand($omniscient);

	foreach my $locusID ( keys %{$sortBySeq}){ # tag_l1 = gene or repeat etc...
		foreach my $tag_l1 ( keys %{$sortBySeq->{$locusID}} ) {

			#create list to keep track of l1
			my %to_check;
			foreach my $feature_l1 ( @{$sortBySeq->{$locusID}{$tag_l1}} ) {
				my $id_l1 = lc($feature_l1->_tag_value('ID'));
				$to_check{$id_l1}++;
			}

			# Go through location from left to right ###
			while ( @{$sortBySeq->{$locusID}{$tag_l1}} ){

				my $feature_l1 = shift @{$sortBySeq->{$locusID}{$tag_l1}};
				my $id_l1 = lc($feature_l1->_tag_value('ID'));
				my @location = ($id_l1, int($feature_l1->start()), int($feature_l1->end())); # This location will be updated on the fly

				# Go through location from left to right ### !!
				foreach my $l1_feature2 ( @{$sortBySeq->{$locusID}{$tag_l1}} ) {
					my $id2_l1 = lc($l1_feature2->_tag_value('ID'));
					my @location_to_check = ($id2_l1, int($l1_feature2->start()), int($l1_feature2->end()));

					#If location_to_check start if over the end of the reference location, we stop
					if($location_to_check[1] > $location[2]) {last;}

					# Let's check at Gene LEVEL
					if(location_overlap(\@location, \@location_to_check)){

						#let's check at CDS level
						if(check_gene_overlap_at_CDSthenEXON($omniscient, $omniscient , $id_l1, $id2_l1)){ #If contains CDS it has to overlap at CDS level to be merged, otherwise any type of feature level3 overlaping is sufficient to decide to merge the level1 together
							#they overlap in the CDS we should give them the same name
							$resume_case++;

							dual_print($log, "$id_l1 and $id2_l1 same locus. We merge them together: Below the two features:\n".$feature_l1->gff_string."\n".$l1_feature2->gff_string."\n", 0); # print only in log
							# remove the level1 of the ovelaping one
							delete $omniscient->{'level1'}{$tag_l1}{$id2_l1};
							# remove the level2 to level1 link stored into the mRNAGeneLink hash. The new links will be added just later after the check to see if we keep the level2 feature or not (we remove it when identical)
							foreach my $l2_type (%{$omniscient->{'level2'}}){
								if(exists_keys($omniscient,('level2', $l2_type, $id2_l1))){
									foreach my $feature_l2 (@{$omniscient->{'level2'}{$l2_type}{$id2_l1}}){
										delete $mRNAGeneLink->{lc($feature_l2->_tag_value('ID'))};
									}
								}
							}

							# Let's change the parent of all the L2 features
							foreach my $l2_type ( keys	%{$omniscient->{'level2'}} ){

								if(exists_keys($omniscient,('level2', $l2_type, $id2_l1))){
									###############################
									# REMOVE THE IDENTICAL ISOFORMS

									# first list uniqs
									my $list_of_uniqs	= keep_only_uniq_from_list2($omniscient, $omniscient->{'level2'}{$l2_type}{$id_l1}, $omniscient->{'level2'}{$l2_type}{$id2_l1}, $verbose); # remove if identical l2 exists


									#Now manage the rest
									foreach my $feature_l2 (@{$list_of_uniqs}){

										create_or_replace_tag($feature_l2,'Parent', $feature_l1->_tag_value('ID')); #change the parent
										# Add the corrected feature to its new L2 bucket
										push (@{$omniscient->{'level2'}{$l2_type}{$id_l1}}, $feature_l2);

										# Attach the new parent into the mRNAGeneLink hash
										$mRNAGeneLink->{lc($feature_l2->_tag_value('ID'))}=$feature_l2->_tag_value('Parent');

									}
									# remove the old l2 key
									delete $omniscient->{'level2'}{$l2_type}{$id2_l1};
								}
							}
							check_level1_positions( { omniscient => $omniscient, feature => $omniscient->{'level1'}{$tag_l1}{$id_l1} } );
						}
					}
				}
	 		}
	 	}
	}
	if($resume_case){
		dual_print($log, "$resume_case overlapping cases found. For each case 2 loci have been merged within a same locus\n", $verbose);
	}
	else{
		dual_print($log, "None found\n", $verbose);
	}
}


# @Purpose: When too feature l2 isoform are identical, we remove one
# @input: 2 =>	hash,	integer for verbosity
# @output: 0
sub _check_identical_isoforms{
	my ($log, $omniscient, $mRNAGeneLink, $verbose) = @_;
	my $resume_case=undef;

	# Go through oall l2 feature
	foreach my $l2_type (keys %{$omniscient->{'level2'}}){
		foreach my $id2_l1 (keys %{$omniscient->{'level2'}{$l2_type}}){
			# If more than 1 related to level1

			if(exists_keys($omniscient,('level2', $l2_type, $id2_l1)) and scalar @{$omniscient->{'level2'}{$l2_type}{$id2_l1}} > 1){ # more than one l2 feature of that type

				my @L2_list_to_remove;
				my %checked;
				foreach my $feature2 (sort {$b->_tag_value('ID') cmp $a->_tag_value('ID')} @{$omniscient->{'level2'}{$l2_type}{$id2_l1}}){
					$checked{lc($feature2->_tag_value('ID'))}{lc($feature2->_tag_value('ID'))}++;

					my $keep = 1;
					foreach my $feature1 (sort {$b cmp $a} @{$omniscient->{'level2'}{$l2_type}{$id2_l1}}){

						# If not itself and not already checked (A -> B is the same as B -> A), and A or B already removed and must now be skiped (skipme key)
						if( (! exists_keys(\%checked, (lc($feature2->_tag_value('ID')), "skipme"))) and (! exists_keys(\%checked, (lc($feature1->_tag_value('ID')), "skipme"))) and ! exists_keys(\%checked, ( lc($feature2->_tag_value('ID')), lc($feature1->_tag_value('ID')) ) ) ){ #
							$checked{lc($feature2->_tag_value('ID'))}{lc($feature1->_tag_value('ID'))}++;
							$checked{lc($feature1->_tag_value('ID'))}{lc($feature2->_tag_value('ID'))}++;

							#check their position are identical
							if($feature1->start().$feature1->end() eq $feature2->start().$feature2->end()){

								#Check their subfeature are	identicals
								if(l2_identical($omniscient, $feature1, $feature2, $verbose )){
									$keep = undef;
									last;
								}
							}
						}
					}
					# We dont keep the l2 feature so we have to remove all related features and itself
					if(! $keep){
						$resume_case++;
						dual_print($log, "Lets remove isoform ".$feature2->_tag_value('ID')."\n",$verbose);
						$checked{lc($feature2->_tag_value('ID'))}{"skipme"}++;# will be removed later do not check anymore this one

						foreach my $tag (keys %{$omniscient->{'level3'}}){
							if(exists_keys($omniscient, ('level3', $tag, lc($feature2->_tag_value('ID'))))){
								delete $omniscient->{'level3'}{$tag}{lc($feature2->_tag_value('ID'))};
							}
						}
						#Has to be removed once we finished to go through the l2 list
						my $ID_to_remove = lc($feature2->_tag_value('ID'));
						push(@L2_list_to_remove,$ID_to_remove);
						delete $mRNAGeneLink->{$ID_to_remove};
					}
				}

				#L2 has to be removed from List
				my @newL2List;
				foreach my $feature ( @{$omniscient->{'level2'}{$l2_type}{$id2_l1}} ){
					my $keep = 1;
					foreach my $id_l2 (@L2_list_to_remove){
						if( lc($feature->_tag_value('ID')) eq lc($id_l2) ){
							$keep = undef;
						}
					}
					if($keep){
						push (@newL2List,$feature)
					}
				}
				@{$omniscient->{'level2'}{$l2_type}{$id2_l1}}=@newL2List;
			}
		}
	}
	if($resume_case){
		dual_print($log, "$resume_case identical isoforms removed\n", $verbose);
	}
	else{dual_print($log,"None found\n", $verbose)}
}

# Sort by locusID and strand
# LocusID_strand->typeFeature = [feature, feature, feature]
# return a hash. Key is position,tag and value is list of feature l1. The list is sorted
sub _gather_and_sort_l1_by_seq_id_and_strand{
	my ($omniscient) = @_;

	my %hash_sortBySeq;
		foreach my $tag_level1 (keys %{$omniscient->{'level1'}}){
				foreach my $level1_id (keys %{$omniscient->{'level1'}{$tag_level1}}){
					my $level1_feature = $omniscient->{'level1'}{$tag_level1}{$level1_id};
					my $position_l1=$level1_feature->seq_id.$level1_feature->strand;
					push (@{$hash_sortBySeq{$position_l1}{$tag_level1}}, $level1_feature);
				}
				foreach my $position_l1 (keys %hash_sortBySeq){
					@{$hash_sortBySeq{$position_l1}{$tag_level1}} = sort { ncmp ($a->start.$a->end.$a->_tag_value('ID'), $b->start.$b->end.$b->_tag_value('ID') ) } @{$hash_sortBySeq{$position_l1}{$tag_level1}};
				}
		}
	return \%hash_sortBySeq;
}

#
sub modelate_utr_and_cds_features_from_exon_features_and_cds_start_stop{

	my ($exon_features, $ORFstart, $ORFend)=@_;

	my @cds_features;
	my @utr3_features;
	my @utr5_features;
	my $strand = $exon_features->[0]->strand;
		my @exon_features_sorted = sort {$a->start <=> $b->start} @{$exon_features}; # be sure that exon list is sorted

		my $cds_counter=1;
		my $utr3_counter=1;
		my $utr5_counter=1;
 	foreach my $exon_feature (@exon_features_sorted){

			# exon overlap fully a CDS
			if( ($exon_feature->end >= $ORFend) and ($exon_feature->start <= $ORFstart) ){

 			my $cds_feature=clone($exon_feature);#create a copy of the feature 					exon		====================================
 			$cds_feature->start($ORFstart); #modify start 											 cds		 ============================
 			$cds_feature->end($ORFend); #modify end
 			$cds_feature->primary_tag('CDS');
 			#get old name
 			my $ID = $cds_feature->_tag_value('ID');
 			create_or_replace_tag($cds_feature,'ID',$ID.'-cds-'.$cds_counter); #modify name
 			push(@cds_features, $cds_feature);#save that cds
 			$cds_counter++;
 			if($exon_feature->start < $ORFstart){
 				my $utr_feature=clone($exon_feature);#create a copy of the feature
 				$utr_feature->end($ORFstart-1); #modify start
 				if ( ($strand == -1) or ($strand eq "-") ) {
 					$utr_feature->primary_tag('three_prime_UTR');
 					create_or_replace_tag($utr_feature,'ID',$ID.'-utr3-'.$utr3_counter); #modify name
	 				push(@utr3_features, $utr_feature);#save that cds
	 				$utr3_counter++;
	 			}else{
	 				$utr_feature->primary_tag('five_prime_UTR');
	 				create_or_replace_tag($utr_feature,'ID',$ID.'-utr5-'.$utr5_counter); #modify name
	 				push(@utr5_features, $utr_feature);#save that cds
	 				$utr5_counter++;
	 			}
 			}
 			if($exon_feature->end > $ORFend){
 				my $utr_feature=clone($exon_feature);#create a copy of the feature
 				$utr_feature->start($ORFend+1); #modify start
 				if ( ($strand == -1) or ($strand eq "-") ) {
 					$utr_feature->primary_tag('five_prime_UTR');
	 				create_or_replace_tag($utr_feature,'ID',$ID.'-utr5-'.$utr5_counter); #modify name
	 				push(@utr5_features, $utr_feature);#save that cds
	 				$utr5_counter++;
	 			}else{
	 				$utr_feature->primary_tag('three_prime_UTR');
 					create_or_replace_tag($utr_feature,'ID',$ID.'-utr3-'.$utr3_counter); #modify name
	 				push(@utr3_features, $utr_feature);#save that cds
	 				$utr3_counter++;
	 			}
 			}
		}
		# cds overlap fully an exon
		elsif( ($exon_feature->end <= $ORFend) and ($exon_feature->start >= $ORFstart) ){
 			my $cds_feature=clone($exon_feature);#create a copy of the feature 						exon		========================
 			$cds_feature->primary_tag('CDS');
 			#get old name 																			cds	===============================
 			my $ID = $cds_feature->_tag_value('ID');
 			create_or_replace_tag($cds_feature,'ID',$ID.'-cds-'.$cds_counter); #modify name
 			push(@cds_features, $cds_feature);#save that cds
 			$cds_counter++;
		}
		# cds overp partially an exon
			elsif( ($exon_feature->end >= $ORFstart) and ($exon_feature->start <= $ORFend) ){ #they overlap

				if($exon_feature->start >= $ORFstart){ # cds overlap start of exon																		exon ===============================
					#Manage CDS
					my $cds_feature=clone($exon_feature);#create a copy of the feature 						cds ===============================
 			$cds_feature->end($ORFend); #modify end
 			$cds_feature->primary_tag('CDS');
 				#get old name
 			my $ID = $cds_feature->_tag_value('ID');
 			create_or_replace_tag($cds_feature,'ID',$ID.'-cds-'.$cds_counter); #modify name
 			push(@cds_features, $cds_feature);#save that cds
 			$cds_counter++;
 			#manage UTR
 			my $utr_feature=clone($exon_feature);#create a copy of the feature
 			$utr_feature->start($ORFend+1); #modify end
 			$ID = $utr_feature->_tag_value('ID');
	 		if ( ($strand == -1) or ($strand eq "-") ) {
	 			$utr_feature->primary_tag('five_prime_UTR');
	 			create_or_replace_tag($utr_feature,'ID',$ID.'-utr5-'.$utr5_counter); #modify name
	 			push(@utr5_features, $utr_feature);#save that cds
	 			$utr5_counter++;

	 		}else{
				$utr_feature->primary_tag('three_prime_UTR');
	 			create_or_replace_tag($utr_feature,'ID',$ID.'-utr3-'.$utr3_counter); #modify name
	 			push(@utr3_features, $utr_feature);#save that cds
	 			$utr3_counter++;
	 		}

				}
				else{ #cds overlap start end exon
				 	#Manage CDS
				 	my $cds_feature=clone($exon_feature);#create a copy of the feature
 			$cds_feature->start($ORFstart); #modify start 										exon ===============================
 			$cds_feature->primary_tag('CDS');
 				#get old name 																					cds =====================================
 			my $ID = $cds_feature->_tag_value('ID');
 			create_or_replace_tag($cds_feature,'ID',$ID.'-cds-'.$cds_counter); #modify name
 			push(@cds_features, $cds_feature);#save that cds
 			$cds_counter++;
	 		 #Manage UTR
 			my $utr_feature=clone($exon_feature);#create a copy of the feature
 			$utr_feature->end($ORFstart-1); #modify start
 			$ID = $utr_feature->_tag_value('ID');
	 		if ( ($strand == -1) or ($strand eq "-") ) {
	 			$utr_feature->primary_tag('three_prime_UTR');
	 			create_or_replace_tag($utr_feature,'ID',$ID.'-utr3-'.$utr3_counter); #modify name
	 			push(@utr3_features, $utr_feature);#save that cds
	 			$utr3_counter++;
	 		}else{
	 			$utr_feature->primary_tag('five_prime_UTR');
	 			create_or_replace_tag($utr_feature,'ID',$ID.'-utr5-'.$utr5_counter); #modify name
	 			push(@utr5_features, $utr_feature);#save that cds
	 			$utr5_counter++;
	 		}
				}
			}###### Only UTR part
			else{ #Does not overlap
				if($exon_feature->end < $ORFstart){ #UTR5 in + strand
					my $utr_feature=clone($exon_feature);#create a copy of the feature 			exon ===============================
	 			#get old name 																											 cds ===============================
	 			my $ID = $utr_feature->_tag_value('ID');
	 			if ( ($strand == -1) or ($strand eq "-") ) {
	 				$utr_feature->primary_tag('three_prime_UTR');
	 				create_or_replace_tag($utr_feature,'ID',$ID.'-utr3-'.$utr3_counter); #modify name
	 				push(@utr3_features, $utr_feature);#save that cds
	 				$utr3_counter++;
	 			}else{
	 				$utr_feature->primary_tag('five_prime_UTR');
	 				create_or_replace_tag($utr_feature,'ID',$ID.'-utr5-'.$utr5_counter); #modify name
	 				push(@utr5_features, $utr_feature);#save that cds
	 				$utr5_counter++;
	 			}

	 		}
	 		else{	#UTR3 in + strand
					my $utr_feature=clone($exon_feature);#create a copy of the feature 													exon ===============================
	 			#get old name
	 			my $ID = $utr_feature->_tag_value('ID'); 									#cds ===============================
	 			if ( ($strand == -1) or ($strand eq "-") ) {
	 				$utr_feature->primary_tag('five_prime_UTR');
	 				create_or_replace_tag($utr_feature,'ID',$ID.'-utr5-'.$utr5_counter); #modify name
	 				push(@utr5_features, $utr_feature);#save that cds
	 				$utr5_counter++;
	 			}else{
	 				$utr_feature->primary_tag('three_prime_UTR');
	 				create_or_replace_tag($utr_feature,'ID',$ID.'-utr3-'.$utr3_counter); #modify name
	 				push(@utr3_features, $utr_feature);#save that cds
	 				$utr3_counter++;
	 			}
	 		}
			}
 	}
my @utr5_features_sorted=sort {$a->start <=> $b->start} @utr5_features;
my @cds_features_sorted=sort {$a->start <=> $b->start} @cds_features;
my @utr3_features_sorted=sort {$a->start <=> $b->start} @utr3_features;
return \@utr5_features_sorted, \@cds_features_sorted, \@utr3_features_sorted; #really utr5 and utr3 that are return
}


# Actually the duplicates have been collected during the parsing process here we just print them.
sub _check_duplicates{
	my ($log, $duplicate, $omniscient, $verbose) = @_	;

		my $keyExist = keys %{$duplicate};
		if($keyExist){#print result
			dual_print ($log, "Duplicates detected! (Same chr/contig/scaffold, same position, same ID)\n", $verbose);

			my $string="";
			foreach my $level (keys %{$duplicate}){ # primary_tag_key_level1 = gene or repeat etc...
				foreach my $primary_tag (keys %{$duplicate->{$level}}){
					my $nb_by_pt=0;
					my $nb_feat_pt=0;
					foreach my $id (keys %{$duplicate->{$level}{$primary_tag}}){
						$nb_feat_pt++;
						foreach my $feature (@{$duplicate->{$level}{$primary_tag}{$id}}){
							$nb_by_pt++;
							print $log $feature->gff_string."\n" if $log;  # print feature only in log
						}
					}
					dual_print( $log, "Removed $nb_feat_pt duplicated $primary_tag features for a total of $nb_by_pt duplicates.\n", $verbose);
				}
			}
		}
		else{
			dual_print( $log, "None found\n", $verbose);
		}
}

# Method to store all headers (before the first feature)
# Input: filename
# Output: string (header lines)
sub get_header_lines{
	my ($file, $verbose, $log) = @_;

	#HANDLE format
	my @headers;

	open(my $fh, '<', $file) or dual_print($log, "cannot open file $file", 1) && die;
	{
		while(<$fh>){
			if($_ =~ /^#/){
				if($_ =~ /##gff-version/){next;}# we do not keep the version line because we will write it ourself
				push @headers, $_;
				dual_print($log, "catch header line: $_", 0); # print in log only
			} #if it is a commented line starting by # we skip it.
			else{
				close($fh);
				return \@headers;
			}
		}
	}
	close($fh);
	return \@headers;
}

#GFF format guess
# Input: filename
# Output: Integer (1,2 or 3)
sub select_gff_format{
		my ($file, $verbose, $log) = @_;

		#HANDLE format
		my %format;
		my $problem3=undef;
		my $nbLineChecked=100; #number line to use to check the formnat
		my $cpt=0;

		open(my $fh, '<', $file) or dual_print($log, "cannot open file $file", 1) && die;
		{
			while(<$fh>){

				if($_ =~ /^#/){next;} #if it is a commented line starting by # we skip it.

				$cpt++;
				if($cpt > $nbLineChecked){
								last;
				}
				if($_ =~ /^.*\t.*\t.*\t.*\t.*\t.*\t.*\t.*\t(.*)/){
					if(length($1) < 1){next;}

					my $Ninethcolum = $1;
					if($Ninethcolum =~ /=/	and $Ninethcolum =~ /;/ ){ $format{3}++;};

					if($Ninethcolum !~ /=/	and $Ninethcolum !~ /;/ ){
									 $format{1}++;
					}
					elsif($Ninethcolum !~ /=/	and $Ninethcolum =~ /;/ ){
													 $format{2}++;
					}
					my $c = () = $Ninethcolum =~ /=/g;
					my $d = () = $Ninethcolum =~ /\ /g;
					if($c > 1 and $d > 1	and $Ninethcolum !~ /;/ ){
								 $problem3=1;
					}
	 			}
			}
		}
		close($fh);

	if($problem3){
				dual_print ($log, surround_text("There is a problem with your GFF format.\nThis format is wrong: tag=value tag=value.\nYou should have: tag=value;tag=value or tag value ; tag value\nThe best parser (gff1) we can use will keep only the first attribute.",100,"!"));
			$format{1}++;
		}

	 if (%format){
			my $number_of_format = scalar keys %format;
			if ($number_of_format > 1){
				my $stringprint = "There is a problem we found several formats in this file:\n";
				$stringprint .= join ",", sort keys %format;
				$stringprint .= "\nLet's see what we can do...\n";
				print $stringprint if ($verbose);
		}
	}
	else{
		dual_print ($log, surround_text("Doesn't look like a GFF file\nLet's see what the Bioperl parser can do with that...(using gff3 parser)",100,"!") );
		$format{3}++;
	}
	if($format{3}){return 3;}
	if($format{2}){return 2;}
	if($format{1}){return 1;}
}

# We modify the attributes: group=gene_id "e_gw1.5.2.1" protein_id 335805 exonNumber 1
# in order to get : gene_id=e_gw1.5.2.1;protein_id=335805;exonNumber=1
sub _gff1_corrector{
	my ($feat, $verbose)=@_;

	if($feat->has_tag('group')){

		my @attribs = $feat->get_tag_values('group');
		my $attribs = join ' ', @attribs;
		my @parsed;
		my $flag = 0; # this could be changed to a bit and just be twiddled

			# run through each character one at a time and check it
			my $previousChar=undef;
			my $string="";
			foreach my $a ( split //, $attribs ) {
					$string.=$a;

					# flag up on entering quoted text, down on leaving it
					if( $a eq '"') { $flag = ( $flag == 0 ) ? 1:0 ;} #active deactive the flag

					if ($previousChar and $previousChar eq '"' and $flag == 0){ # case we have to strip the " characters
						chop $string;
						chop $string;
						$string = reverse($string);
						chop($string);
						$string= reverse($string);
						push @parsed, $string;
						$string="";
					}
					elsif( ( $a eq " " and $flag == 0) and !($string =~ /^ *$/) ){
						chop $string;
						push @parsed, $string;
						$string="";
					}
					$previousChar = $a;
			}
			# ---- Check now last string ----
			# If it was quoted
			if ($previousChar and $previousChar eq '"' and $flag == 0){ # case we have to strip the " characters
				chop $string;
				$string = reverse($string);
				chop($string);
				$string= reverse($string);
				push @parsed, $string;
			}# If it not empty or not only space and not quoted
			elsif( ($string ne "") and !($string =~ /^ *$/)	){
				if($previousChar eq " "){
					chop $string;
				}
				push @parsed, $string;
			}

			while (@parsed){
				my $value = pop @parsed;
				my $tag = pop @parsed;
				$feat->add_tag_value($tag, $value);
			}
		#remove old group attribute
		$feat->remove_tag('group');
		}
}

# @Purpose: Create a hash containing all the name and identifier of an ontology.
# @input: 1 =>	Object Bio::Ontology
# @output: 1 => hash containing all the name and identifier
sub create_term_and_id_hash{
		my ($ontology, $verbose, $log, $debug) = @_;

		my %hash_term_id;

		# Print some information
		dual_print($log, "	Filtering ontology:\n", $verbose);
		dual_print($log, "The feature type (3rd column) is constrained to be either a term from the Sequence ".
		"Ontology or an SO accession number. The latter alternative is distinguished".
		" using the syntax SO:000000. In either case, it must be sequence_feature ".
		"(SO:0000110) or an is_a child of it.\nWe filter the ontology to apply this rule.", 0) ; # print in log only

		#Get top term
		my ($term) = $ontology->find_terms(-name => "sequence_feature");
		# get all descendant and the top term
		my @descendants = $ontology->get_descendant_terms($term);
		dual_print($log, "		We found ".(scalar @descendants)." terms that are sequence_feature or descendant terms of it (is_a constraint not applied).\n", $verbose) if($debug);
		# get a is_a relationship
		my $is_a = Bio::Ontology::RelationshipType->get_instance('IS_A');
		@descendants = $ontology->get_descendant_terms($term, $is_a);
		dual_print($log, "		We found ".(scalar @descendants)." terms that are sequence_feature or is_a child of it.\n", $verbose);

		foreach my $term (@descendants) {
			 $hash_term_id{lc($term->name)} = lc($term->identifier);
			 $hash_term_id{lc($term->identifier)} = lc($term->name);
			 #print $term->name." <=> ".$term->identifier."\n";
		}
		return \%hash_term_id;
}

#Look for gff3 specific header
#@INPUT: 1 => string (a file)
#@OUPUT: 1 => hash of the different header and their values
sub _check_header{
		my ($file, $log) = @_;

		#HANDLE format
		my %headerInfo;

		#check it is a file
		if(-f $file){
			open(my $fh, '<', $file) or dual_print($log, "cannot open file $file", 1) && die;
			{
					while(<$fh>){
							if($_ !~ /^##[^#]/) {
									 last;
							}
							else{
								my @data = split /\s/, $_ ;
								my $type = shift @data;

								if($type eq /^##gff-version/){
									$headerInfo{$type}=$data[0]; #1 element
								}
					if($type eq "##sequence-region"){
						$headerInfo{$type}=@data; # 3 elements
					}
					if($type eq "##feature-ontology"){
						$headerInfo{$type}=$data[0] #1 element
					}
					if($type eq "##attribute-ontology"){
						$headerInfo{$type}=$data[0]; #1 element
					}
					if($type eq "##species"){
						$headerInfo{$type}=$data[0]; #1 element
					}
					if($type eq "##genome-build"){
						$headerInfo{$type}=@data; #2 elements
					}
						}
					}
			}
			close($fh);
	}

		return \%headerInfo;
}

# @Purpose: Read a file from URL
# @input: 2 =>	String URL, String target (Target is not mandatory)
# @output: none
sub fetcher_JD {
	my ($url, $target, $log) = @_;
		my $ua = LWP::UserAgent->new;
		$ua->timeout(10);
		$ua->env_proxy;

		my $response = $ua->get($url);
		if ($response->is_success) {
			if($target){
			 	open my $OUT, '>', $target or dual_print($log, "File error: $! :: $?", 1) && die;
					print $OUT $response->decoded_content;	# or whatever
			}
			else{
				my $string = $response->decoded_content;
				return $string ;
			}
		}
		else {
				dual_print($log, $response->status_line, 1) && die;
		}
}

# @Purpose: retrieve the feature_ontology
# @input: 3 =>	String file, Hash, Int
# @output: 1 => Object Ontology
# @Remark: Do not deal if multiple ontologies (we will use the first one meet)
sub _handle_ontology{
	my ($gff3headerInfo, $verbose, $log) = @_ ;

	dual_print( $log, "=> Accessing Ontology\n", $verbose );
	my $ontology_obj=undef;
	my $internalO=1;

		if(exists_keys($gff3headerInfo, ("##feature-ontology"))){

			dual_print( $log, "	feature-ontology URI defined within the file: ".$gff3headerInfo->{'##feature-ontology'}."\n", $verbose);
			#retrieve the data from URI and save it in a string
			my $stringFILE=undef;
			try{
				$stringFILE = fetcher_JD($gff3headerInfo->{"##feature-ontology"}, undef, $log);
			}
			catch{
				dual_print( $log, "	The URI provided (".$gff3headerInfo->{'##feature-ontology'}.") doesn't work.\n", $verbose);
				dual_print( $log, "error: $_\n", $verbose);
			};

			if($stringFILE){
				#create a filehandler from a string
				open( my $fh_uriOnto, '<', \$stringFILE) or
								dual_print( $log, "Cannot read the string: $! :: $?", 1) && die;

				#parse the ontology saved
		 		my $parser = undef;
		 		try{
		 			$parser = Bio::OntologyIO->new(-format => "obo",
																						 -fh => $fh_uriOnto);
		 			$ontology_obj = $parser->parse();
		 			close $fh_uriOnto;
		 		}
		 		catch{
		 			dual_print( $log, "	The URI provided doesn't point to obo ontology format data.\n", $verbose );
		 			dual_print( $log, "error: $_\n", $verbose);
		 			$parser = undef;
		 		};

				if($parser){ #We got ontology at the URI location, no need to use the internal one
					$internalO=undef;
					dual_print( $log, "	feature-ontology parsed correctly\n", $verbose );
				}
			}
		}

	if($internalO){ #No URI provided for the feature-ontology(file case), or doesn't exist (hash / table case) let's use the interal one
		dual_print( $log, "	No ontology accessible from the gff file header!\n", $verbose);
		try{
			my $sofa_file_path = dist_file('AGAT', 'so.obo');
			dual_print( $log, "	We use the SOFA ontology distributed with AGAT:\n		$sofa_file_path\n", $verbose);

			#parse the ontology
			my $parser = Bio::OntologyIO->new(-format => "obo",
																				-file => $sofa_file_path);
			$ontology_obj = $parser->parse();
			if($verbose) {
				my $nbroot_terms =0;
				foreach my $term ($ontology_obj->get_root_terms) {
					$nbroot_terms++;
				}
				my $nbterms =0;
				foreach my $term ($ontology_obj->get_all_terms) {
					$nbterms++;
				}
				my $nbleaf_terms =0;
				foreach my $term ($ontology_obj->get_leaf_terms) {
					$nbleaf_terms++;
				}
				dual_print( $log, "	Read ontology $sofa_file_path:\n".
						 	"		$nbroot_terms root terms, and ".
						 	"$nbterms total terms, and ".
						 	"$nbleaf_terms leaf terms\n", $verbose );
			}
		}
		catch{
			dual_print($log, "error: $_\n", $verbose );
 			dual_print($log, "	Let's continue without feature-ontology information.\n", $verbose );
		};
	}
	return $ontology_obj;
}

# @Purpose: Handle global warnings to provide more information to the user according
# to problems encountered. Global warning is related to feature type (3rd columm)
# if problem with agat parser, or with ontology
# @input: 3 =>	hash,
# @output: 1 => none (because it will just display infromation)
# @Remark: none
sub _handle_globalWARNS{
	my ($args) = @_;

	# -------------- OUTPUT --------------
	my $result = "";

	# -------------- INPUT --------------
	# Check we receive a hash as ref
	if(ref($args) ne 'HASH'){ warn "Hash Arguments expected for _handle_globalWARNS. Please check the call.\n";exit;	}
	# -- Declare all variables and fill them --
	my ($globalWARNS, $ontology, $log, $type, $verbose);

	if( defined($args->{warning})) {$globalWARNS = $args->{warning};} else{ $globalWARNS = undef;}
	if( defined($args->{ontology})) {$ontology = $args->{ontology};} else{ $ontology = undef;}
	if( defined($args->{log})) {$log = $args->{log};} else{ $log = undef;}
	if( defined($args->{type})) {$type = $args->{type};} else{ $type = undef;}
	if( defined($args->{verbose})) {$verbose = $args->{verbose};} else{ $verbose = undef;}

	if( lc($type) eq "ontology" ){
		dual_print($log, file_text_line({ string => "ontology", char => "-" }), $verbose);
		my $string;
		if(exists_keys( $globalWARNS, ("ontology1") ) ) {
			if( keys %{$ontology} ){
				my %hash	 = map { $_, 1 } @{$globalWARNS->{ontology1}};
				my @unique = sort keys %hash;
				$string = "INFO - Feature types not expected by the GFF3 specification:\n* ".
				join("\n* ", @unique).
				"\nThe feature type (3rd column in GFF3) is constrained to be either a term from the Sequence Ontology ".
				"or an SO accession number. The latter alternative is distinguished using the ".
				"syntax SO:000000. In either case, it must be sequence_feature (SO:0000110) or ".
				"an is_a child of it.".
				"To follow rigorously the gff3 format, please visit this website:\n".
				"https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md\n";
			}
			else{
				$string = "No ontology was available, we haven't checked if the feature types (3rd column) correspond to the gff3 specifactions.";
			}
		}
		else{
			$string = "All feature types in agreement with the Ontology.";
		}
		dual_print ($log, print_wrap_text($string,80), $verbose);
	}
	if( lc($type) eq "agat" ){
		dual_print($log, file_text_line({ string => "agat", char => "-" }), $verbose);
		my $string;
		if(exists_keys($globalWARNS, ("parser1") ) ) {
			my %hash	 = map { $_, 1 } @{$globalWARNS->{parser1}};
			my @unique = sort keys %hash;
			$string = "WARNING - Feature types not expected by AGAT:\n* ".
			join("\n* ", @unique).
			"\nThe feature of these types (3rd column in GFF3) are skipped by the parser!\n".
			"To take them into account you must update the feature json files. To access the json files run:".
			"\n			agat_convert_sp_gxf2gxf.pl --expose\n".
			"In which file to add my feature?\n".
			"* Feature level1 (e.g. gene, match, region):\n  My feature has no parent\n  => features_level1.json\n".
			"* Feature level2 (e.g. mrna, match_part, trna):\n  My feature has one parent and children\n  => features_level2.json.\n".
			"* Feature level3 (e.g. exon, intron, cds):\n  My feature has one parent (the parent has also a parent) and no children\n  => features_level3.json.\n".
			"* Feature level3 discontinuous (e.g. cds, utr):\n  A single feature that exists over multiple genomic locations\n  => features_spread.json.";
		}
		else{
			$string = "AGAT can deal with all the encountered feature types (3rd column)";
		}
		dual_print ($log, print_wrap_text($string,80), $verbose);
	}
}

1;
