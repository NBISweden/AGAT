#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use AGAT::AGAT;

start_script();
my $header = get_agat_header();
# -------------------------------- LOAD OPTIONS --------------------------------
my $mfannot_file;
my $gff_file;
my %startend_hash;     # Stores start and end positions of each feature reported
my %sorted_hash;
my %hash_uniqID;
my %filtered_result;
my $omniscient={}; #Hash where all the features will be saved
my $hashID={}; # ex %miscCount;# Hash to store any counter.

# OPTION MANAGEMENT: partition @ARGV into shared vs script options via library
my ($shared_argv, $script_argv) = split_argv_shared_vs_script(\@ARGV);

# Parse script-specific options from its own list
my $script_parser = Getopt::Long::Parser->new;
$script_parser->configure('bundling','no_auto_abbrev');
if ( !$script_parser->getoptionsfromarray(
				$script_argv,
				'mfannot|m|i=s'  => \$mfannot_file,
				'o|out|output=s' => \$gff_file,
				'h|help'         => sub { pod2usage( -exitstatus=>0, -verbose=>99, -message => "$header\n" ); },
				'man'            => sub { pod2usage(-exitstatus=>0, -verbose=>2); }
			))
{
    pod2usage( { -message => "Failed to parse command line.\n",
                 -verbose => 1,
                 -exitval => 1 } );
}
if ( ! (defined($mfannot_file)) ){
    pod2usage( {
           -message => "$header\nAt least 1 parameter is mandatory:\nInput mfannot file (--mfannot).\n\n",
           -verbose => 0,
           -exitval => 1 } );
}
# Parse shared options without pass_through for strong type errors. CPU and config are handled there.
my ($shared_opts) = parse_shared_options($shared_argv);

# --- Load config file into global CONFIG ---
initialize_agat({config_file_in => ( $shared_opts->{config} ), input => $mfannot_file, shared_opts => $shared_opts });

# ------------------------------------------------------------------------------

## Manage output file
my $gffout = prepare_gffout( $gff_file );

## MAIN ##############################################################
read_mfannot($mfannot_file);

handle_records();

my ($hash_omniscient) = slurp_gff3_file_JD({ input => $omniscient });

print_omniscient( {omniscient => $hash_omniscient, output => $gffout} );

# --- final messages ---
end_script();

## SUBROUTINES #######################################################

sub read_mfannot {
    my $current_contig;         # Track the current contig
    my $current_genetic_code;   # Track current genetic code
    my $current_pos=1;          # Track current position
    my $writeflag=0;
    my $previousDirection=undef;
    my $previousStartEnd=undef;;
	my $previousIntron=undef;
	my $previousRnl=undef;
	my $previousRns=undef;
	my $position=0;

    open(INPUT, "<", "$_[0]") or die ("$!\n");
    # Open Mfannot file for reading
    while (<INPUT>) {
        chomp;
		#print "reading $_  ++\n";
        if ($_ =~ /^>(.*) gc=(\d+)/) {
            # If a header line, update the current contig and genetic code
            ($current_contig, $current_genetic_code) = ($1, $2);
            $current_pos=1; # Reset the position counter
			if (! exists_keys($omniscient,("other","header")) ){
				push @{$omniscient->{"other"}{"header"}}, "##transl_table=".$current_genetic_code;
			}
		}
        elsif ($_ =~ /^\s*(\d+)\s+([ATCGatcgNn]+)/) {
			dual_print2 "DNA sequence line\n";
            # If line is a numbered sequence line
            my ($pos_begin,$seqline) = ($1, $2);   # Sequence position
            $current_pos = length($seqline) + $pos_begin - 1;
        }
        elsif ( ($_ =~ /^;+\s+G-(\w.*)/) or ($_ =~ /^;; mfannot:\s+(\/group=.*)/) or ($_ =~ /^;; mfannot:$/) or ($_ =~ /^;+\s+(rnl.*)/) or ($_ =~ /^;+\s+(rns.*)/) ){
			dual_print2 "Feature line\n";
			if ( ($_ =~ /^;+\s+G-(\w.*)/) or ($_ =~ /^;+\s+(rnl.*)/) or ($_ =~ /^;+\s+(rns.*)/) ){

				# If line is a feature boundary, save that information
				my @splitline = split /\s/, $1;
				my $current_name = $splitline[0];
				my $current_direction = $splitline[1];
				my $current_startend = $splitline[2];
				my $type = undef;

				if ($previousIntron){
					$type = "group_II_intron";
					if (defined $startend_hash{$current_contig}{$previousIntron}{$type}{"end"}{0}) {
							my $i = keys %{$startend_hash{$current_contig}{$previousIntron}{$type}{"end"}};
							$startend_hash{$current_contig}{$previousIntron}{$type}{"end"}{$i} = $current_pos;
							dual_print2 "Feature ". $previousIntron. " already defined. Please manually verify in $mfannot_file\n";
					}
					else { 
						$startend_hash{$current_contig}{$previousIntron}{$type}{"end"}{0} = $current_pos; 
					}
					$previousIntron = undef;
				}

				if ($previousRns){
					$type = "rRNA";
					if (defined $startend_hash{$current_contig}{$previousRns}{$type}{"end"}{0}) {
							my $i = keys %{$startend_hash{$current_contig}{$previousRns}{$type}{"end"}};
							$startend_hash{$current_contig}{$previousRns}{$type}{"end"}{$i} = $current_pos;
							dual_print2 "Feature ". $previousRns. " already defined. Please manually verify in $mfannot_file\n";
					}
					else { $startend_hash{$current_contig}{$previousRns}{$type}{"end"}{0} = $current_pos; }
					$previousRns = undef;
				}

				if ($previousRnl){
					$type = "rRNA";
					if (defined $startend_hash{$current_contig}{$previousRnl}{$type}{"end"}{0}) {
							my $i = keys %{$startend_hash{$current_contig}{$previousRnl}{$type}{"end"}};
							$startend_hash{$current_contig}{$previousRnl}{$type}{"end"}{$i} = $current_pos;
							dual_print2 "Feature ". $previousRnl. " already defined. Please manually verify in $mfannot_file\n";
					}
					else { $startend_hash{$current_contig}{$previousRnl}{$type}{"end"}{0} = $current_pos; }
					$previousRnl = undef;
				}

				# --- NOT PREVIOUSLY DEFINED ---

				# --- RNL RNS lines ---
				if( $current_startend eq ";;"){ #rns rnl cases

					if( $current_name eq "rnl"){	
						if ($previousRnl){
							if (defined $startend_hash{$current_contig}{$previousRnl}{"rRNA"}{"end"}{0}) {
									my $i = keys %{$startend_hash{$current_contig}{$previousRnl}{"rRNA"}{"end"}};
									$startend_hash{$current_contig}{$previousRnl}{"rRNA"}{"end"}{$i} = $current_pos;
									dual_print2 "Feature ". $previousRnl. " already defined. Please manually verify in $mfannot_file\n";
							}
							else { $startend_hash{$current_contig}{$previousRnl}{"rRNA"}{"end"}{0} = $current_pos; }
							$previousRnl = undef;
							next;
						}

						if (defined $startend_hash{$current_contig}{$current_name}{"rRNA"}{"start"}{0} ) {
								my $i = keys %{$startend_hash{$current_contig}{$current_name}{"rRNA"}{"start"}};
								$startend_hash{$current_contig}{$current_name}{"rRNA"}{"start"}{$i} = $current_pos + 1;
								dual_print2 "Feature ". $current_name. " already defined. Please manually verify in $mfannot_file\n";
						}
						else { $startend_hash{$current_contig}{$current_name}{"rRNA"}{"start"}{0} = $current_pos + 1;}
						$previousRnl=$current_name;
					}

					if( $current_name eq "rns"){
						if ($previousRns){
							if (defined $startend_hash{$current_contig}{$previousRns}{"rRNA"}{"end"}{0}) {
									my $i = keys %{$startend_hash{$current_contig}{$previousRns}{"rRNA"}{"end"}};
									$startend_hash{$current_contig}{$previousRns}{"rRNA"}{"end"}{$i} = $current_pos;
									dual_print2 "Feature ". $previousRns. " already defined. Please manually verify in $mfannot_file\n";
							}
							else { $startend_hash{$current_contig}{$previousRns}{"rRNA"}{"end"}{0} = $current_pos; }
							$previousRns = undef;
							next;
						}

						if (defined $startend_hash{$current_contig}{$current_name}{"rRNA"}{"start"}{0} ) {
								my $i = keys %{$startend_hash{$current_contig}{$current_name}{"rRNA"}{"start"}};
								$startend_hash{$current_contig}{$current_name}{"rRNA"}{"start"}{$i} = $current_pos + 1;
								dual_print2 "Feature ". $current_name. " already defined. Please manually verify in $mfannot_file\n";
						}
						else { $startend_hash{$current_contig}{$current_name}{"rRNA"}{"start"}{0} = $current_pos + 1;}
						$previousRns=$current_name;
					}

				}

				#  --- GENE lines ----

				elsif ( ( $current_direction eq "<==" or $current_direction eq "==>" )  && ( $current_startend eq "end" or $current_startend eq "start" ) ) {
					
					# Try to better define name and detect type of the feature
					my @split_current_name = split /-/, $current_name;
					# More than one element in the array
					if ($#split_current_name > 0){
						if ($split_current_name[1] =~ /^I/){ $type="intron";}
						elsif ($split_current_name[1] =~ /^E/){ $type="exon";}
						$current_name=$split_current_name[0]; # remove the -I or -E from the genename
					}
					elsif ($current_name =~ /^orf/){ $type="orf";}
					elsif ($current_name =~ /^trn/){ $type="tRNA";}
					else { $type="mRNA";}


					
					if ($current_direction eq "<==" && $current_startend eq "start" ) {

						if (defined $startend_hash{$current_contig}{$current_name}{$type}{"start"}) {
							
							if ($previousDirection eq $current_direction and $previousStartEnd eq $current_startend){ #keep the first key and the second value
								my $i = keys %{$startend_hash{$current_contig}{$current_name}{$type}{"start"}};
								$startend_hash{$current_contig}{$current_name}{$type}{"start"}{$i-1} = $current_pos;
								dual_print2 "11 - Feature ". $current_name. " already defined. Please manually verify in $mfannot_file\n";
								next;
							}

							my $i = keys %{$startend_hash{$current_contig}{$current_name}{$type}{"start"}};
							$startend_hash{$current_contig}{$current_name}{$type}{"start"}{$i} = $current_pos;
							dual_print2 "1 - Feature ". $current_name. " already defined. Please manually verify in $mfannot_file\n";
						}
						else { $startend_hash{$current_contig}{$current_name}{$type}{"start"}{0} = $current_pos; }
					}
					elsif ($current_direction eq "==>" && $current_startend eq "end" ) {

						if (defined $startend_hash{$current_contig}{$current_name}{$type}{"end"}{0}) {

							if ($previousDirection eq $current_direction and $previousStartEnd eq $current_startend){ #keep the first key and the second value
								my $i = keys %{$startend_hash{$current_contig}{$current_name}{$type}{"end"}};
								$startend_hash{$current_contig}{$current_name}{$type}{"end"}{$i-1} = $current_pos;
								dual_print2 "22 - Feature ". $current_name. " already defined. Please manually verify in $mfannot_file\n";
								next;
							}

							my $i = keys %{$startend_hash{$current_contig}{$current_name}{$type}{"end"}};
							$startend_hash{$current_contig}{$current_name}{$type}{"end"}{$i} = $current_pos;
							dual_print2 "2 - Feature ". $current_name. " already defined. Please manually verify in $mfannot_file\n";
						}
						else { $startend_hash{$current_contig}{$current_name}{$type}{"end"}{0} = $current_pos; }

					}
					elsif ($current_direction eq "==>" && $current_startend eq "start") {
						my $value = ($current_pos == 1) ? $current_pos : $current_pos + 1;
						if (defined $startend_hash{$current_contig}{$current_name}{$type}{"start"}{0}) {

							if ($previousDirection eq $current_direction and $previousStartEnd eq $current_startend){
								dual_print2 "3 - Feature ". $current_name. " already defined. Please manually verify in $mfannot_file\n";
								next;
							} #keep the first key and the first value

							my $i = keys %{$startend_hash{$current_contig}{$current_name}{$type}{"start"}};
							$startend_hash{$current_contig}{$current_name}{$type}{"start"}{$i} = $value;
							dual_print2 "3 - Feature ". $current_name. " already defined. Please manually verify in $mfannot_file\n";
						}
						else { $startend_hash{$current_contig}{$current_name}{$type}{"start"}{0} = $value; }
					}
					elsif ($current_direction eq "<==" && $current_startend eq "end") {
						my $value = ($current_pos == 1) ? $current_pos  : $current_pos + 1;
						if (defined $startend_hash{$current_contig}{$current_name}{$type}{"end"}{0}) {

							if ($previousDirection eq $current_direction and $previousStartEnd eq $current_startend){
							dual_print2 "44 - Feature ". $current_name. " already defined. Please manually verify in $mfannot_file\n";
							next;
							} #keep the first key and the first val

							my $i = keys %{$startend_hash{$current_contig}{$current_name}{$type}{"end"}};
							$startend_hash{$current_contig}{$current_name}{$type}{"end"}{$i} = $value;
							dual_print2 "4 - Feature ". $current_name. " already defined. Please manually verify in $mfannot_file\n";
						}
						else { $startend_hash{$current_contig}{$current_name}{$type}{"end"}{0} = $value; }
					}
				}
				# --- COMBINATION UNKNOW ---
				else { 
					print STDERR "Exception to possible combination of feature boundaries and directions: $_ \n"; 
				}

				$previousDirection=$current_direction;
				$previousStartEnd=$current_startend;
			}

			# --- INTRON lines ---
			if ($_ =~ /^;; mfannot:\s+\/(group=.*)/) {
				my $type = "group_II_intron";

				if ($previousIntron){
					if (defined $startend_hash{$current_contig}{$previousIntron}{$type}{"end"}{0}) {
							my $i = keys %{$startend_hash{$current_contig}{$previousIntron}{$type}{"end"}};
							$startend_hash{$current_contig}{$previousIntron}{$type}{"end"}{$i} = $current_pos;
							dual_print2 "Feature ". $previousIntron. " already defined. Please manually verify in $mfannot_file\n";
					}
					else { 
						$startend_hash{$current_contig}{$previousIntron}{$type}{"end"}{0} = $current_pos; 
					}
					$previousIntron = undef;
					next;
				}

				if (defined $startend_hash{$current_contig}{$1}{$type}{"start"}{0} ) {
						my $i = keys %{$startend_hash{$current_contig}{$1}{$type}{"start"}};
						$startend_hash{$current_contig}{$1}{$type}{"start"}{$i} = $current_pos + 1;
						dual_print2 "Feature ".$1. " already defined. Please manually verify in $mfannot_file\n";
				}
				else { 
					$startend_hash{$current_contig}{$1}{$type}{"start"}{0} = $current_pos + 1;
				}
				$previousIntron=$1;
			}
			elsif ($_ =~ /^;; mfannot:$/) {
				my $type = "group_II_intron";
				if ($previousIntron){
					if (defined $startend_hash{$current_contig}{$previousIntron}{$type}{"end"}{0}) {
							my $i = keys %{$startend_hash{$current_contig}{$previousIntron}{$type}{"end"}};
							$startend_hash{$current_contig}{$previousIntron}{$type}{"end"}{$i} = $current_pos;
							dual_print2 "Feature ". $previousIntron. " already defined. Please manually verify in $mfannot_file\n";
					}
					else { 
						$startend_hash{$current_contig}{$previousIntron}{$type}{"end"}{0} = $current_pos; 
					}
					$previousIntron = undef;
					next;
				}
			}
		}
	}
	close(INPUT);
}

# Gene have exon defined and global feature which we define as mRNA
# tRNA have single exon defined and global feature which we define as level2 tRNA
# orf are single exon feature defined as level2 orf
sub handle_records {

	my $gene_uniqid=0;
	foreach my $contig (sort keys %startend_hash){
		foreach my $name (sort keys %{$startend_hash{$contig}} ){

			my $gene_name = $name;
			if ($name =~ /(_[0-9]+)/) {
				my @splitline = split /_/, $name;
				$gene_name = $splitline[0];
			}

			foreach my $type ( sort keys %{$startend_hash{$contig}{$name}} ){

				# for each subfeature e.g. exon 0,1,2
				foreach my $nb ( sort keys %{$startend_hash{$contig}{$name}{$type}{'start'}} ){
					#get the start and end of the feature
					my $start = $startend_hash{$contig}{$name}{$type}{'start'}{$nb};
					my $end = $startend_hash{$contig}{$name}{$type}{'end'}{$nb};
					# get the strand
					my $featuredir = "+";
					if ( $start > $end) {
							$featuredir = "-";
							my $tmpstart=$start;
							$start = $end;
							$end = $tmpstart;
					}
					
					# Shift to exon for allowing multiple exons genes
					my $realType = $type;											
					if ($type eq "rRNA"){
						$realType = "exon";
					}

					# Create unique ID
					my $id = $realType."_".$name;
					if(! exists_keys($hashID,($id)) ){
						$hashID->{$id}++;
					}
					else {
						$id = $id."_".$hashID->{$id};
						$hashID->{$id}++;
					}
					
					# create feature
					my $feature = AGAT::SeqFeatureLite->new(-seq_id => $contig,
																-source_tag => "AGAT",
																-primary_tag => $realType,
																-start => $start,
																-end => $end ,
																-phase => ".",
																-strand => $featuredir,
																-tag => {'ID' => $id, 'Name' => $gene_name, 'locus_tag' => $name }
																) ;
					
					if ($type eq "rRNA"){
						create_or_replace_tag($feature , "agat_parent_type", $type)
					}

					if ($realType eq "mRNA" or $realType eq "tRNA"  or $realType eq "rRNA" or $realType eq "orf" ){
						push (@{$omniscient->{"level2"}{$realType}{lc($name)}}, $feature);
					} else {
						push (@{$omniscient->{"level3"}{$realType}{lc($name)}}, $feature);

					}	
				}
			}
		}
	}
}

__END__

=head1 NAME

agat_convert_mfannot2gff.pl

=head1 DESCRIPTION

Conversion utility for MFannot "masterfile" annotation produced by the MFannot
pipeline (http://megasun.bch.umontreal.ca/RNAweasel/). Reports GFF3 format.

=head1 SYNOPSIS

    agat_convert_mfannot2gff.pl -m <mfannot> -o <gff>
    agat_convert_mfannot2gff.pl --help

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2015, Brandon Seah (kbseah@mpi-bremen.de)
... GPL-3 ...
modified by jacques dainat 2017-11

=head1 OPTIONS

=over 8

=item B<-m> or B<-i> or B<--mfannot> <file>

The mfannot input file

=item B<-o>, B<--out> or B<--output> <file>

Output GFF file to create. If no output file is specified, the output will be
written to STDOUT.

=item B<-h> or B<--help>

Display this helpful text.

=item B<-v> or B<--verbose> <int>

Verbosity, choice are 0,1,2,3,4. 0 is quiet, 1 is normal, 2,3,4 is more verbose. Default 1.

=back

=head1 FEEDBACK

For questions, suggestions, or general discussions about AGAT, please use the AGAT community forum:
https://github.com/NBISweden/AGAT/discussions

=head1 BUG REPORTING

Bug reports should be submitted through the AGAT GitHub issue tracker:
https://github.com/NBISweden/AGAT/issues

=cut

AUTHOR - Jacques Dainat
