#!/usr/bin/env perl

# Convert Mfannot output file to GFF3 format
# kbseah@mpi-bremen.de      2015-04-01
# modified by jacques dainat: jacques.dainat@nbis.se

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use AGAT::Omniscient;

my $header = get_agat_header();
my $mfannot_file;
my $verbose;
my $gff_file;
my %startend_hash;     # Stores start and end positions of each feature reported
my %sorted_hash;
my %hash_uniqID;
my %filtered_result;
my %gencode_hash;

GetOptions(
    'mfannot|m|i=s' => \$mfannot_file,
    'gff|g|o=s' => \$gff_file,
		'v|verbose!' => \$verbose,
    'help|h' => sub { pod2usage( -exitstatus=>0, -verbose=>99, -message => "$header\n" ); },
    'man' => sub { pod2usage(-exitstatus=>0, -verbose=>2); }
) or pod2usage ( -exitstatus=>2, -verbose=>2 );

if (!defined $mfannot_file) {
    pod2usage( -message=>"Insufficient options supplied", -exitstatus=>2 );
}


## MAIN ##############################################################

read_mfannot($mfannot_file);
sort_result();
write_gff($gff_file);

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

        if ($_ =~ /^>(.*) gc=(\d+)/) {
            # If a header line, update the current contig and genetic code
            ($current_contig, $current_genetic_code) = ($1, $2);
            $current_pos=1; # Reset the position counter
            $gencode_hash{$current_contig} = $current_genetic_code;
        }
        elsif ($_ =~ /^\s*(\d+)\s+([ATCGatcgNn]+)/) {
            # If line is a numbered sequence line
            my ($pos_begin,$seqline) = ($1, $2);   # Sequence position
            $current_pos = length($seqline) + $pos_begin - 1;
        }
        elsif ( ($_ =~ /^;+\s+G-(\w.*)/) or ($_ =~ /^;; mfannot:\s+(\/group=.*)/) or ($_ =~ /^;+\s+(rnl.*)/) or ($_ =~ /^;+\s+(rns.*)/) ){

					if ( ($_ =~ /^;+\s+G-(\w.*)/) or ($_ =~ /^;+\s+(rnl.*)/) or ($_ =~ /^;+\s+(rns.*)/) ){

						# If line is a feature boundary, save that information
            my @splitline = split /\s/, $1;
						my $current_name = $splitline[0];
						my $current_direction = $splitline[1];
						my $current_startend = $splitline[2];

						if ($previousIntron){
							if (defined $startend_hash{$current_contig}{$previousIntron}{"end"}{0}) {
									my $i = keys %{$startend_hash{$current_contig}{$previousIntron}{"end"}};
									$startend_hash{$current_contig}{$previousIntron}{"end"}{$i} = $current_pos;
									print "Feature ". $previousIntron. " already defined. Please manually verify in $mfannot_file\n" if ($verbose);
							}
							else { $startend_hash{$current_contig}{$previousIntron}{"end"}{0} = $current_pos; }
							$previousIntron = undef;
						}

						if ($previousRns){
							if (defined $startend_hash{$current_contig}{$previousRns}{"end"}{0}) {
									my $i = keys %{$startend_hash{$current_contig}{$previousRns}{"end"}};
									$startend_hash{$current_contig}{$previousRns}{"end"}{$i} = $current_pos;
									print "Feature ". $previousRns. " already defined. Please manually verify in $mfannot_file\n" if ($verbose);
							}
							else { $startend_hash{$current_contig}{$previousRns}{"end"}{0} = $current_pos; }
							$previousRns = undef;
						}

						if ($previousRnl){
							if (defined $startend_hash{$current_contig}{$previousRnl}{"end"}{0}) {
									my $i = keys %{$startend_hash{$current_contig}{$previousRnl}{"end"}};
									$startend_hash{$current_contig}{$previousRnl}{"end"}{$i} = $current_pos;
									print "Feature ". $previousRnl. " already defined. Please manually verify in $mfannot_file\n" if ($verbose);
							}
							else { $startend_hash{$current_contig}{$previousRnl}{"end"}{0} = $current_pos; }
							$previousRnl = undef;
						}

						# gene lines
            if ($current_direction eq "<==" && $current_startend eq "start" ) {
                if (defined $startend_hash{$current_contig}{$current_name}{"start"}) {

                    if ($previousDirection eq $current_direction and $previousStartEnd eq $current_startend){ #keep the first key and the second value
                        my $i = keys %{$startend_hash{$current_contig}{$current_name}{"start"}};
                        $startend_hash{$current_contig}{$current_name}{"start"}{$i-1} = $current_pos;
                        print "11 - Feature ". $current_name. " already defined. Please manually verify in $mfannot_file\n" if ($verbose);
                        next;
                    }

                    my $i = keys %{$startend_hash{$current_contig}{$current_name}{"start"}};
                    $startend_hash{$current_contig}{$current_name}{"start"}{$i} = $current_pos;
                    print "1 - Feature ". $current_name. " already defined. Please manually verify in $mfannot_file\n" if ($verbose);
                }
                else { $startend_hash{$current_contig}{$current_name}{"start"}{0} = $current_pos; }
            }
            elsif ($current_direction eq "==>" && $current_startend eq "end" ) {
                if (defined $startend_hash{$current_contig}{$current_name}{"end"}{0}) {

                    if ($previousDirection eq $current_direction and $previousStartEnd eq $current_startend){ #keep the first key and the second value
                        my $i = keys %{$startend_hash{$current_contig}{$current_name}{"end"}};
                         $startend_hash{$current_contig}{$current_name}{"end"}{$i-1} = $current_pos;
                         print "22 - Feature ". $current_name. " already defined. Please manually verify in $mfannot_file\n" if ($verbose);
                         next;
                    }

                    my $i = keys %{$startend_hash{$current_contig}{$current_name}{"end"}};
                    $startend_hash{$current_contig}{$current_name}{"end"}{$i} = $current_pos;
                    print "2 - Feature ". $current_name. " already defined. Please manually verify in $mfannot_file\n" if ($verbose);
                }
                else { $startend_hash{$current_contig}{$current_name}{"end"}{0} = $current_pos; }

            }
            elsif ($current_direction eq "==>" && $current_startend eq "start") {
                if (defined $startend_hash{$current_contig}{$current_name}{"start"}{0}) {

                    if ($previousDirection eq $current_direction and $previousStartEnd eq $current_startend){
                        print "3 - Feature ". $current_name. " already defined. Please manually verify in $mfannot_file\n" if ($verbose);
                        next;
                    } #keep the first key and the first value

                    my $i = keys %{$startend_hash{$current_contig}{$current_name}{"start"}};
                    $startend_hash{$current_contig}{$current_name}{"start"}{$i} = $current_pos + 1;
                    print "3 - Feature ". $current_name. " already defined. Please manually verify in $mfannot_file\n" if ($verbose);
                }
                else { $startend_hash{$current_contig}{$current_name}{"start"}{0} = $current_pos + 1; }
            }
            elsif ($current_direction eq "<==" && $current_startend eq "end") {
                if (defined $startend_hash{$current_contig}{$current_name}{"end"}{0}) {

                    if ($previousDirection eq $current_direction and $previousStartEnd eq $current_startend){
                    print "44 - Feature ". $current_name. " already defined. Please manually verify in $mfannot_file\n" if ($verbose);
                    next;
                    } #keep the first key and the first val

                    my $i = keys %{$startend_hash{$current_contig}{$current_name}{"end"}};
                    $startend_hash{$current_contig}{$current_name}{"end"}{$i} = $current_pos + 1;
                    print "4 - Feature ". $current_name. " already defined. Please manually verify in $mfannot_file\n" if ($verbose);
                }
                else { $startend_hash{$current_contig}{$current_name}{"end"}{0} = $current_pos + 1; }
            }

						# rnl rns lines
						elsif( $current_startend eq ";;"){ #rns rnl cases
							if( $current_name eq "rnl"){
								if ($previousRnl){
									if (defined $startend_hash{$current_contig}{$previousRnl}{"end"}{0}) {
											my $i = keys %{$startend_hash{$current_contig}{$previousRnl}{"end"}};
											$startend_hash{$current_contig}{$previousRnl}{"end"}{$i} = $current_pos;
											print "Feature ". $previousRnl. " already defined. Please manually verify in $mfannot_file\n" if ($verbose);
									}
									else { $startend_hash{$current_contig}{$previousRnl}{"end"}{0} = $current_pos; }
									$previousRnl = undef;
									next;
								}

								if (defined $startend_hash{$current_contig}{$current_name}{"start"}{0} ) {
										my $i = keys %{$startend_hash{$current_contig}{$current_name}{"start"}};
										$startend_hash{$current_contig}{$current_name}{"start"}{$i} = $current_pos + 1;
										print "Feature ". $current_name. " already defined. Please manually verify in $mfannot_file\n" if ($verbose);
								}
								else { $startend_hash{$current_contig}{$current_name}{"start"}{0} = $current_pos + 1;}
								$previousRnl=$current_name;
							}

							if( $current_name eq "rns"){
								if ($previousRns){
									if (defined $startend_hash{$current_contig}{$previousRns}{"end"}{0}) {
											my $i = keys %{$startend_hash{$current_contig}{$previousRns}{"end"}};
											$startend_hash{$current_contig}{$previousRns}{"end"}{$i} = $current_pos;
											print "Feature ". $previousRns. " already defined. Please manually verify in $mfannot_file\n" if ($verbose);
									}
									else { $startend_hash{$current_contig}{$previousRns}{"end"}{0} = $current_pos; }
									$previousRns = undef;
									next;
								}

								if (defined $startend_hash{$current_contig}{$current_name}{"start"}{0} ) {
										my $i = keys %{$startend_hash{$current_contig}{$current_name}{"start"}};
										$startend_hash{$current_contig}{$current_name}{"start"}{$i} = $current_pos + 1;
										print "Feature ". $current_name. " already defined. Please manually verify in $mfannot_file\n" if ($verbose);
								}
								else { $startend_hash{$current_contig}{$current_name}{"start"}{0} = $current_pos + 1;}
								$previousRns=$current_name;
							}

						}
            else { print STDERR "Exception to possible combination of feature boundaries and directions: $_ \n"; }
            $previousDirection=$current_direction;
            $previousStartEnd=$current_startend;
					}

					# intron lines
					if ($_ =~ /^;; mfannot:\s+\/(group=.*)/) {

						if ($previousIntron){
							if (defined $startend_hash{$current_contig}{$previousIntron}{"end"}{0}) {
									my $i = keys %{$startend_hash{$current_contig}{$previousIntron}{"end"}};
									$startend_hash{$current_contig}{$previousIntron}{"end"}{$i} = $current_pos;
									print "Feature ". $previousIntron. " already defined. Please manually verify in $mfannot_file\n" if ($verbose);
							}
							else { $startend_hash{$current_contig}{$previousIntron}{"end"}{0} = $current_pos; }
							$previousIntron = undef;
							next;
						}

						if (defined $startend_hash{$current_contig}{$1}{"start"}{0} ) {
								my $i = keys %{$startend_hash{$current_contig}{$1}{"start"}};
								$startend_hash{$current_contig}{$1}{"start"}{$i} = $current_pos + 1;
								print "Feature ". $1. " already defined. Please manually verify in $mfannot_file\n" if ($verbose);
						}
						else { $startend_hash{$current_contig}{$1}{"start"}{0} = $current_pos + 1;}
						$previousIntron=$1;
	        }
				}
    }
    close(INPUT);
}

sub sort_result {

	my $gene_uniqid=0;
	foreach my $contig (keys %startend_hash){
		foreach my $name (keys %{$startend_hash{$contig}} ){
			foreach my $nb ( keys %{$startend_hash{$contig}{$name}{'start'}} ){
				my $start = $startend_hash{$contig}{$name}{'start'}{$nb};
				my $end = $startend_hash{$contig}{$name}{'end'}{$nb};
				my $featuredir = "+";
				if ( $start > $end) {
						$featuredir = "-";
						my $tmpstart=$start;
						$start = $end;
						$end = $tmpstart;
				}

				my $parent = undef;
				my $type = undef;
				my $gene_id = undef;
				my $gene_name = $name;
				if ($name =~ /^rnl/ | $name =~ /^rns/) { $type="rRNA"; }
				elsif ($name =~ /^trn/) { $type = "tRNA"; }
				elsif ($name =~ /^group/){$type = "group_II_intron";}
				elsif ($name =~ /^(\w+)-I\w+/){$type="intron"; $parent=$1; $gene_name=$1; $gene_id = $1;}
				elsif ($name =~ /^(\w+)-E\w+/){$type="exon"; $parent=$1; $gene_name=$1; $gene_id = $1;}
				else {$type="mRNA"; $gene_name=$name; $gene_id = $1;}

				if (! $gene_id ){$gene_id = $gene_uniqid++;}

				my %hash_value = (
					start => $start,
					end => $end,
					strand => $featuredir,
					type  => $type,
					parent => $parent,
					name => $name,
					gene_name => $gene_name
				);



				push ( @{$filtered_result{ $contig }{ $gene_id }{ lc($type) }}, {%hash_value} );

				if ($type ne "intron" and $type ne "exon"){
					$sorted_hash{$contig}{"$start$end$name"} =  $gene_id; # to print the features sorted
				}
			}
		}
	}
}

sub write_gff {
		my %id_hash;
		if ($_[0]){
	    open(GFF, ">", "$_[0]") or die ("$!\n");
		}
		else{ # print to STDOUT
			*GFF = *STDOUT;
		}
  	#use Data::Dumper; print Dumper(\%filtered_result);exit;
    print GFF "##gff-version 3\n";  # header line

		foreach my $current_contig ( sort keys %filtered_result ){
      foreach my $uniqid ( sort { (($a =~ /^(\d+)/)[0] || 0) <=> (($b =~ /^(\d+)/)[0] || 0) } keys %{$sorted_hash{$current_contig}}) {
				my $gene_name = $sorted_hash{$current_contig}{$uniqid};

					# mRNA can have exon or not (If none we create one)
					if (exists_keys (\%filtered_result, ($current_contig, $gene_name, 'mrna')) ){
						write_feature($current_contig, $filtered_result{$current_contig}{$gene_name}{'mrna'});
						if (exists_keys (\%filtered_result, ($current_contig, $gene_name, 'exon')) ){
							write_feature($current_contig, $filtered_result{$current_contig}{$gene_name}{'exon'});
						}
						# create exon because none exists
						else{
							my $mrna_hash = $filtered_result{$current_contig}{$gene_name}{'mrna'}[0];

							my %hash_value = (
								start => $mrna_hash->{'start'},
								end => $mrna_hash->{'end'},
								strand => $mrna_hash->{'strand'},
								type  => 'exon',
								parent => $mrna_hash->{'name'},
								name => $mrna_hash->{'name'},
								gene_name => $mrna_hash->{'gene_name'}
							);
							write_feature($current_contig, [\%hash_value] );
						}
						if (exists_keys (\%filtered_result, ($current_contig, $gene_name, 'intron')) ){
							write_feature($current_contig, $filtered_result{$current_contig}{$gene_name}{'intron'});
						}
					}
					# Other than mRNA
					else{
						foreach my $type ( keys %{$filtered_result{$current_contig}{$gene_name}}) {

							write_feature($current_contig, $filtered_result{$current_contig}{$gene_name}{$type});

							#create exon for other feature than group_ii_intron
							if ( lc($type) ne "group_ii_intron" ) {
								my @list_hashes;
								foreach my $other_hash ( @{$filtered_result{$current_contig}{$gene_name}{$type} }){

									my %hash_value = (
										start => $other_hash->{'start'},
										end => $other_hash->{'end'},
										strand => $other_hash->{'strand'},
										type  => 'exon',
										parent => $other_hash->{'id'},
										name => $other_hash->{'name'},
										gene_name => $other_hash->{'gene_name'}
									);
									push @list_hashes, {%hash_value};
								}

								write_feature($current_contig, \@list_hashes );
							}
						}
					}
      }
		}

    close (GFF);
}


sub write_feature{
	my ($contig, $list)=@_;

	foreach my $hash ( sort {$a->{'start'} <=> $b->{'start'}} @{$list} ) {

		# deal with frame
		my $frame;
		if ($hash->{'type'} eq "CDS") { $frame="0"; }
		else { $frame = "."; }

		#ID and Parent
		my $uniqID = create_uniq_id($hash->{'type'});
		$hash->{'id'} = $uniqID;
		my $mandatory = undef;
		if( defined ($hash->{'parent'} ) ){
			$mandatory = "ID=$uniqID;Parent=$hash->{'parent'}";
		}
		else{
			$mandatory = "ID=$uniqID";
		}

		my @gff3_line = ($contig,
										 "mfannot",
										 $hash->{'type'},
										 $hash->{'start'},
										 $hash->{'end'},
										 ".",
										 $hash->{'strand'},
										 $frame,
										 "$mandatory;Name=$hash->{'name'};transl_table=$gencode_hash{$contig};gene=$hash->{'gene_name'}"
										 );
		print GFF join ("\t", @gff3_line)."\n";
	}
}


sub create_uniq_id{
	my ( $tag ) = @_;

	my $uniqID;

	if(! exists_keys(\%hash_uniqID,($tag) ) ){
		$uniqID=$tag."_1";
		$hash_uniqID{$tag}=1;
	}
	else{
		$hash_uniqID{$tag}++;
		$uniqID=$tag."_".$hash_uniqID{$tag};
	}

	return $uniqID;
}


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

=item B<-m> or B<-i> or B<--mfannot>

The mfannot input file

=item B<-g> or B<-o> or B<--gff>

the gff output file

=item B<-h> or B<--help>

Display this helpful text.

=back

=head1 FEEDBACK

=head2 Did you find a bug?

Do not hesitate to report bugs to help us keep track of the bugs and their
resolution. Please use the GitHub issue tracking system available at this
address:

            https://github.com/NBISweden/AGAT/issues

 Ensure that the bug was not already reported by searching under Issues.
 If you're unable to find an (open) issue addressing the problem, open a new one.
 Try as much as possible to include in the issue when relevant:
 - a clear description,
 - as much relevant information as possible,
 - the command used,
 - a data sample,
 - an explanation of the expected behaviour that is not occurring.

=head2 Do you want to contribute?

You are very welcome, visit this address for the Contributing guidelines:
https://github.com/NBISweden/AGAT/blob/master/CONTRIBUTING.md

=cut
