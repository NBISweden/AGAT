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

sub usage {
    print STDERR "Convert Mfannot Masterfile to GFF3 format\n";
    print STDERR "\n";
    print STDERR "Usage: perl mfannot2gff.pl -m input.new -g output.gff \n";
    print STDERR "\n";
    exit();
}

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

	foreach my $contig (keys %startend_hash){
		foreach my $name (keys %{$startend_hash{$contig}} ){
			foreach my $nb ( keys %{$startend_hash{$contig}{$name}{'start'}} ){
				my $start = $startend_hash{$contig}{$name}{'start'}{$nb};
				$filtered_result{$contig}{"$start$name"}{$name}{'start'}=$start;
				$filtered_result{$contig}{"$start$name"}{$name}{'end'}=$startend_hash{$contig}{$name}{'end'}{$nb};;
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
  	#use Data::Dumper; #print Dumper(\%filtered_result);exit;
    print GFF "##gff-version 3\n";  # header line

		foreach my $current_contig ( sort keys %filtered_result ){
      foreach my $thefeature_uniq ( sort { (($a =~ /^(\d+)/)[0] || 0) <=> (($b =~ /^(\d+)/)[0] || 0) } keys %{$filtered_result{$current_contig}}) {

				foreach my $thefeature ( keys %{$filtered_result{$current_contig}{$thefeature_uniq}}) {
					my $featuretype;
					#print $thefeature."\n";
          if ($thefeature =~ /^rnl/ | $thefeature =~ /^rns/) { $featuretype="rRNA"; }
          elsif ($thefeature =~ /^trn/) { $featuretype = "tRNA"; }
					elsif ($thefeature =~ /^group/){$featuretype = "group_II_intron";}
          else {$featuretype="CDS";}

          my $featuredir;
          my $frame;

          my $start = $filtered_result{$current_contig}{$thefeature_uniq}{$thefeature}{'start'};
          my $end = $filtered_result{$current_contig}{$thefeature_uniq}{$thefeature}{'end'};
          if ( $start > $end) {
              $featuredir = "-";
							my $tmpstart=$start;
              $start = $end;
              $end = $tmpstart;
          }
					else {
              $featuredir="+";
          }
          if ($featuretype eq "CDS") { $frame="0"; } else { $frame = "."; }

					#create uniq ID
					my $uniqID;
					if ( defined($id_hash{$featuretype}) ){
						$id_hash{$featuretype}++;
						$uniqID = $featuretype."".$id_hash{$featuretype};
					}
					else{
						$uniqID = $featuretype."1";
						$id_hash{$featuretype}=1;
					}


					my @gff3_line = ($current_contig,
                           "mfannot",
                           $featuretype,
                           $start,
                           $end,
                           ".",
                           $featuredir,
                           $frame,
                           "ID=$uniqID;Name=$thefeature;transl_table=$gencode_hash{$current_contig};gene=$thefeature"
                           );
          print GFF join ("\t", @gff3_line)."\n";
				}
      }
		}

    close (GFF);
}

=head1 NAME

gaas_convert_mfannot2gff.pl

=head1 DESCRIPTION

Conversion utility for MFannot "masterfile" annotation produced by the MFannot
pipeline (http://megasun.bch.umontreal.ca/RNAweasel/). Reports GFF3 format.

=head1 SYNOPSIS

    gaas_convert_mfannot2gff.pl -m <mfannot> -o <gff>
    gaas_convert_mfannot2gff.pl --help

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
