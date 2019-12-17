#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use List::MoreUtils qw(uniq);
use Sort::Naturally;
use Bio::Tools::GFF;
use AGAT::Omniscient;

my $header = get_agat_header();
my $outfile = undef;
my $gff1 = undef;
my $gff2 = undef;
my $verbose = undef;
my $opt_help= 0;

my @copyARGV=@ARGV;
if ( !GetOptions(
    "help|h"                  => \$opt_help,
    "gff1=s"                  => \$gff1,
    "gff2=s"                  => \$gff2,
    "v!"                      => \$verbose,
    "output|outfile|out|o=s"  => \$outfile))

{
    pod2usage( { -message => 'Failed to parse command line',
                 -verbose => 1,
                 -exitval => 1 } );
}

# Print Help and exit
if ($opt_help) {
    pod2usage( { -verbose => 99,
                 -exitval => 0,
                 -message => "$header\n" } );
}

if ( ! $gff1 or ! $gff2){
    pod2usage( {
           -message => "$header\nAt least 2 parameters are mandatory:\n".
                       "Input reference gff file1 (--gff1)\n".
                       "Input reference gff file2 (--gff2)\n\n",
           -verbose => 0,
           -exitval => 1 } );
}

######################
# Manage output file #
my $report = IO::File->new();
if ($outfile) {
  open($report, '>', $outfile) or die "Could not open file $outfile $!";
}
else{
  $report->fdopen( fileno(STDOUT), 'w' );
}

                #####################
                #     MAIN          #
                #####################


######################
### Parse GFF input #
my ($omniscient1, $hash_mRNAGeneLink1) = slurp_gff3_file_JD({ input => $gff1
                                                              });
my ($omniscient2, $hash_mRNAGeneLink2) = slurp_gff3_file_JD({ input => $gff2
                                                              });
print ("GFF3 files parsed\n");


my $sortBySeq1 = gather_and_sort_l1_location_by_seq_id_and_strand_chimere($omniscient1);
my $sortBySeq2 = gather_and_sort_l1_location_by_seq_id_and_strand_chimere($omniscient2);
print ("GFF3 files sorted\n");

#get top feature first
my $hash = get_levels_info(); # get from the file
my $top_features = $hash->{'other'}{'level'}{'topfeature'};

# ----- Remove $top_features ------
foreach my $sortBySeq ($sortBySeq1, $sortBySeq2){
  foreach my $locusID1 ( sort keys %{$sortBySeq}){ # tag_l1 = gene or repeat etc...

    # Skip to features
    foreach my $type_top_feature (keys %{$top_features}){
      if (exists_keys( $sortBySeq, ($locusID1, $type_top_feature) ) ){
        delete  $sortBySeq->{$locusID1}{$type_top_feature};
      }
    }
  }
}
# ----- END Remove $top_features ------

# ------------------------------------------------------------------------------
# ---- NOW COMPARE LOCATIONS from file 1 and file 2
# ------------------------------------------------------------------------------
my %overlap_info; # < Will save all results about if it overlap or not
foreach my $locusID (  keys %{$sortBySeq1} ){
  foreach my $chimere_type ( keys %{$sortBySeq1->{$locusID}} ){
    print "\n\n\nFirst round chimere_type location1: <$locusID><$chimere_type>\n" if ($verbose);
    my $one_more_round_signalA = 0;
    my $overlap_A = 0;
    my $overlap_B = 0;
    my $current_flattened_locations = [];
    my $locations1;
    my $type_deeper;
    my $l2_type1;
    my $l1_type;

    my $location1_back =[];

    while ( my $location1 = shift @{$sortBySeq1->{$locusID}{$chimere_type}} ){
      my $l1_id1 = lc($location1->[0]);
      $l1_type = lc($location1->[3]);
      ($locations1, $type_deeper, $l2_type1) = get_locations_deeper($omniscient1, $l1_id1);
      if ($verbose) { print "GENERAL list of location1 <$type_deeper>: "; foreach my $array ( @{$locations1}){print "@{$array} - "; } print "\n";}
      if ($verbose) { print "Lets work with $l2_type1 @$location1\n"};
      if ($verbose and $current_flattened_locations) { print "current_flattened_locations investigated: "; foreach my $array ( @{$current_flattened_locations}){print "@{$array} - "; } print "\n";}

      # First Round
      if (! @$current_flattened_locations){
        $overlap_A=1;
        $current_flattened_locations = $locations1;
        print "First round location1!\n" if ($verbose);
      }

      # NOT First Round
      else{
        my $continue = 1;
        $one_more_round_signalA = undef;
        while ( $continue ){
          #  location A  -------------------------
          #  location A bis                                     -------------------------
          if(  $current_flattened_locations->[$#$current_flattened_locations][1] < $location1->[1] ){ # End Location A < Start Location A bis
            $continue = undef;

            if(! $one_more_round_signalA){
                print "Let save the result 1: ( A $overlap_A => $overlap_B)\n" if ($verbose);
                $overlap_info{$l1_type}{$l2_type1}{$overlap_A}{$overlap_B}++;
                $current_flattened_locations = $locations1;
                $overlap_A=1;
                $overlap_B=0;
                print "Re-initialize location A search\n" if ($verbose);
                if ($verbose and $current_flattened_locations) { print "current_flattened_locations investigated: "; foreach my $array ( @{$current_flattened_locations}){print "@{$array} - "; } print "\n";}
            }
            else{
              print "Push back Location1!!<--------\n" if ($verbose);
              push @{$location1_back}, $location1;
            }

          }
          else{
            my ($overlap, $flattened_locations) = flatten_locations_and_merge($current_flattened_locations, $locations1);
            if( $overlap ){
              print "Location1 still overlap! Lets append current_flattened_locations with location1\n" if ($verbose);
              $current_flattened_locations = $flattened_locations;
              $overlap_A++;
              $one_more_round_signalA=1;

              print "take next location A\n" if ($verbose);

              if(scalar @{$sortBySeq1->{$locusID}{$chimere_type}} != 0){
                $location1 = shift @{$sortBySeq1->{$locusID}{$chimere_type}};
              }
              else{print "Print no more location A\n" if ($verbose); $continue = undef;}
            }
            else{
              print "Push back Location1-2!!<--------\n" if ($verbose);
              print "Get next Location1.\n" if ($verbose);

              push @{$location1_back}, $location1;
              $location1 = shift @{$sortBySeq1->{$locusID}{$chimere_type}};
            }
          }
        }
      }

      if (@$location1_back){
          if ($verbose) { print "location A to put back into the list: ";  foreach my $array ( @{$location1_back}){print "@{$array} - "; } print "\n";}
          unshift @{$sortBySeq1->{$locusID}{$chimere_type}}, @{$location1_back};
          $location1_back =[];
      }

#      if ($verbose) { print "current_flattened_locations investigated Again: "; foreach my $array ( @{$current_flattened_locations}){print "@{$array} - "; } print "\n";}

      if ( exists_keys ($sortBySeq2, ($locusID,$chimere_type) ) ){

        # No more locationB, Store the locations A
        if(scalar @{$sortBySeq2->{$locusID}{$chimere_type}} == 0){
          print "  no more location B\n" if ($verbose);
          $overlap_info{$l1_type}{$l2_type1}{$overlap_A}{$overlap_B}++;
          print "Let save the result 2 ( A $overlap_A => $overlap_B)\n" if ($verbose);
          $overlap_A = 0;
          $overlap_B = 0;
        }

        my $location2_back =[];
        while ( my $location2 = shift @{$sortBySeq2->{$locusID}{$chimere_type}} ){
          my $l1_id2 = lc($location2->[0]);
          print " location2 investigated:  @$location2\n" if ($verbose);


          #  flatenned locations A                       ----  -----  -------
          #  location B             ---------------
          if ($location2->[2] < $current_flattened_locations->[0][0]) {
            print " location2 before location1!\n" if ($verbose);

            if($overlap_A == 1){
              $overlap_info{$l1_type}{$l2_type1}{0}{1}++;
              print "  Let save the result 3 ( A 0 => 1)\n" if ($verbose); # uniq to annotationB
            }
            else{
              print "  Push back Location2-1!!<--------\n" if ($verbose);
              push @{$location2_back}, $location2;
            }
            next;
          }


          # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
          #      ------------ OVERLAP  AT EXON or CDS level -----------
          ## # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
          my ($locations2, $type_deeper2, $l2_type2) = get_locations_deeper($omniscient2, $l1_id2);

          my ($overlap, $flattened_locations) = flatten_locations_and_merge($current_flattened_locations, $locations2);
          if( $overlap ){ #This check is a deeplevel location check !!!!
            print " locations2 overlap current_flattened_locations!\n" if ($verbose);
            if ($verbose and $current_flattened_locations) { print "current_flattened_locations investigated: "; foreach my $array ( @{$current_flattened_locations}){print "@{$array} - "; } print "\n";}
            if ($verbose and $current_flattened_locations) { print "locations2 investigated: "; foreach my $array ( @{$locations2}){print "@{$array} - "; } print "\n";}
            $current_flattened_locations = $flattened_locations;
            $overlap_B++;
          }
          #  location A  -------------------------
          #  location B                                     -------------------------
          elsif( $location2->[1] > $current_flattened_locations->[$#$current_flattened_locations][1]  ){
            print " Location2 after Location1\n" if ($verbose);
            print " Push back Location2-2!!<--------\n" if ($verbose);
            push @{$location2_back}, $location2;
            if ($verbose) { print " location B to put back into the list";  foreach my $array ( @{$location2_back}){print "@{$array} - "; } print "\n";}
            unshift @{$sortBySeq2->{$locusID}{$chimere_type}}, @{$location2_back};
            $location2_back =[];
            last; # Go back to the list of locationA
          }
          #  location deeplevel A  ----------             --------                   (deep level means we don't look at top feature but rather at CDS / exons level)
          #  location deeplevel B                -------             -----------
          else{
              push @{$location2_back}, $location2;
          }
        }
      }
    }

    print "no more location A\n" if ($verbose);
    if($overlap_A){
      print "Let save the result 4 (A $overlap_A => $overlap_B)\n" if ($verbose);
      $overlap_info{$l1_type}{$l2_type1}{$overlap_A}{$overlap_B}++;
      $overlap_A = 0;
      $overlap_B = 0;
    }
  }
}

# ---- NOw deal with what is remaining in annotationB
 foreach my $locusID (  keys %{$sortBySeq2} ){
   foreach my $chimere_type ( keys %{$sortBySeq2->{$locusID}} ){
     while ( my $location2 = shift @{$sortBySeq2->{$locusID}{$chimere_type}} ){
       my $l1_id1 = lc($location2->[0]);
       my $l1_type = lc($location2->[3]);
       my ($tothrow1, $tothrow2, $l2_type1) = get_locations_deeper($omniscient2, $l1_id1);
       $overlap_info{$l1_type}{$l2_type1}{0}{1}++; # uniq to annotationB
       print "Let save the result 5 ( A 0 => 1)\n" if ($verbose);
     }
   }
 }



##############
# STATISTICS #
if($verbose){
  print "stat file1:\n";
  my ($stat, $distri) = gff3_statistics($omniscient1);
  foreach my $infoList (@$stat){
    foreach my $info (@$infoList){
      print "$info";
    }
    print "\n";
  }
  print "stat file2:\n";
  ($stat, $distri) = gff3_statistics($omniscient2);
  foreach my $infoList (@$stat){
    foreach my $info (@$infoList){
      print "$info";
    }
    print "\n";
  }
}

# ------------------------------------------------------------------------------
# ------------------------- Now print the Results -------------------------
# ------------------------------------------------------------------------------
my %total;
my $separator_table = join('', '-') x 94;
$separator_table .= "\n";
my ($filename1) = fileparse($gff1,qr/\.[^.]*/);
my ($filename2) = fileparse($gff2,qr/\.[^.]*/);
my $string_to_print = "usage: $0 @copyARGV\nResults of number of genes from file1 that overlap genes from file2:\n\n";


foreach my $type_l1 ( sort keys %overlap_info ){
  foreach my $type_l2 ( sort keys %{$overlap_info{$type_l1}} ){

    $string_to_print .= "$separator_table|".sizedPrint("$type_l1 with $type_l2",92)."|\n";
    $string_to_print .= "$separator_table|".sizedPrint($filename1,30)."|".sizedPrint($filename2,30)."|".sizedPrint("Number of cases",30)."|\n$separator_table";

    foreach my $value1 ( sort {$a <=> $b} keys %{$overlap_info{$type_l1}{$type_l2}} ){
      foreach my $value2 ( sort {$a <=> $b}  keys %{$overlap_info{$type_l1}{$type_l2}{$value1}} ){
        $string_to_print .= "|".sizedPrint($value1,30)."|".sizedPrint($value2,30)."|".sizedPrint($overlap_info{$type_l1}{$type_l2}{$value1}{$value2},30)."|\n";
        if ($value1 != 0){
          $total{$type_l1}{$type_l2}{'A'} += $value1 * $overlap_info{$type_l1}{$type_l2}{$value1}{$value2};
        }
        if ($value2 != 0){
          $total{$type_l1}{$type_l2}{'B'} += $value2 * $overlap_info{$type_l1}{$type_l2}{$value1}{$value2};
        }
      }
    }
    $string_to_print .=  $separator_table;
    $string_to_print .= "Number gene in $filename1: $total{$type_l1}{$type_l2}{'A'}\nNumber gene in $filename2: $total{$type_l1}{$type_l2}{'B'}\n";
    $string_to_print .= "\n\n"
  }
}
$string_to_print .= "\n";

if ($outfile){
  print $report $string_to_print;
}
print $string_to_print;
print "Bye Bye.\n";
#######################################################################################################################
        ####################
         #     METHODS    #
          ################
           ##############
            ############
             ##########
              ########
               ######
                ####
                 ##

# Get location from level3 or level2 is no level3
sub get_locations_deeper{
  my ($omniscient, $l1_id, $verbose) = @_;
  print "get_locations_deeper...\n" if ($verbose);
  my @locations;
  my @locations_l2;
  my $type;

  foreach my $type_l2 (keys %{$omniscient->{'level2'}}){
    print "type_l2 $type_l2\n" if ($verbose);
    if(exists_keys($omniscient,('level2', $type_l2, $l1_id) ) ){
      print "exist l2" if ($verbose);
      foreach my $l2_f ( @{$omniscient->{'level2'}{$type_l2}{$l1_id} } ){
        # Define location l2
        my $l2_id = lc($l2_f->_tag_value('ID'));
        push @locations_l2,[ int($l2_f->start()), int($l2_f->end()) ];
        ################
        # Go to level3 #
        ################
        if(exists_keys($omniscient,('level3', 'cds', $l2_id) ) ){
          print "exits cds" if ($verbose);
          foreach my $l3_f ( @{$omniscient->{'level3'}{'cds'}{$l2_id} } ){
            # Define location l3
            my $current_location_l3 = [int($l3_f->start()), int($l3_f->end()) ];
            push @locations, $current_location_l3 ;
            $type = 'cds';
          }
        }
        if(! @locations){
          if(exists_keys($omniscient,('level3', 'exon', $l2_id) ) ){
            print "exits exon" if ($verbose);
            foreach my $l3_f ( @{$omniscient->{'level3'}{'exon'}{$l2_id} } ){
              # Define location l3
              my $current_location_l3 = [int($l3_f->start()), int($l3_f->end()) ];
              push @locations, $current_location_l3 ;
              $type = 'exon';
            }
          }
        }
      }
      # No L3 location
      if(! @locations){
        @locations = @locations_l2;
        $type = $type_l2;
      }
      my ($tothrow, $locations_fixed) = flatten_locations_and_merge(\@locations, \@locations); # merge (flatten) locations coming from isoforms
      return $locations_fixed, $type, $type_l2;
    }
  }
}

# merge overlapping locations in a list of locations
# lcation [[int, int],[int, int],[int, int]]
sub flatten_locations_and_merge{
  my ($locations1, $locations2) = @_;

  my $locations = [@$locations1, @$locations2];
  my @newlocations;
  my $overlap =undef;
  my $previous_location = undef;
  foreach my $location ( sort {$a->[0] <=> $b->[0]} @{$locations} ){

    # first round
    if (! $previous_location){
       push @newlocations, $location;
       $previous_location = $location;
    }
    # Not first round
    else{
      #  OVERLAP
      if ( ($previous_location->[0] <= $location->[1]) and ($previous_location->[1] >= $location->[0])){
        $overlap=1;
        #  location A  -------------------------
        #  location B           -------------------------
        if($previous_location->[1] <= $location->[1]){
          # take back last location pushed in the array
          $previous_location = pop @newlocations;
          $previous_location  = [$previous_location->[0], $location->[1]];
          # push back into the array the previous location that has been modified
          push @newlocations, $previous_location ;
        }
      }
      else{
        push @newlocations, $location ;
        $previous_location = $location;
      }
    }
  }
  return $overlap, \@newlocations ;
}

 # @Purpose: Create a hash of level1 location (location = [level1ID,start,end]) sorted by feature type and localisation. A localisation is the sequence_id appended by the strand
 # @input: 1 => hash omniscient
 # @output: 1 => hash => LocusID->typeFeatureChimere =[ID,start,end, type]
 # TagChimere allows to divide the L1 into relatedted l2 type. (e.g like do split by level 2 feature)
 sub gather_and_sort_l1_location_by_seq_id_and_strand_chimere{
 	my ($omniscient) = @_;

 	my %hash_sortBySeq;

 	foreach my $tag_level1 (keys %{$omniscient->{'level1'}} ){
   	foreach my $level1_id (keys %{$omniscient->{'level1'}{$tag_level1}} ){
	    my $level1_feature = $omniscient->{'level1'}{$tag_level1}{$level1_id};

      my $tag_chimere = $tag_level1;
      my $l2_type_keep=undef;
      foreach my $l2_type (keys %{$omniscient->{'level2'}} ){
        if( exists_keys($omniscient,('level2', $l2_type, $level1_id ) ) ) {
          $l2_type_keep = $l2_type;
        }
      }
      if($l2_type_keep){$tag_chimere .= $l2_type_keep;}

    	my $ID = $level1_feature->_tag_value('ID');
	    my $strand="+";
	    if($level1_feature->strand != 1){$strand = "-";}
	    my $position_l1=$level1_feature->seq_id."".$strand;
	    push ( @{$hash_sortBySeq{$position_l1}{$tag_chimere}}, [$ID, int($level1_feature->start), int($level1_feature->end), $tag_level1] );
    }

    foreach my $position_l1 (keys %hash_sortBySeq){
      foreach my $tag_chimere (keys %{$hash_sortBySeq{$position_l1}} ){
        @{$hash_sortBySeq{$position_l1}{$tag_chimere}} = sort { ncmp ( $a->[1], $b->[1] ) } @{$hash_sortBySeq{$position_l1}{$tag_chimere}};
      }
    }
 	}
 	return \%hash_sortBySeq;
}

__END__

=head1 NAME

agat_sp_compare_two_annotations.pl

=head1 DESCRIPTION

The script aims to compare two annotation of the same assembly. It provided
information about split/fusion of genes between the two annotations.
The most common case are:
1 => 0 ( gene uniq to file1)
0 => 1 ( gene uniq to file2)
1 => 1 ( 1 gene from file 1 overlaps only 1 gene from file2)
1 => <many> ( 1 gene from file 1 overlaps <many> genes from file2) => split case (with file 1 as reference)
<many> => 1 ( <many> genes from file 1 overlap only 1 gene from file2) => fusion case (with file 1 as reference)

Then you can get more complex cases:
<many> => <many>  (<many> genes from file 1 overlap <many> genes from file2)

=head1 SYNOPSIS

    agat_sp_compare_two_annotations.pl -gff1 infile.gff [ -o outfile ]
    agat_sp_compare_two_annotations.pl --help

=head1 OPTIONS

=over 8

=item B<-gff1>

Input GTF/GFF file1.

=item B<-gff2>

Input GTF/GFF file2.

=item B<-o> , B<--output> , B<--out> or B<--outfile>

Output GFF file.  If no output file is specified, the output will be
written to STDOUT.

=item B<-v>

Verbose option, make it easier to follow what is going on for debugging purpose.

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

AUTHOR - Jacques Dainat
