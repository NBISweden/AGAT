#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long;
use IO::File;
use Pod::Usage;
use AGAT::AGAT;

my $header = get_agat_header();
my $config;
my $cpu;

my @opt_files;
my $opt_output=undef;
my $opt_plot;
my $opt_breaks;
my $Xpercent=1;
my $opt_help = 0;

my @copyARGV=@ARGV;
if ( !GetOptions( 'f|gff|ref|reffile=s' => \@opt_files,
                  'o|out|output=s'      => \$opt_output,
                  'w|window|b|break|breaks=i'  => \$opt_breaks,
                  'x|p=f'               => \$Xpercent,
                  'plot!'               => \$opt_plot,
                  'c|config=s'               => \$config,
                    'thread|threads|cpu|cpus|core|cores|job|jobs=i' => \$cpu,
                  'h|help!'             => \$opt_help ) )
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

if ( ! ( $#opt_files  >= 0) ) {
    pod2usage( {
           -message => "$header\nMust specify at least 1 parameters:\nReference data gff3 file (--gff)\n",
           -verbose => 0,
           -exitval => 1 } );
}

# --- Manage config ---
initialize_agat({ config_file_in => $config, input => $opt_files[0] });
$CONFIG->{cpu} = $cpu if defined($cpu);

# #######################
# # START Manage Option #
# #######################
my $ostreamReport_file;
if (defined($opt_output) ) {
  if (-d $opt_output){
    print "The output directory choosen already exists. Please geve me another Name.\n";exit();
  }
  else{
    mkdir $opt_output;
  }
  $ostreamReport_file = $opt_output."/report.txt";
}

my $ostreamReport = prepare_fileout($ostreamReport_file);

my $string1 .= "usage: $0 @copyARGV\n\n";

print $ostreamReport $string1;
if($opt_output){print $string1;}


#############################
####### Manage R option #####
#############################

#Choose breaks value:
if(! $opt_breaks){
  $opt_breaks="1000";
}

#############################
####### Manage output #######
#############################
my $outputPDF_prefix;
if (defined($opt_output) ) {
  if (-f $opt_output){
      print "Cannot create a directory with the name $opt_output because a file with this name already exists.\n";exit();
  }
 $outputPDF_prefix=$opt_output."/intronPlot_";
}
else{
  $outputPDF_prefix="intronPlot_";
}

# Check if dependencies for plot are available
if($opt_plot){
  if ( ! may_i_plot() ) {
    $opt_plot = undef;
  }
}

# #####################################
# # END Manage OPTION
# #####################################

#                         #######################
# >>>>>>>>>>>>>>>>>>>>>>>>#        MAIN         #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#                         #######################

#PART 1
###################################
# Read input gff3 files one by one and save value in hash of list

my %introns;
foreach my $file (@opt_files){

  print "Reading ".$file,"\n";
  my $log = create_log_file({input => $file});
	$LOGGING->{'log'} = $log ;
  ######################
  ### Parse GFF input #
  my ($hash_omniscient, $hash_mRNAGeneLink) = slurp_gff3_file_JD({ input => $file });
  ### END Parse GFF input #
  #########################

  #print statistics
	print "Compute statistics\n";
	print_omniscient_statistics ({ input => $hash_omniscient,
																 output => $ostreamReport
															 });

  ######################
  ### Parse GFF input #
  # get nb of each feature in omniscient;
  foreach my $tag_l2 (keys %{$hash_omniscient->{'level2'}}){
    foreach my $id_l1 (keys %{$hash_omniscient->{'level2'}{$tag_l2}}){
      my $one_f2 = $hash_omniscient->{'level2'}{$tag_l2}{$id_l1}[0];

      #######################
      #get feature1 and info
      my $feature_l1=undef;
      my $tag_l1;
      foreach my $tag_level1 (keys %{$hash_omniscient->{'level1'}}){
        if (exists ($hash_omniscient->{'level1'}{$tag_level1}{$id_l1})){
          $feature_l1=$hash_omniscient->{'level1'}{$tag_level1}{$id_l1};
          $tag_l1=$tag_level1;
          last;
        }
      }
      if(! $feature_l1){print "Problem ! We didnt retrieve the level1 feature with id $id_l1\n";exit;}

      #####
      # get all level2
      my $All_l2_single=1;
      my $counterL2_match=-1;
      foreach my $feature_l2 ( @{$hash_omniscient->{'level2'}{$tag_l2}{$id_l1}} ){

        #MATCH CASE - We ahve to count the L2 match features
        if($tag_l2 =~ "match"){
          my @sortedList = sort {$a->start <=> $b->start} @{$hash_omniscient->{'level2'}{$tag_l2}{$id_l1}};
          my $indexLastL2 = $#{$hash_omniscient->{'level2'}{$tag_l2}{$id_l1}};
          $counterL2_match++;

          if($counterL2_match > 0 and $counterL2_match <= $indexLastL2){
            my $intronSize = $sortedList[$counterL2_match]->start - $sortedList[$counterL2_match-1]->end;
            push @{$introns{$tag_l2}}, $intronSize;
          }
        }

        ######
        #get all level3
        my $id_l2=lc($feature_l2->_tag_value('ID'));

        foreach my $tag_l3 ( keys %{$hash_omniscient->{'level3'}} ){

          if(exists_keys($hash_omniscient,('level3',$tag_l3, $id_l2))){

          my $counterL3=-1;
          #Initialize intron to 0 to avoid error during printing results
          my $indexLast = $#{$hash_omniscient->{'level3'}{$tag_l3}{$id_l2}};

          my @sortedList = sort {$a->start <=> $b->start} @{$hash_omniscient->{'level3'}{$tag_l3}{$id_l2}};

            foreach my $feature_l3 ( @sortedList ){

              #count number feature of tag_l3 type
              $counterL3++;

              ################
              #Manage Introns#
              # from the second intron to the last (from index 1 to last index of the table sortedList)
              # We go inside this loop only if we have more than 1 feature.
              if($counterL3 > 0 and $counterL3 <= $indexLast){
                my $intronSize = $sortedList[$counterL3]->start - $sortedList[$counterL3-1]->end;
                push @{$introns{$tag_l3}}, $intronSize;
              }
            }# END FOREACH L3
          }
        }
      }
    }
  }
}

# PART 2

foreach  my $tag (sort keys %introns){
  ###############################
  my $biggest_value=0;
  my $pathIntron="tmp_intron_".$tag.".txt";
  my @sorted_intron = (sort { $a <=> $b } @{$introns{$tag}});
  #########################
  # Write value in tmp files

  # Manage Output
  my $ostreamAED   = IO::File->new();
  $ostreamAED->open( $pathIntron, 'w' ) or
        croak(
          sprintf( "Can not open '%s' for writing %s", $pathIntron, $! )
        );
  foreach  my $value ( @sorted_intron ){
    print $ostreamAED "$value\n";
    if($value > $biggest_value){
      $biggest_value=$value;
    }
  }
  $ostreamAED->close();


  # Part 3
  #########################################
  #Calcul longest after remove X percent  #
  my $lastIndex = $#sorted_intron;
  my $nbValueToRemove = int(($Xpercent*($lastIndex+1))/100);
  my $resu =  $sorted_intron[$lastIndex-$nbValueToRemove];

  my $stringPrint =  "Introns in feature $tag: Removing $Xpercent percent of the highest values ($nbValueToRemove values) gives you $resu bp as the longest intron in $tag.\n";

  print $ostreamReport $stringPrint;
  if($opt_output){print $stringPrint;}


  # Part 4
  #########
  # PLOT  #
  if($opt_plot){
    #chose output file name
    my $outputPDF=$outputPDF_prefix.$tag.".pdf";
    #Choose a main title:
    my $title="Intron distribution in $tag";
    #choose x title
    my $xlab="size bp";

    ## check using R
    my $R = Statistics::R->new() or die "Problem with R : $!\n";

    #calculate the breaks
    my $breaks_ok=int($biggest_value/$opt_breaks);
    # try with shorter breaks
    if ($breaks_ok == 0){ $breaks_ok = int($biggest_value/100)};

    if($breaks_ok){
      #R command
      $R->run(qq`

            listValues1=as.matrix(read.table("$pathIntron", sep="\t", he=F))
            pdf("$outputPDF")
            hist1<-hist(listValues1,breaks=$breaks_ok,main="$title", xlab="$xlab")
            plot(hist1\$mids,hist1\$counts)
            mylims <- par("usr")
            # Add Title second plot
            title(main="$title")`
      );

    }
    # If no breaks but we have values
    elsif($biggest_value){
      $R->run(qq`

            listValues1=as.matrix(read.table("$pathIntron", sep="\t", he=F))
            pdf("$outputPDF")
            hist1<-hist(listValues1,main="$title", xlab="$xlab")
            plot(hist1\$mids,hist1\$counts)
            mylims <- par("usr")
            # Add Title second plot
            title(main="$title")`
      );
      # Close the bridge
      $R->stopR();
    }
  }

  # remove temporary files
  unlink $pathIntron;
}

      #########################
      ######### END ###########
      #########################


#######################################################################################################################
        ####################
         #     methods    #
          ################
           ##############
            ############
             ##########
              ########
               ######
                ####
                 ##



__END__


=head1 NAME

agat_sp_manage_introns.pl

=head1 DESCRIPTION

The script provides information about introns (longest, shortest size mean ...) using the statistic method,
then plot all the intron size values to get an overview of the introns size distribution.
It gives you as well the value of the longest intron after removing X percent(s) of the longest (removing potential biais / false positive).

=head1 SYNOPSIS

    agat_sp_manage_introns.pl --gff infile [--out outFile]
    agat_sp_manage_introns.pl --help

=head1 OPTIONS

=over 8

=item B<--gff>, B<-f>, B<--ref> or B<-reffile>

Input GTF/GFF file. You can use several input files by doing: -f file1 -f file2 -f file3

=item  B<-w>, B<--window>, B<--break>, B<--breaks> or B<-b>

It the number of break used within the histogram plot. By default it's 1000. You can modify the value to get something more or less precise.

=item  B<-x>, B<--p>

Allows to modify the X values to calculate the percentage of the longest introns to remove. By default the value is 1 (We remove 1 percent of the longest).

=item  B<--plot>

Allows to create an histogram in pdf of intron sizes distribution.

=item  B<--out>, B<--output> or B<-o>

Output gff3 file where the gene incriminated will be write.

=item B<-c> or B<--config>

String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any, 
otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

=item B<--help> or B<-h>

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
