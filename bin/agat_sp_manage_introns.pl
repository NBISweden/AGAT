#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long;
use IO::File;
use Pod::Usage;
use AGAT::AGAT;

my $header = get_agat_header();
start_script();
# ---------------------------- OPTIONS ----------------------------
my @opt_files;
my $opt_output=undef;
my $opt_plot;
my $opt_breaks;
my $Xpercent=1;
my $opt_help = 0;

#############################
# >>>>>>>>>>>>> OPTIONS <<<<<<<<<<<<
#############################
my ($shared_argv, $script_argv) = split_argv_shared_vs_script(\@ARGV);

my $script_parser = Getopt::Long::Parser->new();
if ( ! $script_parser->getoptionsfromarray(
    $script_argv,
    'f|gff|ref|reffile=s' => \@opt_files,
    'o|out|output=s'      => \$opt_output,
    'w|window|b|break|breaks=i'  => \$opt_breaks,
    'x|p=f'               => \$Xpercent,
    'plot!'               => \$opt_plot,
    'h|help!'             => \$opt_help,
  ) )
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

#############################
# >>>>>>> Manage config <<<<<<<
#############################
my $shared_opts = parse_shared_options($shared_argv);
initialize_agat({ config_file_in => ( $shared_opts->{config} ), input => $opt_files[0], shared_opts => $shared_opts });

# #######################
# # START Manage Option #
# #######################
my $ostreamReport_file;
if (defined($opt_output) ) {
  if (-d $opt_output){
    die "The output directory choosen already exists. Please geve me another Name.\n";
  }
  else{
    mkdir $opt_output;
  }
  $ostreamReport_file = $opt_output."/report.txt";
}

my $ostreamReport = prepare_fileout($ostreamReport_file);

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
      die "Cannot create a directory with the name $opt_output because a file with this name already exists.\n";
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

  dual_print1 "Reading ".$file."\n";
  my $log = create_log_file({input => $file});
	$LOGGING->{'log'} = $log ;
  ######################
  ### Parse GFF input #
  my ($hash_omniscient) = slurp_gff3_file_JD({ input => $file });
  ### END Parse GFF input #
  #########################

  #print statistics
  dual_print1 "Compute statistics\n";
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
      if(! $feature_l1){ die "Problem ! We didnt retrieve the level1 feature with id $id_l1\n"; }

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
  if($opt_output){ dual_print1 "$stringPrint"; }


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

# --- final messages ---
end_script( $ostreamReport );

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

=item B<--gff>, B<-f>, B<--ref> or B<-reffile> <file>

Input GTF/GFF file. You can use several input files by doing: -f file1 -f file2 -f file3

=item  B<-w>, B<--window>, B<--break>, B<--breaks> or B<-b> <int>

It the number of break used within the histogram plot. By default it's 1000. You can modify the value to get something more or less precise.

=item  B<-x>, B<--p> <string> <int>

Allows to modify the X values to calculate the percentage of the longest introns to remove. By default the value is 1 (We remove 1 percent of the longest).

=item  B<--plot>

Allows to create an histogram in pdf of intron sizes distribution.

=item  B<-o>, B<--out> or B<--output> <file>

Output file to create (default GFF3 - see config to modify output format).
If no output file is specified, the output will be written to STDOUT.

=item B<--help> or B<-h>

Display this helpful text.

=back

=head1 SHARED OPTIONS

Shared options are defined in the AGAT configuration file and can be overridden via the command line for this script only.
Common shared options are listed below; for the full list, please refer to the AGAT agat_config.yaml.

=over 8

=item B<--config> <file>

Path to a custom AGAT configuration file.  
By default, AGAT uses `agat_config.yaml` from the working directory if present, otherwise the default file shipped with AGAT
(available locally via `agat config --expose`).

=item B<--cpu>, B<--core>, B<--job> or B<--thread> <int>

Number of parallel processes to use for file input parsing (via forking).

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
