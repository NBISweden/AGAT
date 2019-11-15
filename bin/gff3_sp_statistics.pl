#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Statistics::R;
use IO::File;
use List::MoreUtils qw(uniq);
use Bio::Tools::GFF;
use Agat::Omniscient;
use Agat::PlotR;

my $header = get_agat_header();
my $gff = undef;
my $opt_output = undef;
my $opt_genomeSize = undef;
my $opt_plot = undef;
my $help= 0;

if ( !GetOptions(
    "help|h" => \$help,
    'o|output=s'      => \$opt_output,
    'd|p'   => \$opt_plot,
    'g|gs=s' => \$opt_genomeSize,
    "gff|f=s" => \$gff))

{
    pod2usage( { -message => "Failed to parse command line",
                 -verbose => 1,
                 -exitval => 1 } );
}

# Print Help and exit
if ($help) {
    pod2usage( { -verbose => 1,
                 -exitval => 0,
                 -message => "$header\n" } );
}

if ( ! (defined($gff)) ){
    pod2usage( {
           -message => "$header\nAt least 1 parameter is mandatory:\nInput reference gff file (--gff) \n\n",
           -verbose => 0,
           -exitval => 2 } );
}

#### IN / OUT
my $out = IO::File->new();
if ($opt_output) {

  if (-f $opt_output){
      print "Cannot create a directory with the name $opt_output because a file with this name already exists.\n";exit();
  }
  if (-d $opt_output){
      print "The output directory choosen already exists. Please geve me another Name.\n";exit();
  }

  open($out, '>', $opt_output) or die "Could not open file '$opt_output' $!";
  }
else{
  $out->fdopen( fileno(STDOUT), 'w' );
}

#Manage plot folder output
if($opt_plot){
  if ($opt_output){
    $opt_plot = $opt_output."_distribution_plots";
  }
  else{
    $opt_plot = "distribution_plots";

    if (-f $opt_plot){
      print "Cannot create a directory with the name $opt_plot because a file with this name already exists.\n";exit();
    }
    if (-d $opt_plot){
      print "The default output directory $opt_plot use to save the distribution plots already exists. Please give me another folder name.\n";exit();
    }
  }

  # Check R is available. If not we try to load it through Module software
  if ( system("R --version 1>/dev/null 2>/dev/null") == 0 ) {
    print "R is available. We can continue\n";
  }
  else {
    print "R no available. We cannot perform any plot\n";
    $opt_plot = undef;
  }
}

                #####################
                #     MAIN          #
                #####################



######################
### Parse GFF input #
print "Reading file $gff\n";
my ($hash_omniscient, $hash_mRNAGeneLink) =  slurp_gff3_file_JD({
                                                               input => $gff,
                                                               verbose => 1
                                                               });
print "Parsing Finished\n";
### END Parse GFF input #
#########################

#check number of level1
my $nbLevel1 = 0;
foreach my $tag_l1 (keys %{$hash_omniscient->{'level1'}}){
  $nbLevel1 += keys %{$hash_omniscient->{'level1'}{$tag_l1}};
}

#chech number of level2
my $nbLevel2 = keys %$hash_mRNAGeneLink;

##############
# STATISTICS #
my $stat;
my $distri;
if($opt_genomeSize){
  ($stat, $distri) = gff3_statistics($hash_omniscient, $opt_genomeSize);
}
else{
  ($stat, $distri) = gff3_statistics($hash_omniscient);
}

#print statistics
foreach my $infoList (@$stat){
  foreach my $info (@$infoList){
    print $out "$info";
  }
  print $out "\n";
}

#Check if we have isoforms
if($nbLevel1 != $nbLevel2){

  #print distribution before removing isoforms
  if($opt_plot){
    print_distribution($opt_plot, "with_isoforms", $distri);
  }

  print $out "\nApparently we have isoforms : Number of level1 features: $nbLevel1 / Number of level2 features: $nbLevel2\n";
  print $out "We will proceed to the statistics analysis using only the mRNA with the longest cds\n";

  #create list of level2 where we kept only level2 that have cds and only the longest isoform !
  my $list_id_l2 = get_longest_cds_level2($hash_omniscient);

  # create a new omniscient with only one mRNA isoform per gene
  my $omniscientNew = create_omniscient_from_idlevel2list($hash_omniscient, $hash_mRNAGeneLink, $list_id_l2);

  # print stats
  my $stat;
  my $distri;
  if($opt_genomeSize){
    ($stat, $distri) = gff3_statistics($omniscientNew, $opt_genomeSize);
  }else{
    ($stat, $distri) = gff3_statistics($omniscientNew);
  }

  #print statistics
  foreach my $infoList (@$stat){
    foreach my $info (@$infoList){
      print $out "$info";
    }
    print $out "\n";
  }

  #print distribution after having removed the isoforms
  if($opt_plot){
    print_distribution($opt_plot, "without_isoforms", $distri);
  }
}
else{ #No isoforms
  if($opt_plot){
    print_distribution($opt_plot, "without_isoforms", $distri);
  }
}

# END STATISTICS #
##################
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

sub print_distribution{
  my ($folder, $subfolder, $distri)=@_;

  foreach my $type (keys %{$distri} ) {

    foreach my $level (keys %{$distri->{$type}} ) {
      foreach my $tag ( keys %{$distri->{$type}{$level}} ) {
        if( exists_keys ($distri,($type, $level, $tag, 'whole') ) ){

          if(! -d $folder){
            mkdir $folder;
          }

          if(! -d $folder."/".$subfolder){
            mkdir $folder."/".$subfolder;
          }

          my $outputPDF = $folder."/".$subfolder."/".$type."Class_".$tag.".pdf";

          #CREATE THE R COMMAND
          my $nbValues = @{$distri->{$type}{$level}{$tag}{'whole'}};
          my $R_command = rcc_plot_from_list($distri->{$type}{$level}{$tag}{'whole'}, "", "histogram", "$tag"." size (nt)", "Number of $tag", "Distribution of $tag sizes\nMade with $nbValues $tag"."s", $outputPDF);
          #EXECUTE THE R COMMAND
          execute_R_command($R_command);
        }

        if( exists_keys ($distri,($type, $level, $tag, 'piece') ) ){
        }

      }
    }
  }
}

__END__

=head1 NAME

gff3_statistics.pl -
The script take a gff3 file as input. -
The script give basic statistics of a gff file.
Remark: identical feature from level1 or level2 with identical ID will be merged as well as their subsequent features (Level2 or level3).

=head1 SYNOPSIS

    ./gff3_statistics.pl --gff file.gff  [ -o outfile ]
    ./gff3_statistics.pl --help

=head1 OPTIONS

=over 8

=item B<--gff> or B<-f>

Input GFF3 file that will be read (and sorted)

=item B<--gs> or B<-g>

This option inform about the genome size in oder to compute more statistics. You can give the size in Nucleotide or directly the fasta file.

=item B<-d> or B<-p>

When this option is used, an histogram of distribution of the features will be printed in pdf files. (d means distribution, p means plot).

=item B<--output> or B<-o>

File where will be written the result. If no output file is specified, the output will be written to STDOUT.

=item B<-h> or B<--help>

Display this helpful text.

=back

=cut
