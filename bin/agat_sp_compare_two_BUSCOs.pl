#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Clone 'clone';
use Pod::Usage;
use Getopt::Long;
use IO::File ;
use List::Util 'first';
use AGAT::AGAT;

my $header = get_agat_header();
my $config;
my $threads;
my $start_run = time();
my $folderIn1=undef;
my $folderIn2=undef;
my $outfolder=undef;
my $verbose=undef;
my $opt_help = 0;


Getopt::Long::Configure ('bundling');
if ( !GetOptions ('f1=s'        => \$folderIn1,
                  "f2=s"        => \$folderIn2,
                  'o|output=s'  => \$outfolder,
                  'v|verbose=i' => \$verbose,
                  'c|config=s'  => \$config,
                  'h|help!'     => \$opt_help ) )
{
    pod2usage( { -message => 'Failed to parse command line',
                 -verbose => 1,
                 -exitval => 1 } );
}

if ($opt_help) {
    pod2usage( { -verbose => 99,
                 -exitval => 0,
                 -message => "$header\n" } );
}

if ( !defined($folderIn1) or  !defined($folderIn2) ){
   pod2usage( {  -message => "$header\nAt least 2 parameters are mandatory: --f1 and --f2",
                 -verbose => 0,
                 -exitval => 2 } );
}

# --- Manage config ---
initialize_agat({ config_file_in => $config, input => $folderIn1 });
$CONFIG->{threads} = $threads if defined($threads);

# Manage input folder1
my $fh1;
$folderIn1 = remove_slash_path_folder($folderIn1);
opendir(DIR, "$folderIn1")  or die "Unable to read Directory : $!";
my @files_table1 = grep(/^full_table/,readdir(DIR));
if (! @files_table1){print "full_table[_abinitio].tsv file missing in $folderIn1\n"; exit;}
my $path1=$folderIn1."/".$files_table1[0];
open($fh1, '<', $path1) or die "Could not open file '$path1' $!";

#Manage input folder2
my $fh2;
$folderIn2 = remove_slash_path_folder($folderIn2);
opendir(DIR, "$folderIn2") or die "Unable to read Directory $folderIn2 : $!";;
my @files_table2 = grep(/^full_table/,readdir(DIR));
if (! @files_table2){print "full_table[_abinitio].tsv file missing in $folderIn2\n"; exit;}
my $path2=$folderIn2."/".$files_table2[0];
open($fh2, '<', $path2) or die "Could not open file '$path2' $!";


#Manage output folder
if ($outfolder) {
  $outfolder = remove_slash_path_folder($outfolder);
  if(! -d $outfolder ){
    mkdir $outfolder;
  }
  else{
    print "$outfolder output folder already exists !\n"; exit;
  }
}

# Manage Output gff files
my $gffout_complete_path;
my $gffout_fragmented_path;
my $gffout_duplicated_path;
if ($outfolder) {
  $gffout_complete_path = $outfolder."/"."f1_complete.gff";
	$gffout_fragmented_path = $outfolder."/"."f1_fragmented.gff";
	$gffout_duplicated_path = $outfolder."/"."f1_duplicated.gff";
}
my $gffout_complete = prepare_gffout( $gffout_complete_path);
my $gffout_fragmented = prepare_gffout( $gffout_fragmented_path);
my $gffout_duplicated = prepare_gffout( $gffout_duplicated_path);

my %gff_out;
$gff_out{'complete'}=$gffout_complete;
$gff_out{'fragmented'}=$gffout_fragmented;
$gff_out{'duplicated'}=$gffout_duplicated;

#############################################################
#                         MAIN
#############################################################

#Read busco1 file
my %busco1;
while( my $line = <$fh1>)  {

  if( $line =~ m/^\w+\s{1}Complete/){
    my @list = split(/\s/,$line);
    $busco1{'complete'}{$list[0]}=$line;
  }
  if( $line =~ m/^\w+\s{1}Missing/){
    my @list = split(/\s/,$line);
    $busco1{'missing'}{$list[0]}=$line;
  }
  if( $line =~ m/^\w+\s{1}Fragmented/){
    my @list = split(/\s/,$line);
    $busco1{'fragmented'}{$list[0]}=$line;
  }
  if( $line =~ m/^\w+\s{1}Duplicated/){
    my @list = split(/\s/,$line);
    $busco1{'duplicated'}{$list[0]}=$line;
  }
}

#Read busco2 file
my %busco2;
while( my $line = <$fh2>)  {

  if( $line =~ m/^\w+\s{1}Complete/){
    my @list = split(/\s/,$line);
    $busco2{'complete'}{$list[0]}=$line;
  }
  if( $line =~ m/^\w+\s{1}Missing/){
    my @list = split(/\s/,$line);
    $busco2{'missing'}{$list[0]}=$line;
  }
  if( $line =~ m/^\w+\s{1}Fragmented/){
    my @list = split(/\s/,$line);
    $busco2{'fragmented'}{$list[0]}=$line;
  }
  if( $line =~ m/^\w+\s{1}Duplicated/){
    my @list = split(/\s/,$line);
    $busco2{'duplicated'}{$list[0]}=$line;
  }
}

my %hashCases;
my %streamOutputs;
#compare busco1 and busco2
foreach my $type1 (keys %busco1){
  foreach my $id1 (keys %{$busco1{$type1}} ){

    foreach my $type2 (keys %busco2){
      if($type1 ne $type2){
        if(exists_keys (\%busco2,($type2,$id1)  ) ){

          my $name=$type1."2".$type2;
          $hashCases{$id1}=$name;
          # create streamOutput
          if($outfolder){
            if (! exists_keys (\%streamOutputs,($name)) ){
              my $ostream = IO::File->new();
              $ostream->open( $outfolder."/$name.txt", 'w' ) or croak( sprintf( "Can not open '%s' for writing %s", $outfolder."/$name.txt", $! ) );
              $streamOutputs{$name}=$ostream;
            }
            my $streamOut=$streamOutputs{$name};
            print $streamOut  $busco1{$type1}{$id1};
          }
          else{
            print "$id1 was $type1 and it is now $type2\n";
          }
        }
      }
    }
    if(! exists_keys(\%hashCases,($id1) ) ){
      $hashCases{$id1}=$type1."2".$type1;
    }
  }
}

#extract gff from folder1
my $full_omniscient={};
my $loop = 0;
my $list_uID_new_omniscient=undef;
my $augustus_gff_folder=$folderIn1."/augustus_output/predicted_genes";

if (-d $augustus_gff_folder){
  opendir(DH, $augustus_gff_folder);
  my @files = readdir(DH);

  my %track_found;
  my @list_cases=("complete","fragmented","duplicated");
  foreach my $type (@list_cases){
    print "extract gff for $type cases\n" if $verbose;
    foreach my $id (sort keys %{$busco1{$type}}){
      my @list = split(/\s/,$busco1{$type}{$id});
      my $seqId = $list[2];
      my $start = $list[3];
      my $end = $list[4];

      my @matches = grep { /\Q$id/ } @files;
      if( @matches){
        foreach my $match (sort @matches){
          my $path = $augustus_gff_folder."/".$match;
          if (-f $path ){
            my  $found=undef;
            print $path."\n" if $verbose;

            my ($hash_omniscient, $hash_mRNAGeneLink) = slurp_gff3_file_JD({ input => $path,
                                                                             config => $config
                                                                        });
            if (!keys %{$hash_omniscient}){
              print "No gene found for $path\n";exit;
            }

            my @listIDl1ToRemove;
            if( exists_keys ($hash_omniscient,('level1','gene'))){
              foreach my $id_l1 (keys %{$hash_omniscient->{'level1'}{'gene'}}){
                my $feature = $hash_omniscient->{'level1'}{'gene'}{$id_l1};
                if ($feature->seq_id() eq $seqId and  $feature->start == $start and $feature->end == $end){
                  $found=1;
                  $track_found{$type}{$id}++;

                  #Add the OG name to the feature, to be displayed in WA
                  foreach my $tag_l2 (keys %{$hash_omniscient->{'level2'}}){
                    if( exists_keys($hash_omniscient,('level2', $tag_l2, $id_l1))){
                      foreach my $feature_l2 ( @{$hash_omniscient->{'level2'}{$tag_l2}{$id_l1}} ){
                        my $value=$id."-".$hashCases{$id};
                        $feature_l2->add_tag_value('description', $value);
                      }
                    }
                  }
                }
                else{push(@listIDl1ToRemove,$id_l1);}
              }

              if ($found){
                if(@listIDl1ToRemove){
                  print "lets remove those supernumary annotation: @listIDl1ToRemove \n" if $verbose;
                  remove_omniscient_elements_from_level1_id_list($hash_omniscient, \@listIDl1ToRemove);
                }

                if($loop == 0){
                  $full_omniscient = clone($hash_omniscient);
                  $loop++;
                }
                elsif($loop == 1){
                  ( $full_omniscient, $list_uID_new_omniscient) = merge_omniscients($full_omniscient, $hash_omniscient);
                  $loop++;
                }
                else{
                  ( $full_omniscient, $list_uID_new_omniscient) = merge_omniscients($full_omniscient, $hash_omniscient, $list_uID_new_omniscient);
                }
              }
              else{
                print "No annotation as described in the tsv file found in the gff file $path\n" if $verbose;
              }
            }
            else{
              print "No annotation in the file $path, lets look the next one.\n" if $verbose;
            }
          }
          else{
            print "A) file $id not found among augustus gff output\n" if $verbose;
          }
        }
      }
      else{
        print "file $id not found among augustus gff output\n" if $verbose;
      }
      if(! exists_keys(\%track_found,($type,$id))){
        print "WARNING After reading all the files related to id $id we didn't found any annotation matching its described in the tsv file.\n";
      }
    }
    my $out = $gff_out{$type};
    print_omniscient( {omniscient => $full_omniscient, output => $out} );
    %$full_omniscient = (); # empty hash
    $list_uID_new_omniscient=undef; #Empty Id used;
    my $nb = keys %{$track_found{$type}};
    $loop = 0;
    print "We found $nb annotations from $type busco\n";
  }

}
else{ print "$augustus_gff_folder folder doesn't exits\n"; exit;}




##Last round
my $end_run = time();
my $run_time = $end_run - $start_run;
print "Job done in $run_time seconds\n";
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

sub remove_slash_path_folder{
  my ($folder_path)=@_;
  if ( $folder_path =~ /\/$/){
    return  $folder_path = substr $folder_path, 0, -1;
  }
  else{
    return $folder_path;
  }
}

__END__


=head1 NAME

agat_sp_compare_two_BUSCOs.pl

=head1 DESCRIPTION

The tool compares the results from two BUSCO runs (genome and proteome mode) in order to pinpoint the differences.
It compares the BUSCOs classification (complete,fragmented, duplicated) of the 1st run (genome mode)
against the classification found in the second run. It will report the results in txt files, and
extracts the complete,fragmented and duplicated annotated BUSCOs from the 1st run in gff files.
We add in the gff an attribute specifying the cases e.g. description=EOG090W00UK-complete2duplicated.
Where EOG090W00UK is the BUSCO name/label/group investigated, and complete2duplicated the case we found
(was complete in run1 and duplicated in run2).
By loading these gff tracks in a web browser and helped by other tracks (e.g the genome annotation/prediction)
can help to understand why the BUSCO have been classified differently from run1 to run2.
In other term it allows to catch potential problems in an annotation.
agat_sp_compare_two_BUSCOs.pl has been tested with results from BUSCO version 3 and 4.
/!\ The tool expects a BUSCO run in genome mode as input folder 1 and a BUSCO run in proteins mode
as input folder 2. You can also decide to provide twice (--f1 --f2) the same BUSCO run in genome mode,
the tool will only extract the annotation of the complete,fragmented and duplicated annotated BUSCOs from the 1st run in gff.

=head1 SYNOPSIS

    agat_sp_compare_two_BUSCOs.pl --f1 <input busco folder1> --f2 <input busco folder2> [-o <output folder>]
    agat_sp_compare_two_BUSCOs.pl --help

=head1 OPTIONS

=over 8

=item B<--f1>

STRING: Input busco folder1

=item B<--f2>

STRING: Input busco folder2

=item B<-v> or B<--verbose>

Integer: For displaying extra information use -v 1.

=item B<-o> or B<--output>

STRING: Output folder.

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
