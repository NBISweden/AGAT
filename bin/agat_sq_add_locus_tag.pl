#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Pod::Usage;
use Getopt::Long::Descriptive;
use IO::File;
use AGAT::AGAT;

my $start_run = time();
my $header = get_agat_header();
my @copyARGV = @ARGV;
my ( $opt, $usage, $config ) = AGAT::AGAT::describe_script_options( $header,
    [ 'file|input|gff=s', 'Input GFF file', { required => 1 } ],
    [ 'to|lo=s',          'Name of locus tag attribute', { default => 'locus' } ],
    [ 'ti|li=s',          'Input attribute to derive locus tag' ],
    [ 'p|type|l=s',       'Comma-separated list of feature types' ],
);

my $inputFile  = $opt->file;
my $locus_tag  = $opt->to;
my $tag_in     = $opt->ti;
my $primaryTag = $opt->p;
my $locus_cpt  = 1;

my $log;
if ( my $log_name = $config->{log_path} ) {
    open( $log, '>', $log_name ) or die "Can not open $log_name for printing: $!";
    dual_print( $log, $header, 0 );
}

# Manage input gff file
my $format = $config->{force_gff_input_version};
if(! $format ){ $format = select_gff_format($inputFile); }
my $ref_in = AGAT::BioperlGFF->new(-file => $inputFile, -gff_version => $format);

# Manage Output
my $outfile = $config->{output};
my $gffout = prepare_gffout($config, $outfile);

#define the locus tag
if(! $locus_tag){
  $locus_tag="locus_tag";
}

# Manage $primaryTag
my @ptagList;
my $hash_levels= get_levels_info();
my $hash_level1 = $hash_levels->{'other'}{'level'}{'level1'};

if(! $primaryTag){
  dual_print($log, "We will work on attributes from all Level1 features.\n", $config->{verbose});
  push(@ptagList, "all");
}
else{
   @ptagList= split(/,/, $primaryTag);
   foreach my $tag (@ptagList){
      if ( exists_keys ( $hash_level1, ( lc($tag) ) ) ){
        dual_print($log, "We will work on attributes from <$tag> feature.\n", $config->{verbose});
      }
      else{
        dual_print($log, "<$tag> feature is not a level1 feature. Current accepted value are:\n", $config->{verbose});
        foreach my $key ( keys %{$hash_level1}){
          dual_print($log, $key." ", $config->{verbose});
        }
        dual_print($log, "\n", $config->{verbose});
        exit;
      }
   }
}

#time to calcul progression
my $startP=time;
my $nbLine=`wc -l < $inputFile`;
$nbLine =~ s/ //g;
chomp $nbLine;
dual_print($log, "$nbLine line to process...\n", $config->{verbose});

my $line_cpt=0;
my $locus=undef;
while (my $feature = $ref_in->next_feature() ) {
  $line_cpt++;

  my $ptag = lc($feature->primary_tag());

  if ( exists_keys( $hash_level1, ( $ptag ) ) and  ( $hash_level1->{$ptag} ne "topfeature" ) ){

    # initialize locus_tag
    if ( grep( /^$ptag/, @ptagList ) or   grep( /^all/, @ptagList ) ) {

      # if locus_tag has to be the value of an existing attribute.
      if( $tag_in){
        if( $feature->has_tag($tag_in)){
          $locus = $feature->_tag_value($tag_in);
        }
        else{
          dual_print($log, "No attribute $tag_in for the following feature:\n".$feature->gff_string()."\n", $config->{verbose});
          $locus = $locus_tag.$locus_cpt;$locus_cpt++;
          dual_print($log, "We will use the created locus_tag value: $locus instead to name the locus!\n", $config->{verbose});
        }
      }
      else{
        $locus = $locus_tag.$locus_cpt;$locus_cpt++;
      }
    }
    else{
      $locus=undef;
    }
    # if level1 and part of those to provide locus_tag
    if($locus){
      create_or_replace_tag($feature,$locus_tag, $locus);
    }
  }
  else{
  # if not level 1 and we have to spread locus_tag to sub feature.
    if($locus){
      create_or_replace_tag($feature,$locus_tag, $locus);
    }
  }

  $gffout->write_feature($feature);

  #Display progression
  if ((30 - (time - $startP)) < 0) {
    my $done = ($line_cpt*100)/$nbLine;
    $done = sprintf ('%.0f', $done);
        dual_print($log, "\rProgression : $done % processed.\n", $config->{verbose});
    $startP= time;
  }
}


##Last round
my $end_run = time();
my $run_time = $end_run - $start_run;
dual_print($log, "Job done in $run_time seconds\n", $config->{verbose});


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

agat_sq_add_locus_tag.pl

=head1 DESCRIPTION

Add a shared locus tag per record. A record is all features linked by each other
by parent/children relationship (e.g Gene,mRNA,exon, CDS).

=head1 SYNOPSIS

    agat_sq_add_locus_tag.pl --gff <input file> [-o <output file>]
    agat_sq_add_locus_tag.pl --help

=head1 OPTIONS

=over 8

=item B<--gff>, B<--file> or B<--input>

STRING: Input GTF/GFF file.

=item B<-p>,  B<--type> or  B<-l>

Primary tag option, case insensitive, list. Allow to specied the Level1 feature types that will be handled.
By default all feature Level1 are taken into account.

=item B<--lo> or B<--to>

Locus tag output, by defaut it will be called locus_tag, but using this option you can specied the name of this attribute.

=item B<--li> or B<--ti>

Tag input, by default the value of the locus tag attribute will be locusX where X is an incremented number.
You can use the values of an existing attribute instead e.g the ID value: --li ID.

=item B<-o> or B<--output>

STRING: Output file.  If no output file is specified, the output will be written to STDOUT. The result is in tabulate format.

=item B<-q> or B<--quiet>

To remove verbosity.

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
