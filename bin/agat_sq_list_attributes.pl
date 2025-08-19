#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Clone 'clone';
use Getopt::Long;
use Pod::Usage;
use IO::File;
use List::MoreUtils qw(uniq);
use AGAT::AGAT;

my $header = get_agat_header();
my $config;
my $start_run = time();
my %handlers;
my $gff = undef;
my $opt_help= 0;
my $primaryTag=undef;
my $outfile=undef;

if ( !GetOptions(
    'c|config=s'             => \$config,
    "h|help"                 => \$opt_help,
    "gff|f=s"                => \$gff,
    "p|t|l=s"                => \$primaryTag,
    "output|outfile|out|o=s" => \$outfile))

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

if ( ! (defined($gff)) ){
    pod2usage( {
           -message => "$header\nAt least 1 parameter is mandatory:\nInput reference gff file (--gff) \n\n",
           -verbose => 0,
           -exitval => 2 } );
}

# --- Manage config ---
$config = get_agat_config({config_file_in => $config});

my $log;
if ($config->{log}) {
  my ($file) = $0 =~ /([^\/]+)$/;
  my $log_name = $file . ".agat.log";
  open($log, '>', $log_name) or die "Can not open $log_name for printing: $!";
  dual_print($log, $header, 0);
}

# Manage $primaryTag
my @ptagList;
if(! $primaryTag or $primaryTag eq "all"){
  dual_print($log, "We will work on attributes from all features\n");
  push(@ptagList, "all");
}
else{
   @ptagList= split(/,/, $primaryTag);
   foreach my $tag (@ptagList){
      dual_print($log, "We will work on attributes from $tag feature.\n");
   }
}

# Manage input fasta file
my $format = $config->{force_gff_input_version};
if(! $format ){ $format = select_gff_format($gff); }
my $ref_in = AGAT::BioperlGFF->new(-file => $gff, -gff_version => $format);


                #####################
                #     MAIN          #
                #####################

my %all_attributes;
my %attributes_per_level;
######################
### Parse GFF input #

#time to calcul progression
my $startP=time;
my $nbLine=`wc -l < $gff`;
$nbLine =~ s/ //g;
chomp $nbLine;
dual_print($log, "$nbLine line to process...\n");

my $geneName=undef;
my $line_cpt=0;
while (my $feature = $ref_in->next_feature() ) {
  $line_cpt++;


    manage_attributes($feature, \@ptagList, \%all_attributes, \%attributes_per_level);

  #Display progression
  if ((30 - (time - $startP)) < 0) {
    my $done = ($line_cpt*100)/$nbLine;
    $done = sprintf ('%.0f', $done);
        dual_print($log, "\rProgression : $done % processed.\n");
    $startP= time;
  }
}

#print "We added $nbNameAdded Name attributes\n";

my $out = prepare_fileout($outfile);

# Print information by feature
my $nbFeat = scalar keys %attributes_per_level;
print $out "\nWe met ".$nbFeat." different feature types.";
foreach my $feature_type ( sort keys %attributes_per_level){
  my $nbAtt = scalar keys %{$attributes_per_level{$feature_type}};
  print $out "\nHere the list of all the attributes tags met for the feature type <".$feature_type."> (".$nbAtt." attributes):\n";
  foreach my $attribute ( sort keys %{$attributes_per_level{$feature_type}} ){
    print $out $attribute."\n";
  }
}

# Print Global information
my $nbAtt = scalar keys %all_attributes;
print $out "\nHere the list of all the attributes tags met (".$nbAtt." attributes):\n";
foreach my $attribute ( sort keys %all_attributes){
  print $out $attribute."\n";
}


##Last round
my $end_run = time();
my $run_time = $end_run - $start_run;
dual_print($log, "\nJob done in $run_time seconds\n");

close $log if $log;

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

sub  manage_attributes{
  my  ($feature, $ptagList, $all_attributes, $attributes_per_level)=@_;

  my $primary_tag=$feature->primary_tag;

  # check primary tag (feature type) to handle
  foreach my $ptag (@$ptagList){

    if($ptag eq "all"){
      tag_from_list($feature,$all_attributes, $attributes_per_level);
    }
    elsif(lc($ptag) eq lc($primary_tag) ){
      tag_from_list($feature,$all_attributes, $attributes_per_level);
    }
  }
}

sub tag_from_list{
  my  ($feature, $all_attributes, $attributes_per_level)=@_;

  foreach my $tag ($feature->get_all_tags) {
      # create handler if needed (on the fly)
      if(! exists_keys( $all_attributes,($tag) ) ) {
        $all_attributes{$tag}++;
      }
      if(! exists_keys ( $attributes_per_level,($feature->primary_tag,$tag) ) ) {
        $attributes_per_level{$feature->primary_tag}{$tag}++;
      }

  }
}


__END__


=head1 NAME

agat_sq_list_attributes.pl

=head1 DESCRIPTION

The script take a gff3 file as input. -
The script give information about attribute tags used within you file.

=head1 SYNOPSIS

    agat_sq_list_attributes.pl -gff file.gff -p level2,cds,exon [ -o outfile ]
    agat_sq_list_attributes.pl --help

=head1 OPTIONS

=over 8

=item B<--gff> or B<-f>

Input GTF/GFF file.

=item B<-p>,  B<-t> or  B<-l>

primary tag option, case insensitive, list. Allow to specied the feature types that will be handled.
You can specified a specific feature by given its primary tag name (column 3) as: cds, Gene, MrNa
You can specify directly all the feature of a particular level:
      level2=mRNA,ncRNA,tRNA,etc
      level3=CDS,exon,UTR,etc
By default all feature are taking in account. fill the option by the value "all" will have the same behaviour.

=item B<-o> , B<--output> , B<--out> or B<--outfile>

Output GFF file.  If no output file is specified, the output will be
written to STDOUT.

=item B<-c> or B<--config>

String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any, 
otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

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
