#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Clone 'clone';
use Getopt::Long;
use Pod::Usage;
use List::MoreUtils qw(uniq);
use AGAT::AGAT;

my $header = get_agat_header();
my $config;
my $gff = undef;
my $opt_help= 0;
my $attribute='transcript_id';
my $start_run = time();
my $outfile=undef;
my $cpt_case=0;

if ( !GetOptions(
    'c|config=s'               => \$config,
    "h|help"      => \$opt_help,
    "gff|f=s"     => \$gff,
    "tag|att=s"   => \$attribute,
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

if ( ! $gff ){
    pod2usage( {
           -message => "$header\nAt least 1 parameter is mandatory:\nInput reference gff file (--gff)\n\n",
           -verbose => 0,
           -exitval => 2 } );
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>     MAIN     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

print "Looking to $attribute attribute.\n";

# Manage input fasta file
my $format = $config->{gff_output_version};
if(! $format ){ $format = select_gff_format($gff); }
my $ref_in = Bio::Tools::GFF->new(-file => $gff, -gff_version => $format);

#time to calcul progression
my $startP=time;
my $nbLine=`wc -l < $gff`;
$nbLine =~ s/ //g;
chomp $nbLine;
print "$nbLine line to process...\n";

my $line_cpt=0;
my %hash_values;
my $nb_attributes=0;
while (my $feature = $ref_in->next_feature() ) {
  $line_cpt++;

  if ($feature->has_tag($attribute) ){
    $hash_values{$feature->_tag_value($attribute)}++;
    $nb_attributes++;
  }


  #####################
  #Display progression
  if ((30 - (time - $startP)) < 0) {
    my $done = ($line_cpt*100)/$nbLine;
    $done = sprintf ('%.0f', $done);
        print "\rProgression : $done % processed.\n";
    $startP= time;
  }
}

##Last round
my $end_run = time();
my $run_time = $end_run - $start_run;


my $result = scalar keys %hash_values;

print "$line_cpt features read. Among them, $nb_attributes has the $attribute attribute.\n".
"There is $result unique value within $attribute attribute\n";
print "Job done in $run_time seconds\n";

__END__

=head1 NAME

agat_sq_count_attributes.pl

=head1 DESCRIPTION

The script count the number of a choosen attribute and also count the number of
unique value of this attribute.

=head1 SYNOPSIS

    agat_sq_count_attributes.pl --gff file.gff  --att gene_id [ -o outfile ]
    agat_sq_count_attributes.pl --help

=head1 OPTIONS

=over 8

=item B<--gff> or B<-f>

Input GTF/GFF file.

=item B<--tag>, B<--att>

The name of the attribute that will be investigated.

=item B<-o> , B<--output> , B<--out> or B<--outfile>

Output GFF file.  If no output file is specified, the output will be
written to STDOUT.

=item B<-c> or B<--config>

String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any, 
otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
The --config option gives yo the possibility to use your own AGAT config file (located elsewhere or named differently).

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
