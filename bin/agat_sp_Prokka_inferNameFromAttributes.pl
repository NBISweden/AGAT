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
my $force=undef;
my $outfile=undef;

if ( !GetOptions(
    'c|config=s'               => \$config,
    "h|help" => \$opt_help,
    "gff|f=s" => \$gff,
    "force" => \$force,
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
           -exitval => 1 } );
}

my $gffout = prepare_gffout($config, $outfile);

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>     MAIN     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

######################
### Parse GFF input #
my ($hash_omniscient, $hash_mRNAGeneLink) = slurp_gff3_file_JD({ input => $gff,
                                                                 config => $config
                                                              });
print ("GFF3 file parsed\n");

my $nbNameAdded=0;

foreach my $tag (keys %{$hash_omniscient->{'level1'}}){
  foreach my $id (keys %{$hash_omniscient->{'level1'}{$tag}}){

    my $feature=$hash_omniscient->{'level1'}{$tag}{$id};

    if($feature->has_tag('name')){
      my $name=$feature->_tag_value('name');
      create_or_replace_tag($feature,'Name', $name);
      $feature->remove_tag('name');
    }

    #Name already contained in the gene attribute.
    if($feature->has_tag('gene')){
      # we get the Name
      my $name=$feature->_tag_value('gene');

      # If no attribute Name or if we have to replace it
      if(! $feature->has_tag('Name') or ($force)){
        create_or_replace_tag($feature,'Name', $name);
        $nbNameAdded++;
      }
      elsif($feature->has_tag('Name') and ( ! $force)){
        print "Feature contains already an attribute Name. You can force it replacement by using the option --force\n";
      }
      print "Name found in gene attribute = $name\n";
    }# Name not found in gene attribute. So we try to get the name included in the inference attribute.
    elsif($feature->has_tag('inference')){
      my @inferenceAtt=$feature->get_tag_values('inference');
      if ($#inferenceAtt > 0){

        my @tab = split /\|/,$inferenceAtt[$#inferenceAtt]; # split the last value by the character |
        my $name = $tab[$#tab];

        # SKIP case
        if($name =~ /protein motif:Pfam:/i){
          next;
        }
        if($name =~ /protein motif:CLUSTERS:/i){
          next;
        }
        if($name =~ /similar to AA sequence:UniProtKB:/i){
          next;
        }
        # ELSE name contains the Uniprot header of the protein coming from "--proteins" option ( Fasta file of trusted proteins to first annotate from ).


        if(! $feature->has_tag('Name') or ($force)){
          create_or_replace_tag($feature,'Name', $name);
          $nbNameAdded++;
        }
        elsif($feature->has_tag('Name') and ( ! $force)){
          print "Feature contains already an attribute Name. You can force it replacement by using the option --force\n";
        }
        print "My Name get in inference attribute = $name\n";
      }
    }
  }
}
print "We added $nbNameAdded Name attributes\n";

print_omniscient( {omniscient => $hash_omniscient, output => $gffout} );

__END__

=head1 NAME

agat_sp_Prokka_inferNameFromAttributes.pl

=head1 DESCRIPTION

The script aims to fill a Name attribute based on <gene> attribute in a prokka gff
annotation file. If no gene attribute is present it take if from the <inference>
attribute.

=head1 SYNOPSIS

    agat_sp_Prokka_inferNameFromAttributes.pl -gff file.gff  [ -o outfile ]
    agat_sp_Prokka_inferNameFromAttributes.pl --help

=head1 OPTIONS

=over 8

=item B<--gff> or B<-f>

Input GTF/GFF file.

=item B<--force>

If Name attribute already exists, they will be replaced if a new one is found

=item B<-o> , B<--output> , B<--out> or B<--outfile>

Output GFF file.  If no output file is specified, the output will be
written to STDOUT.

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
