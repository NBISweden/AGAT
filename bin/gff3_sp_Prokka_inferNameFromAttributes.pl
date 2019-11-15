#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Clone 'clone';
use Getopt::Long;
use Pod::Usage;
use List::MoreUtils qw(uniq);
use Bio::Tools::GFF;
use Agat::Omniscient;

my $header = get_agat_header();
my $gff = undef;
my $help= 0;
my $force=undef;
my $outfile=undef;

if ( !GetOptions(
    "help|h" => \$help,
    "gff|f=s" => \$gff,
    "force" => \$force,
    "output|outfile|out|o=s" => \$outfile))

{
    pod2usage( { -message => 'Failed to parse command line',
                 -verbose => 1,
                 -exitval => 1 } );
}

# Print Help and exit
if ($help) {
    pod2usage( { -verbose => 2,
                 -exitval => 2,
                 -message => "$header\n" } );
}

if ( ! (defined($gff)) ){
    pod2usage( {
           -message => "$header\nAt least 1 parameter is mandatory:\nInput reference gff file (--gff) \n\n",
           -verbose => 0,
           -exitval => 2 } );
}

my $gffout;
if ($outfile) {
  $outfile=~ s/.gff//g;
  open(my $fh, '>', $outfile.".gff") or die "Could not open file '$outfile' $!";
  $gffout= Bio::Tools::GFF->new(-fh => $fh, -gff_version => 3 );
}
else{
  $gffout = Bio::Tools::GFF->new(-fh => \*STDOUT, -gff_version => 3);
}



                #####################
                #     MAIN          #
                #####################


######################
### Parse GFF input #
my ($hash_omniscient, $hash_mRNAGeneLink) = slurp_gff3_file_JD({ input => $gff
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


        #foreach my $val(@inferenceAtt){
        #  print "ok".$val."\n";
        #}
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
      #else{
      #  print "We skip: ".$feature->gff_string."\n";
      #}
    }
  }
}
print "We added $nbNameAdded Name attributes\n";

print_omniscient($hash_omniscient, $gffout); #print gene modified
__END__

=head1 NAME

gff3_addProkkaNameFromInferenceAttribute.pl -
The script take a gff3 file as input. -
The script give basic statistics of a gff file.
Remark: identical feature from level1 or level2 with identical ID will be merged as well as their subsequent features (Level2 or level3).

=head1 SYNOPSIS

    ./gff3_addProkkaNameFromInferenceAttribute.pl -gff file.gff  [ -o outfile ]
    ./gff3_addProkkaNameFromInferenceAttribute.pl --help

=head1 OPTIONS

=over 8

=item B<--gff> or B<-f>

Input GFF3 file that will be read (and sorted)

=item B<--force>

If Name attribute already exists, they will be replaced if a new one is found

=item B<-o> , B<--output> , B<--out> or B<--outfile>

Output GFF file.  If no output file is specified, the output will be
written to STDOUT.

=item B<-h> or B<--help>

Display this helpful text.

=back

=cut
