#!/usr/bin/env perl

use Carp;
use strict;
use Data::Dumper;
use Getopt::Long;
use File::Basename;
use Statistics::R;
use Pod::Usage;
use Bio::Tools::GFF;
use List::MoreUtils qw(uniq);

my $header = qq{
########################################################
# BILS 2023 - Sweden                                   #
# lucile.soler\@nbis.se                               #
# Please cite NBIS (www.nbis.se) when using this tool. #
########################################################
};

my $output = undef;
my $gff = undef;
my $help= 0;
my $output_gff;
my $ID;
my $description;
my $product;


if ( !GetOptions(
    "help|h" => \$help,
    "gff|g=s" => \$gff,
    "output|out|o=s" => \$output_gff)
)

{
    pod2usage( { -message => "$header"."Failed to parse command line.",
                 -verbose => 1,
                 -exitval => 1} );
}

# Print Help and exit
if ($help) {
    pod2usage( { -message => "$header",
                 -verbose => 2,
                 -exitval => 2 } );
}

if ( ! (defined($gff))  ){
    pod2usage( {
           -message => "$header\nAt least 1 parameters are mandatory.\n",
           -verbose => 0,
           -exitval => 1 } );
}

my $ostream = IO::File->new();
if(defined($output_gff)){
  $ostream->open($output_gff, 'w' ) or
  croak(
    sprintf( "Can not open '%s' for reading: %s", $output_gff, $! ) );
  }
else{
  $ostream->fdopen( fileno(STDOUT), 'w' ) or
    croak( sprintf( "Can not open STDOUT for writing: %s", $! ) );
}

    my $input_fh1 = IO::File->new( $gff, '<' ) or croak "Can't read file '$gff' : $!";

    while (my $line1 = $input_fh1->getline()){ #print {$STDERR} "$.\t$line";
        chomp $line1;

        my @col_gffs1		= split(/\t/, $line1);
    my $chr2			= $col_gffs1[0];
    my $origin = $col_gffs1[1];
    my $feature2		= $col_gffs1[2];
    my $pos_start2	= $col_gffs1[3];
    my $pos_end2		= $col_gffs1[4];
    my $col5 = $col_gffs1[5];
    my $strand2	= $col_gffs1[6];
    my $col7 = $col_gffs1[7];
    my $desc2		= $col_gffs1[8];
    my $id2;

#print $feature2."\n";

    if ($feature2 eq "mRNA"){
       $desc2 =~/ID=(\S+);Parent=\S+;(Dbxref=\S+;)_AED=\S+;makerName=\S+;(.*)/; #;Parent=*;(*)maker_name=*;(*)/;
       $ID = $1;
       $description = $2;
       $product = $3;

       print $ostream $line1."\n";

} elsif ($feature2 eq "CDS"){
      $desc2 =~/ID=\S+;Parent=(\S+);makerName=\S+/; #;Parent=*;(*)maker_name=*;(*)/;

      if ($ID eq $1) {

        print $ostream $line1.";".$description.$product."\n";
      }else {
    #    print "ID".$ID."\n";
        print $ostream $line1.";product=hypothetical protein"."\n";
      }
    }
    else {
            print $ostream $line1."\n";
    }



  }



=head1 NAME

agat_create_ENA_compliant_gff_vmaker.pl

This script will create an ENA compliant GFF coming from a maker functional annotation gff (annotated by the pipeline)
It copies the functional information from the mRNA into the CDS and keeping the same ID between the two gffs

=head1 SYNOPSIS

    ./agat_create_ENA_compliant_gff_vmaker.pl -gff annotated.gff -o outputname [Options]
    ./agat_create_ENA_compliant_gff_vmaker.pl --help

=head1 OPTIONS

=over 8

=item B<--gff>, B<--gfffile> or B<-f>

Input GFF3 file corresponding to the file annotated by the functional annotation pipeline and coming from maker.


=item B<--output> or B<-o>

Output name that will be used to create the output files.

=item B<--help> or B<-h>

Getting help.
Display the full information.

=back

=cut
