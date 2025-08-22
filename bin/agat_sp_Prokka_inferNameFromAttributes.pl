#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Clone 'clone';
use Getopt::Long::Descriptive;
use Pod::Usage;
use List::MoreUtils qw(uniq);
use AGAT::AGAT;

my $header = get_agat_header();
my ( $opt, $usage, $config ) = AGAT::AGAT::describe_script_options(
    $header,
    [ 'gff|f=s', 'Input GFF file', { required => 1 } ],
    [ 'force',   'Replace existing Name attributes' ],
);

my $gff     = $opt->gff;
my $force   = $opt->force;
my $outfile = $opt->out;

my $log;
if ( my $log_name = $config->{log_path} ) {
    open( $log, '>', $log_name ) or die "Can not open $log_name for printing: $!";
    dual_print( $log, $header, 0 );
}

# Prepare output
my $gffout = prepare_gffout($config, $outfile);

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>     MAIN     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

######################
### Parse GFF input #
my ($hash_omniscient, $hash_mRNAGeneLink) = slurp_gff3_file_JD({ input => $gff,
                                                                 config => $config
                                                              });
dual_print($log, "GFF3 file parsed\n");

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
        dual_print($log, "Feature contains already an attribute Name. You can force it replacement by using the option --force\n");
      }
      dual_print($log, "Name found in gene attribute = $name\n");
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
          dual_print($log, "Feature contains already an attribute Name. You can force it replacement by using the option --force\n");
        }
        dual_print($log, "My Name get in inference attribute = $name\n");
      }
    }
  }
}
dual_print($log, "We added $nbNameAdded Name attributes\n");

close $log if $log;

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
