#!/usr/bin/perl -w

package AGAT::Utilities;

use strict;
use warnings;
use Bio::Tools::GFF;
use Bio::Seq;
use Clone 'clone';
use Sort::Naturally;
use Exporter;

our @ISA = qw(Exporter);
our @EXPORT = qw(exists_keys exists_undef_value get_proper_codon_table);
sub import {
  AGAT::Utilities->export_to_level(1, @_); # to be able to load the EXPORT functions when direct call; (normal case)
  AGAT::Utilities->export_to_level(2, @_); # to be able to load the EXPORT functions when called from one level up;
}

=head1 SYNOPSIS

 package containing utility tools

=head1 DESCRIPTION



=head1 AUTHOR

    Jacques Dainat - jacques.dainat@nbis.se

=cut

#check if reference exists in hash. Deep infinite : hash{a} or hash{a}{b} or hash{a}{b}{c}, etc.
# usage example: exists_keys($hash_omniscient,('level3','cds',$level2_ID)
sub exists_keys {
    my ($hash, @keys) = @_;

    for my $key (@keys) {
    	if (ref $hash ne 'HASH' or ! exists $hash->{$key}) {
    		return '';
    	}
		  $hash = $hash->{$key};
    }
    return 1;
}


# @Purpose: check if a value is undef in hash recursively (potentialy due to autovivification)
# @input: 1 =>  hash reference
# @output 1 => bolean
sub exists_undef_value {
    my ($hash) = @_;

    if (ref $hash ne 'HASH'){ print "It is not a hash you provided\n";return ''; }
    my $result = undef;

    foreach my $key (keys %{$hash}) {
      if (ref $hash->{$key} eq 'HASH'){
        if ( exists_undef_value($hash->{$key})){return 1;};
      }
      else{
		    if(! defined $hash->{$key}){
          return 1;
        }
      }
    }
    return '';
}

# @Purpose: check if the table codon is available in bioperl
# @input: 1 =>  integer
# @output 1 => integer
sub get_proper_codon_table {
  my ($codon_table_id_original) = @_;
  my $codonTable = Bio::Tools::CodonTable->new( -id => $codon_table_id_original);
  my $codon_table_id_bioperl = $codonTable->id;

  if ($codon_table_id_original == 0 and  $codon_table_id_original != $codon_table_id_bioperl){
    $codonTable->warn("Your version of bioperl do not handle codon table 0\n".
    "see https://github.com/bioperl/bioperl-live/pull/315\n".
    "It uses codon table $codon_table_id_bioperl instead.");
    return $codon_table_id_bioperl;
  }

  return $codon_table_id_original;
}
1;
