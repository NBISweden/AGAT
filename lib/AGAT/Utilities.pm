#!/usr/bin/perl -w

package AGAT::Utilities;

use strict;
use warnings;
use Time::Piece;
use Time::Seconds;
use Exporter;

our @ISA = qw(Exporter);
our @EXPORT = qw(exists_keys exists_undef_value get_proper_codon_table surround_text
sizedPrint activate_warning_limit print_time dual_print file_text_line print_wrap_text
string_sep_to_hash);

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
# @input: 1 =>  integer, 2 => verbose
# @output 1 => integer
sub get_proper_codon_table {
  my ($codon_table_id_original, $log, $verbose) = @_;
  my $codonTable = Bio::Tools::CodonTable->new( -id => $codon_table_id_original);
  my $codon_table_id_bioperl = $codonTable->id;

  # To deal with empty result in version of bioperl < april 2024 when asking with table 0 (it was reutrning an empty string)
  if (! defined($codon_table_id_bioperl)){
        $codon_table_id_bioperl = 1 ; # default codon table
  }

  if ($codon_table_id_original == 0 and  $codon_table_id_original != $codon_table_id_bioperl){
    $codonTable->warn("Your version of bioperl do not handle codon table 0\n".
    "see https://github.com/bioperl/bioperl-live/pull/315\n".
    "It uses codon table $codon_table_id_bioperl instead.");
  }

  dual_print(
    $log,
    "Codon table ".$codon_table_id_bioperl.
      " in use. You can change it using the appropriate parameter.\n",
    $verbose
  );
  return $codon_table_id_bioperl;
}

# the warning message will be filtered to be printed only $nb_warnings times
# To use it in a script do: my %warnings; activate_warning_limit(\%warnings, $nb_warnings);
# @input: 2 => empty hash, Int
# @output 0 => None
sub activate_warning_limit{
	my ($warnings_hash,$nb_warnings) = @_;

	if (! $nb_warnings){ $nb_warnings = 10;} # Handle to not print to much warning

		$SIG{__WARN__} = sub {
		my $message = shift;

		$warnings_hash->{$message}++;

		if ($warnings_hash->{$message} <= $nb_warnings){
			print "WARNING: ".$message;
		}
		if($warnings_hash->{$message} == $nb_warnings){
			print "************** Too much WARNING of this type we skip the next **************\n";
		}
	};
}

# --------------------------------\\ PRINT //-----------------------------------

# @Purpose: allows to add a frame to a string to print
# @input: 4 =>  String (What has to be printed), int (size of the canevas),
#               Char (character to use to make the frame),
#								String (extra to print at the end after the frame)
# @output 1 => String
#e.g. surround_text("- Start parsing -",80,"*")
#  ********************************************************************************
#  *                              - Start parsing -                               *
#  ********************************************************************************
#
sub surround_text{
  my ($term, $size, $char, $extra) = @_;

   my $frame= $char x $size;
  $frame.="\n";

  my $result = $frame;

  my @lines = split(/\n/,$term);

  	foreach my $line (@lines){

			while ( defined($line) ){
				$line=~ s/^\s+//; #removing leading white spaces
				$line=~ s/\s+$//; #removing trailing white spaces
				$result .= "$char";

	  		my $sizeTerm=length($line)+2; # add 2 extra characters beginning and end of line
		  	if ($sizeTerm > $size ){
			    my $lineout = substr($line, 0,($size - 2),"");
					$lineout .= "$char\n";
					$result .= $lineout;
		 	 	}
		 		else{
			    my $nbBlancBefore=int(($size-$sizeTerm) / 2);
			    my $nbBlancAfter = ($size-$sizeTerm) - $nbBlancBefore;
			    my $lineout =  " " x $nbBlancBefore;
			    $lineout .= $line;
			    $lineout .= " " x $nbBlancAfter;
					$lineout .= "$char\n";
					$result .= $lineout;
					$line  = undef;
		  	}
			}
	}
	$result .= "$frame";
	if($extra){$result .= "$extra";}
	return $result;
}

# e.g.  file_text_line({ string => "parse options and metadata", char => "-" })
# -------------------------- parse options and metadata --------------------------
#
sub file_text_line{
  my ($args) = @_;

  # -------------- OUTPUT --------------
	my $result = "";

	# -------------- INPUT --------------
	# Check we receive a hash as ref
	if(ref($args) ne 'HASH'){ warn "Hash Arguments expected for file_text_line. Please check the call.\n";exit;	}
	# -- Declare all variables and fill them --
	my ($term, $size, $char, $extra);
	# string to print
	if( defined($args->{string})) {$term = $args->{string};} else{ $term=" ";}
	#size line
	if( ! defined($args->{size}) ) { $size = 80;} else{ $size = $args->{size}; }
	# character to fill the line with
	if( ! defined($args->{char}) ) { $char = "-";} else{ $char = $args->{char}; }
	# to add at the end
	if( ! defined($args->{extra}) ) { $extra = undef;} else{ $extra = $args->{extra}; }
	# to add at the beginning
	if( defined($args->{prefix}) ) { $result = $args->{prefix};}

 # -------------- MAIN --------------
  my @lines = split(/\n/,$term);

  	foreach my $line (@lines){
			$line=~ s/^\s+//; #removing leading white spaces
			$line=~ s/\s+$//; #removing trailing white spaces
			$line =" $line ";

			while ( defined($line) ){

	  		my $sizeTerm=length($line);
		  	if ($sizeTerm > $size ){
			    my $lineout = substr($line, 0,($size),"");
					$lineout .= "\n";
					$result .= $lineout;
		 	 	}
		 		else{
			    my $nbCharBefore=int(($size-$sizeTerm) / 2);
			    my $nbCharAfter = ($size-$sizeTerm) - $nbCharBefore;
			    my $lineout =  $char x $nbCharBefore;
			    $lineout .= $line;
			    $lineout .= $char x $nbCharAfter;
					$lineout .= "\n";
					$result .= $lineout;
					$line  = undef;
		  	}
			}
	}
	# add extra
	if($extra){$result .= "$extra";}

	return $result;
}

sub print_wrap_text{
  my ($term, $size, $extra) = @_;

  # -------------- OUTPUT --------------
	my $result = "";
 # -------------- MAIN --------------
  my @lines = split(/\n/,$term);

  	foreach my $line (@lines){
			while ( defined($line) ){

	  		my $sizeTerm=length($line);
		  	if ($sizeTerm > $size ){
			    $result .= substr($line, 0,($size),"");
					$result .= "\n";
		 	 	}
		 		else{
					$result .= "$line\n";
					$line  = undef;
		  	}
			}
	}
	# add extra
	if($extra){$result .= "$extra";}

	return $result;
}

# @Purpose: Print a String in a specified String size. Add space before and after String to
# center it in a decided String Size. If String is longer that Int, we shrink it.
# @input: 2 => String (to be print), Int (size to print the string)
# @output 1 => String
sub sizedPrint{
  my ($term,$size,$extra) = @_;
  my $result;
  my $sizeTerm = (defined($term)) ? length($term) : 0; #defined to deal with string/int 0
  if ($sizeTerm > $size ){
    $result=substr($term, 0,$size);
    return $result;
  }
  else{
    my $nbBlanc=$size-$sizeTerm;

    my $float = $nbBlanc/2;
    my $nbBlanc_before = sprintf "%.0f", $float;
    my $nbBlanc_after = $nbBlanc - $nbBlanc_before;

    $result="";
    for (my $i = 0; $i < $nbBlanc_before; $i++){
      $result.=" ";
    }
    $result.=$term;
    for (my $i = 0; $i < $nbBlanc_after; $i++){
      $result.=" ";
    }
		# add extra
		if($extra){$result .= "$extra";}

    return $result;
  }
}

# add stamotime as suffux before printing
# print like [10:31:23] string
sub print_time{
  my $t = localtime;
  my $line = "[".$t->hms."] @_\n";
  print $line;
}

# @purpose: print to screen when verbose not 0 and always print to log if any log
# @input: 3 => fh, string, integer
# @output 0 => None
sub dual_print{
my ($fh, $string, $min) = @_;
my $verbose = defined $AGAT::AGAT::CONFIG->{verbose} ? $AGAT::AGAT::CONFIG->{verbose} : 1;
$min = 1 unless defined $min;
if ($min > 0 && $verbose >= $min) {
print $string;
}
print $fh $string if $fh;
}

# @Purpose: transform a String with separator into hash
# @input: 2 =>  string, char (the char is the separator)
# @output 1 => hash
sub string_sep_to_hash {
	my $sub_name = (caller(0))[3];
	# -------------- INPUT --------------
	my ($args) = @_;
	# Check we receive a hash as ref
	if(ref($args) ne 'HASH'){ warn "Hash Arguments expected for $sub_name. Please check the call.\n";exit;	}
	# Fill the parameters
	my ($string, $separator);
	if( defined($args->{string})) {$string = $args->{string};} else{ print "String parameter mandatory to use $sub_name!"; exit; }
	if( defined($args->{separator})) {$separator = $args->{separator};} else{ $separator = " ";}

	my %hash_result;
	my @values = split(/$separator/, $string);
   	foreach my $value (@values){
		$hash_result{$value}++;
	}
	return \%hash_result;
}

1;
