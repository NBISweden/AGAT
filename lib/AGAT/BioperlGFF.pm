# JD modified 2022 to deal with features having a single value in the last column (9th)
#
# BioPerl module for Bio::Tools::GFF
#
# Please direct questions and support issues to <bioperl-l@bioperl.org>
#
# Cared for by the Bioperl core team
#
# Copyright Matthew Pocock
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::GFF - A Bio::SeqAnalysisParserI compliant GFF format parser

=head1 SYNOPSIS

    use Bio::Tools::GFF;

    # specify input via -fh or -file
    my $gffio = Bio::Tools::GFF->new(-fh => \*STDIN, -gff_version => 2);
    my $feature;
    # loop over the input stream
    while($feature = $gffio->next_feature()) {
        # do something with feature
    }
    $gffio->close();

    # you can also obtain a GFF parser as a SeqAnalasisParserI in
    # HT analysis pipelines (see Bio::SeqAnalysisParserI and
    # Bio::Factory::SeqAnalysisParserFactory)
    my $factory = Bio::Factory::SeqAnalysisParserFactory->new();
    my $parser = $factory->get_parser(-input => \*STDIN, -method => "gff");
    while($feature = $parser->next_feature()) {
        # do something with feature
    }

=head1 DESCRIPTION

This class provides a simple GFF parser and writer. In the sense of a
SeqAnalysisParser, it parses an input file or stream into SeqFeatureI
objects, but is not in any way specific to a particular analysis
program and the output that program produces.

That is, if you can get your analysis program spit out GFF, here is
your result parser.

=head1 GFF3 AND SEQUENCE DATA

GFF3 supports sequence data; see

http://www.sequenceontology.org/gff3.shtml

There are a number of ways to deal with this -

If you call

  $gffio->ignore_sequence(1)

prior to parsing the sequence data is ignored; this is useful if you
just want the features. It avoids the memory overhead in building and
caching sequences

Alternatively, you can call either

  $gffio->get_seqs()

Or

  $gffio->seq_id_by_h()

At the B<end> of parsing to get either a list or hashref of Bio::Seq
objects (see the documentation for each of these methods)

Note that these objects will not have the features attached - you have
to do this yourself, OR call

  $gffio->features_attached_to_seqs(1)

PRIOR to parsing; this will ensure that the Seqs have the features
attached; ie you will then be able to call

  $seq->get_SeqFeatures();

And use Bio::SeqIO methods

Note that auto-attaching the features to seqs will incur a higher
memory overhead as the features must be cached until the sequence data
is found

=head1 TODO

Make a Bio::SeqIO class specifically for GFF3 with sequence data

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Support

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and
reponsive experts will be able look at the problem and quickly
address it. Please include a thorough description of the problem
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Matthew Pocock

Email mrp-at-sanger.ac.uk

=head1 CONTRIBUTORS

Jason Stajich, jason-at-biperl-dot-org
Chris Mungall, cjm-at-fruitfly-dot-org
Steffen Grossmann [SG], grossman at molgen.mpg.de
Malcolm Cook, mec-at-stowers-institute.org
Jacques Dainat

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package AGAT::BioperlGFF;

use vars qw($HAS_HTML_ENTITIES);
use strict;

use AGAT::OmniscientTool;
use AGAT::AGAT;
use Bio::Seq::SeqFactory;
use Bio::LocatableSeq;
use Bio::SeqFeature::Generic;

use base qw(Bio::Root::Root Bio::SeqAnalysisParserI Bio::Root::IO);

my $i = 0;
my %GFF3_ID_Tags = map { $_ => $i++ } qw(ID Parent Target);

# for skipping data that may be represented elsewhere; currently, this is
# only the score
my %SKIPPED_TAGS = map { $_ => 1 } qw(score);

my %GFF_VERSION = map { $_ => $i++ } qw(1 2 2.5 3);
my %GTF_VERSION = map { $_ => $i++ } qw(1 2 2.1 2.2 2.5 3 relax);

# Set GTF version definitions
my %GTF_FEATURES = ( 
    3 => ["gene", "transcript", "exon", "CDS", "Selenocysteine", "start_codon", "stop_codon", "three_prime_utr", "five_prime_utr"],
    2.5 => ["gene", "transcript", "exon", "CDS", "UTR", "start_codon", "stop_codon", "Selenocysteine"],
    2.2 => ["CDS", "start_codon", "stop_codon", "5UTR", "3UTR", "inter", "inter_CNS", "intron_CNS", "exon"],
    2.1 => ["CDS", "start_codon", "stop_codon", "exon", "5UTR", "3UTR"],
    2 => ["CDS", "start_codon", "stop_codon", "exon"],
    1 => ["CDS", "start_codon", "stop_codon", "exon", "intron"]
);

=head2 new

 Title   : new
 Usage   : my $parser = Bio::Tools::GFF->new(-gff_version => 2,
                                             -file        => "filename.gff");
           or
           my $writer = Bio::Tools::GFF->new(-gff_version => 3,
                                             -file        => ">filename.gff3");
 Function: Creates a new instance. Recognized named parameters are -file, -fh,
           and -gff_version.
 Returns : a new object
 Args    : named parameters
           -gff_version => [1,2,3]

=cut

{   # make a class variable such that we can generate unique ID's over
    # a session, no matter how many instances of GFF.pm we make
    # since these have to be unique within the scope of a GFF file.

    my $gff3_featureID = 0;

    sub _incrementGFF3ID {
        my ($self) = @_;
        return ++ $gff3_featureID;
    }
}


sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);

    my ($gff_version, $version, $type, $noparse) = $self->_rearrange([qw(GFF_VERSION VERSION TYPE NOPARSE)],@args);

    # initialize IO
    $self->_initialize_io(@args);
    # parse header or not
    $self->_parse_header() unless $noparse;
    # check if GFF or GTF
    $self->check_type($type);
    # check version of the GFF or GTF
    $self->check_version($gff_version, $version);
    # ?
    $self->{'_first'} = 1;
    return $self;
}


=head2 _parse_header

 Title   : _parse_header
 Usage   : $gffio->_parse_header()
 Function: used to turn parse GFF header lines.  currently
           produces Bio::LocatableSeq objects from ##sequence-region
           lines
 Returns : 1 on success
 Args    : none


=cut

sub _parse_header{
    my ($self) = @_;

    my @unhandled;
    local $^W = 0; # hide warnings when we try and parse from a file opened
                   # for writing - there isn't really a better way to do
                   # AFAIK - cannot detech if a FH is read or write.
    while(my $line = $self->_readline()){
        my $handled = 0;
        next if /^\s+$/;
        if($line =~ /^\#\#sequence-region\s+(\S+)\s+(\S+)\s+(\S+)\s*/){
            my($seqid,$start,$end) = ($1,$2,$3);
            push @{ $self->{'segments'} }, Bio::LocatableSeq->new(
                -id    => unescape($seqid),
                -start => $start,
                -end   => $end,
                -length => ($end - $start + 1),  ## make the length explicit
           );
           $handled = 1;
        } elsif($line =~ /^(\#\#feature-ontology)/) {
            #to be implemented
            $self->warn("$1 header tag parsing unimplemented");
        } elsif($line =~ /^(\#\#attribute-ontology)/) {
            #to be implemented
            $self->warn("$1 header tag parsing unimplemented");
        } elsif($line =~ /^(\#\#source-ontology)/) {
            #to be implemented
            $self->warn("$1 header tag parsing unimplemented");
        } elsif($line =~ /^(\#\#\#)/) {
            #to be implemented
            $self->warn("$1 header tag parsing unimplemented");
        } elsif($line =~ /^(\#\#FASTA)/) {
            # initial ##FASTA is optional - artemis does not use it
            $line = $self->_readline();
            if ($line !~ /^\>(\S+)/) {
                $self->throw("##FASTA directive must be followed by fasta header, not: $line");
            }
        }

        if ($line =~ /^\>(.*)/) {
            # seq data can be at header or footer
            my $seq = $self->_parse_sequence($line);
            if ($seq) {
                $self->_seq_by_id_h->{$seq->primary_id} = $seq;
            }
        }


        if(!$handled){
            push @unhandled, $line;
        }

        #looks like the header is over!
        last unless $line =~ /^\#/;
    }

    foreach my $line (@unhandled){
        $self->_pushback($line);
    }

    return 1;
}

sub _parse_sequence {
    my ($self, $line) = @_;

    if ($line =~ /^\>(.*)/) {

        my $seqid = $1;
        $seqid =~ s/\s+$//;
        my $desc = '';
        if ($seqid =~ /(\S+)\s+(.*)/) {
            ($seqid, $desc) = ($1,$2);
        }
        my $res = '';
        while (my $line = $self->_readline) {
            if ($line =~ /^\#/) {
                last;
            }
            if ($line =~ /^\>/) {
                $self->_pushback($line);
                last;
            }
            $line =~ s/\s//g;
            $res .= $line;
        }
        return if $self->ignore_sequence;

        my $seqfactory = Bio::Seq::SeqFactory->new('Bio::Seq');
        my $seq = $seqfactory->create(-seq=>$res,
                                      -id=>$seqid,
                                      -desc=>$desc);
        $seq->accession_number($seqid);
        if ($self->features_attached_to_seqs) {
            my @feats =
              @{$self->_feature_idx_by_seq_id->{$seqid}};
            $seq->add_SeqFeature($_) foreach @feats;
            @{$self->_feature_idx_by_seq_id->{$seqid}} = ();
        }
        return $seq;
    }
    else {
        $self->throw("expected fasta header, not: $line");
    }
}


=head2 next_segment

 Title   : next_segment
 Usage   : my $seq = $gffio->next_segment;
 Function: Returns a Bio::LocatableSeq object corresponding to a
           GFF "##sequence-region" header line.
 Example :
 Returns : A Bio::LocatableSeq object, or undef if
           there are no more sequences.
 Args    : none


=cut

sub next_segment{
    my ($self,@args) = @_;
    return shift @{ $self->{'segments'} } if defined $self->{'segments'};
    return;
}


=head2 next_feature

 Title   : next_feature
 Usage   : $seqfeature = $gffio->next_feature();
 Function: Returns the next feature available in the input file or stream, or
           undef if there are no more features.
 Example :
 Returns : A Bio::SeqFeatureI implementing object, or undef if there are no
           more features.
 Args    : none

=cut

sub next_feature {
    my ($self) = @_;

    my $gff_string;

    # be graceful about empty lines or comments, and make sure we return undef
    # if the input's consumed
    while(($gff_string = $self->_readline()) && defined($gff_string)) {
        if ($gff_string =~ /^\#\#\#/) {
            # all forward refs have been seen; TODO
        }
        next if($gff_string =~ /^\#/ || $gff_string =~ /^\s*$/ ||
                $gff_string =~ m{^//});

        while ($gff_string =~ /^\>(.+)/) {
            # fasta can be in header or footer
            my $seq = $self->_parse_sequence($gff_string);
            if ($seq) {
                $self->_seq_by_id_h->{$seq->primary_id} = $seq;
                $gff_string = $self->_readline;
                last unless $gff_string;
            }
        }
        last;
    }
    return unless $gff_string;

    my $feat = Bio::SeqFeature::Generic->new();
    $self->from_gff_string($feat, $gff_string);

    if ($self->features_attached_to_seqs) {
        push(@{$self->_feature_idx_by_seq_id->{$feat->seq_id}},
             $feat);
    }

    return $feat;
}

sub _feature_idx_by_seq_id {
    my $self = shift;
    $self->{__feature_idx_by_seq_id} = shift if @_;
    $self->{__feature_idx_by_seq_id} = {}
      unless $self->{__feature_idx_by_seq_id};
    return $self->{__feature_idx_by_seq_id};
}


=head2 from_gff_string

 Title   : from_gff_string
 Usage   : $gff->from_gff_string($feature, $gff_string);
 Function: Sets properties of a SeqFeatureI object from a GFF-formatted
           string. Interpretation of the string depends on the version
           that has been specified at initialization.

           This method is used by next_feature(). It actually dispatches to
           one of the version-specific (private) methods.
 Example :
 Returns : void
 Args    : A Bio::SeqFeatureI implementing object to be initialized
           The GFF-formatted string to initialize it from

=cut

sub from_gff_string {
    my ($self, $feat, $gff_string) = @_;

    if ($self->{'TYPE'} eq "GFF"){
        if($self->{'VERSION'} == 1)  {
            return $self->_from_gff1_string($feat, $gff_string);
        } elsif( $self->{'VERSION'} == 3 ) {
            return $self->_from_gff3_string($feat, $gff_string);
        } else { # GFF2 and 2.5
            return $self->_from_gff2_string($feat, $gff_string);
        }
    } else { # all GTF cases
        return $self->_from_gff2_string($feat, $gff_string);
    }
}


=head2 _from_gff1_string

 Title   : _from_gff1_string
 Usage   :
 Function:
 Example :
 Returns : void
 Args    : A Bio::SeqFeatureI implementing object to be initialized
           The GFF-formatted string to initialize it from

=cut

sub _from_gff1_string {
    my ($gff, $feat, $string) = @_;
    chomp $string;
    my ($seqname, $source, $primary, $start, $end, $score,
        $strand, $frame, @group) = split(/\t/, $string);

    if ( !defined $frame ) {
        $feat->throw("[$string] does not look like GFF to me");
    }
    $frame = 0 unless( $frame =~ /^\d+$/);
    $seqname = 'SEQ' if ! length($seqname);
    $feat->seq_id($seqname);
    $feat->source_tag($source);
    $feat->primary_tag($primary);
    $feat->start($start);
    $feat->end($end);
    $feat->frame($frame);
    if ( $score eq '.' ) {
        #$feat->score(undef);
    } else {
        $feat->score($score);
    }
    if ( $strand eq '-' ) { $feat->strand(-1); }
    if ( $strand eq '+' ) { $feat->strand(1); }
    if ( $strand eq '.' ) { $feat->strand(0); }
    foreach my $g ( @group ) {
        if ( $g =~ /(\S+)=(\S+)/ ) {
            my $tag = $1;
            my $value = $2;
            $feat->add_tag_value($1, $2);
        } else {
            $feat->add_tag_value('group', $g);
        }
    }
}


=head2 _from_gff2_string

 Title   : _from_gff2_string
 Usage   :
 Function:
 Example :
 Returns : void
 Args    : A Bio::SeqFeatureI implementing object to be initialized
           The GFF2-formatted string to initialize it from


=cut

sub _from_gff2_string {
    my ($gff, $feat, $string) = @_;
    chomp($string);

    # according to the Sanger website, GFF2 should be single-tab
    # separated elements, and the free-text at the end should contain
    # text-translated tab symbols but no "real" tabs, so splitting on
    # \t is safe, and $attribs gets the entire attributes field to be
    # parsed later

    # sendu: but the tag value pair can (should?) be separated by a tab. The
    # 'no tabs' thing seems to apply only to the free text that is allowed for
    # the value

    my ($seqname, $source, $primary, $start,
        $end, $score, $strand, $frame, @attribs) = split(/\t+/, $string);
    my $attribs = join ' ', @attribs;

    if ( !defined $frame ) {
        $feat->throw("[$string] does not look like GFF2 to me");
    }
    $seqname = 'SEQ' if ! length($seqname);
    $feat->seq_id($seqname);
    $feat->source_tag($source);
    $feat->primary_tag($primary);
    $feat->start($start);
    $feat->end($end);
    $feat->frame($frame);
    if ( $score eq '.' ) {
        # $feat->score(undef);
    } else {
        $feat->score($score);
    }
    if ( $strand eq '-' ) { $feat->strand(-1); }
    if ( $strand eq '+' ) { $feat->strand(1); }
    if ( $strand eq '.' ) { $feat->strand(0); }


    #  <Begin Inefficient Code from Mark Wilkinson>
    # this routine is necessay to allow the presence of semicolons in
    # quoted text Semicolons are the delimiting character for new
    # tag/value attributes.  it is more or less a "state" machine, with
    # the "quoted" flag going up and down as we pass thorugh quotes to
    # distinguish free-text semicolon and hash symbols from GFF control
    # characters

    my $flag = 0; # this could be changed to a bit and just be twiddled
    my @parsed;

    # run through each character one at a time and check it
    # NOTE: changed to foreach loop which is more efficient in perl
    # --jasons
    for my $a ( split //, $attribs ) {
        # flag up on entering quoted text, down on leaving it
        if( $a eq '"') { $flag = ( $flag == 0 ) ? 1:0 }
        elsif( $a eq ';' && $flag ) { $a = "INSERT_SEMICOLON_HERE"}
        elsif( $a eq '#' && ! $flag ) { last }
        push @parsed, $a;
    }
    $attribs = join "", @parsed; # rejoin into a single string

    # <End Inefficient Code>
    # Please feel free to fix this and make it more "perlish"

    my @key_vals = split /;/, $attribs;   # attributes are semicolon-delimited

    foreach my $pair ( @key_vals ) {
        # replace semicolons that were removed from free-text above.
        $pair =~ s/INSERT_SEMICOLON_HERE/;/g;
        # separate the key from the value
        my ($blank, $key, $values) = split  /^\s*([\w\d]+)\s/, $pair;
        if( defined $values ) {
            my @values;
            # free text is quoted, so match each free-text block
            # and remove it from the $values string
            while ($values =~ s/"(.*?)"//){
                # and push it on to the list of values (tags may have
                # more than one value... and the value may be undef)
                push @values, $1;
            }

            # and what is left over should be space-separated
            # non-free-text values

            my @othervals = split /\s+/, $values;
            foreach my $othervalue(@othervals){
                # get rid of any empty strings which might
                # result from the split
                if (CORE::length($othervalue) > 0) {push @values, $othervalue}
            }

            foreach my $value(@values){
                $feat->add_tag_value($key, $value);
            }
        }
        elsif( $blank ){
           $feat->add_tag_value("ID", $blank );
        }
    }
}


sub _from_gff3_string {
    my ($gff, $feat, $string) = @_;
    chomp($string);

    # according to the now nearly final GFF3 spec, columns should
    # be tab separated, allowing unescaped spaces to occur in
    # column 9

    my ($seqname, $source, $primary, $start, $end,
        $score, $strand, $frame, $attribs) = split(/\t/, $string);

    if ( ! defined $frame ) {
        $feat->throw("[$string] does not look like GFF3 to me");
    }
    $seqname = 'SEQ' if ! length($seqname);
    $feat->seq_id($seqname);
    $feat->source_tag($source);
    $feat->primary_tag($primary);
    $feat->start($start);
    $feat->end($end);
    $feat->frame($frame);
    if ( $score eq '.' ) {
        #$feat->score(undef);
    } else {
        $feat->score($score);
    }
    if ( $strand eq '-' ) { $feat->strand(-1); }
    if ( $strand eq '+' ) { $feat->strand(1); }
    if ( $strand eq '.' ) { $feat->strand(0); }
    my @attribs = split(/\s*;\s*/, $attribs);

    my $size_attribs = scalar(@attribs);
    my $attrib_nb=0;
    my $tag_ID_seen=0;
    my $potential_ID;
    for my $attrib (@attribs) {
        next if( ! $attrib); # avoid issue 528 when two semicolons in a row it will create empty value "" - ID=blabla;;Name=blabla => ID=blabla,"";Name=blabla
        $attrib_nb++;

        my ($tag,$value) = split /=/,$attrib;
        my @values       = map {$_} split /,/,$value;
        $tag_ID_seen++ if ( lc($tag) eq "id" );
        # Case where attribute column contain only one value, no tag/attribute structure 
        # e.g. augustus/tsebra, use the value as ID
        if ( scalar(@values) == 0 ) {
            if ( $size_attribs == 1 ){
                warn "No attribute tag available. A single value provided in 9th column. Using it as ID: ID=<value> (augustus/tsebra) @ - $tag";
                $value = $tag;
                $feat->add_tag_value("ID",$value);
            }
            elsif ( $attrib_nb == 1 ){ 
                warn "No attribute tag found. Since other attributes exist and this is the first one, it is assumed to be the ID (ID=<value>) when the ID attribute is missing (augustus/tsebra). @ Value: $tag";
                $potential_ID = $tag;
            }
        }
        # other cases
        else{
            for my $v ( @values ) { $feat->add_tag_value($tag,$v); }
        }
    }
     # Add (augustus/tsebra) L2 is missing ID attribute, so we create one
    if ($potential_ID && !$tag_ID_seen){
        $feat->add_tag_value("ID",$potential_ID);
    }
}

# taken from Bio::DB::GFF
sub unescape {
  my $v = shift;
  $v =~ s/%([0-9a-fA-F]{2})/chr hex($1)/ge;
  return $v;
}


sub escape{
    my $v = shift;
    # Only decode if url_decode_in option is enabled
    if ($AGAT::AGAT::CONFIG->{url_encode_out}) {
       $v = escape_gff3($v);
    } else {
       $v = unescape_gff3($v);
    }
    return $v;
}

# taken from Bio::DB::GFF
sub unescape_gff3 {
  my $v = shift;
  # Only decode tab (09), comma (2C), semicolon (3B), and equals (3D)
  $v =~ s/%(09|2C|3B|3D)/chr hex($1)/gie;
  return $v;
}

# escape only tab,=; characters
sub escape_gff3 {
  my $v = shift;
  $v =~ s/([\t,=;])/sprintf("%%%X",ord($1))/ge;
  return $v;
}

=head2 write_feature

 Title   : write_feature
 Usage   : $gffio->write_feature($feature);
 Function: Writes the specified SeqFeatureI object in GFF format to the stream
           associated with this instance.
 Returns : none
 Args    : An array of Bio::SeqFeatureI implementing objects to be serialized

=cut

sub write_feature {
    my ($self, @features) = @_;
    return unless @features;
    if( $self->{'_first'} ) {
        my $line = "##";
        if ( $self->{'TYPE'} eq "GTF"){
            $line .= "gtf-version ";
        } else {
            $line .= "gff-version ";
        }
        $line .= $self->{'VERSION'}."\n";
        $self->_print($line);
    }
    $self->{'_first'} = 0;

    # deflate multi-values attribute if asked by config
    if ($AGAT::AGAT::CONFIG->{'deflate_attribute'}){
        foreach my $feature ( @features ) {
            my @list_tags= $feature->get_all_tags();
			foreach my $tag (@list_tags){
				my @tag_values = $feature->get_tag_values($tag);
                if ($#tag_values >= 1){
                    my $tag_counter=-1;
					foreach my $tag_value (@tag_values){
                        $tag_counter++;
                        if ($tag_counter == 0){
                            create_or_replace_tag($feature, $tag , $tag_value);
                        } else {
                            create_or_replace_tag($feature, $tag."_".$tag_counter , $tag_value);
                        }	
					}
				}
			}
        }
    }

    foreach my $feature ( @features ) {
        $self->_print($self->gxf_string($feature)."\n");
    }
}


=head2 gxf_string

 Title   : gff_string
 Usage   : $gffstr = $gffio->gxf_string($feature);
 Function: Obtain the GFF-formatted representation of a SeqFeatureI object.
           The formatting depends on the version specified at initialization.

           This method is used by write_feature(). It actually dispatches to
           one of the version-specific (private) methods.
 Example :
 Returns : A GFF-formatted string representation of the SeqFeature
 Args    : A Bio::SeqFeatureI implementing object to be GFF-stringified

=cut

sub gxf_string{
    my ($self, $feature) = @_;

    if ($self->{'TYPE'} eq "GFF"){
        my $version = $self->check_version();
        if( $version == 1) {
            return $self->_gff1_string($feature);
        } elsif( $version == 3 ) {
            return $self->_gff3_string($feature);
        } elsif( $version == 2.5 ) {
            return $self->_gff25_string($feature);
        } else {
            return $self->_gff2_string($feature);
        }
    } else { # GTF case
        return $self->_gff25_string($feature);
    }
}


=head2 _gff1_string

 Title   : _gff1_string
 Usage   : $gffstr = $gffio->_gff1_string
 Function:
 Example :
 Returns : A GFF1-formatted string representation of the SeqFeature
 Args    : A Bio::SeqFeatureI implementing object to be GFF-stringified

=cut

sub _gff1_string{
    my ($gff, $feat) = @_;
    my ($str,$score,$frame,$name,$strand);

    if( $feat->can('score') ) {
        $score = $feat->score();
    }
    $score = '.' unless defined $score;

    if( $feat->can('frame') ) {
        $frame = $feat->frame();
    }
    $frame = '.' unless defined $frame;

    $strand = $feat->strand();
    if(! $strand) {
        $strand = ".";
    } elsif( $strand == 1 ) {
        $strand = '+';
    } elsif ( $feat->strand == -1 ) {
        $strand = '-';
    }

    if( $feat->can('seqname') ) {
        $name = $feat->seq_id();
    }
    $name = 'SEQ' if ! length($name);

    $str = join("\t",
                $name,
                $feat->source_tag,
                $feat->primary_tag,
                $feat->start,
                $feat->end,
                $score,
                $strand,
                $frame);

    foreach my $tag ( $feat->get_all_tags ) {
        next if exists $SKIPPED_TAGS{$tag};
        foreach my $value ( $feat->get_tag_values($tag) ) {
            $value=escape($value);
            $str .= " $tag=$value" if $value;
        }
    }

    return $str;
}


=head2 _gff2_string

 Title   : _gff2_string
 Usage   : $gffstr = $gffio->_gff2_string
 Function:
 Example :
 Returns : A GFF2-formatted string representation of the SeqFeature
 Args    : A Bio::SeqFeatureI implementing object to be GFF2-stringified

=cut

sub _gff2_string{
    my ($gff, $origfeat) = @_;
    my $feat;
    if ($origfeat->isa('Bio::SeqFeature::FeaturePair')){
        $feat = $origfeat->feature2;
    } else {
        $feat = $origfeat;
    }
    my ($str1, $str2,$score,$frame,$name,$strand);

    if( $feat->can('score') ) {
        $score = $feat->score();
    }
    $score = '.' unless defined $score;

    if( $feat->can('frame') ) {
        $frame = $feat->frame();
    }
    $frame = '.' unless defined $frame;

    $strand = $feat->strand();
    if(! $strand) {
        $strand = ".";
    } elsif( $strand == 1 ) {
        $strand = '+';
    } elsif ( $feat->strand == -1 ) {
        $strand = '-';
    }

    if( $feat->can('seqname') ) {
        $name = $feat->seq_id();
    }
    $name = 'SEQ' if ! length($name);

    $str1 = join("\t",
                 $name,
                 $feat->source_tag(),
                 $feat->primary_tag(),
                 $feat->start(),
                 $feat->end(),
                 $score,
                 $strand,
                 $frame);
    # the routine below is the only modification I made to the original
    # ->gff_string routine (above) as on November 17th, 2000, the
    # Sanger webpage describing GFF2 format reads: "From version 2
    # onwards, the attribute field must have a tag value structure
    # following the syntax used within objects in a .ace file,
    # flattened onto one line by semicolon separators. Tags must be
    # standard identifiers ([A-Za-z][A-Za-z0-9_]*).  Free text values
    # must be quoted with double quotes".
    # MW

    my @group;

    foreach my $tag ( $feat->get_all_tags ) {
        next if exists $SKIPPED_TAGS{$tag};
        my @v;
        foreach my $value ( $feat->get_tag_values($tag) ) {
            unless( defined $value && length($value) ) {
                # quote anything other than valid tag/value characters
                $value = '""';
            } elsif ($value =~ /[^A-Za-z0-9_]/){
                # substitute tab and newline chars by their UNIX equivalent
                $value =~ s/\t/\\t/g;
                $value =~ s/\n/\\n/g;
                $value = '"' . $value . '" ';
            }
            push @v, escape($value);
            # for this tag (allowed in GFF2 and .ace format)
        }
        push @group, "$tag ".join(" ", @v);
    }

    $str2 .= join(' ; ', @group);
    # Add Target information for Feature Pairs
    if( ! $feat->has_tag('Target') && # This is a bad hack IMHO
        ! $feat->has_tag('Group')  &&
        $origfeat->isa('Bio::SeqFeature::FeaturePair') ) {
        $str2 = sprintf("Target %s %d %d", $origfeat->feature1->seq_id,
                       ( $origfeat->feature1->strand < 0 ?
                         ( $origfeat->feature1->end,
                           $origfeat->feature1->start) :
                         ( $origfeat->feature1->start,
                           $origfeat->feature1->end)
                         )) . ($str2?" ; ".$str2:"");  # need to put Target information before other tag/value pairs - mw
    }
    return $str1."\t".$str2;
}

=head2 _gff25_string

 Title   : _gff25_string
 Usage   : $gffstr = $gffio->_gff2_string
 Function: To get a format of GFF that is peculiar to Gbrowse/Bio::DB::GFF
 Example : 9th column: ID "gene-1"; Name "name 1" name2;
 Returns : A GFF2.5-formatted string representation of the SeqFeature
 Args    : A Bio::SeqFeatureI implementing object to be GFF2.5-stringified
 Comments: GFF2.5 is suposed to be similar as GTF (with semicolon at the end).
=cut

sub _gff25_string {
    my ($gff, $origfeat) = @_;

    # for skipping data that may be represented elsewhere; currently, this is
		# only the score
		my %SKIPPED_TAGS = map { $_ => 1 } qw(score);

    my $feat;
    if ($origfeat->isa('Bio::SeqFeature::FeaturePair')){
        $feat = $origfeat->feature2;
    } else {
        $feat = $origfeat;
    }

    if ($gff->{'TYPE'} eq "GTF" and $gff->{'VERSION'} ne "relax"){
        if(!  grep {lc($feat->primary_tag()) eq lc($_)} ( @{$GTF_FEATURES {$gff->{'VERSION'}}} )) {
            next;
        }
    }

    my ($str1, $str2,$score,$frame,$name,$strand);

    if( $feat->can('score') ) {
        $score = $feat->score();
    }
    $score = '.' unless defined $score;

    if( $feat->can('frame') ) {
        $frame = $feat->frame();
    }
    $frame = '.' unless defined $frame;

    $strand = $feat->strand();
    if(! $strand) {
        $strand = ".";
    } elsif( $strand == 1 ) {
        $strand = '+';
    } elsif ( $feat->strand == -1 ) {
        $strand = '-';
    }

    if( $feat->can('seqname') ) {
        $name = $feat->seq_id();
    }
    $name = 'SEQ' if ! length($name);

    $str1 = join("\t",
                 $name,
                 $feat->source_tag(),
                 $feat->primary_tag(),
                 $feat->start(),
                 $feat->end(),
                 $score,
                 $strand,
                 $frame);

    my @all_tags = $feat->all_tags;
    my @group; my @firstgroup;

    if (@all_tags) {   # only play this game if it is worth playing...
        foreach my $tag ( @all_tags ) {
            next if exists $SKIPPED_TAGS{$tag};
            my @v;
            foreach my $value ( $feat->get_tag_values($tag) ) {
                unless( defined $value && length($value) ) {
                    $value = '""';
                } else{ # quote all type of values
                    $value =~ s/\t/\\t/g; # substitute tab and newline
                    # characters
                    $value =~ s/\n/\\n/g; # to their UNIX equivalents
                    $value = '"' . $value . '"';
                }
                push @v, escape($value);
            }
            $v[$#v] =~ s/\s+$//; #remove left space of the last value
            if (($tag eq 'gene_id') || ($tag eq 'transcript_id')){ # hopefully we won't get both...
                push @firstgroup, "$tag ".join(" ", @v);
            } else {
                push @group, "$tag ".join(" ", @v);
            }
        }
    }
     @firstgroup = sort @firstgroup if @firstgroup;
    $str2 = join('; ', (@firstgroup, @group));
    $str2 .= ";";
    # Add Target information for Feature Pairs
    if( ! $feat->has_tag('Target') && # This is a bad hack IMHO
        ! $feat->has_tag('Group') &&
        $origfeat->isa('Bio::SeqFeature::FeaturePair') ) {
        $str2 = sprintf("Target %s ; tstart %d ; tend %d", $origfeat->feature1->seq_id,
                        ( $origfeat->feature1->strand < 0 ?
                          ( $origfeat->feature1->end,
                            $origfeat->feature1->start) :
                          ( $origfeat->feature1->start,
                            $origfeat->feature1->end)
                        )) . ($str2?" ; ".$str2:""); # need to put the target info before other tag/value pairs - mw
    }
    return $str1 . "\t".  $str2;
}

=head2 _gff3_string

  Title   : _gff3_string
  Usage   : $gffstr = $gffio->_gff3_string
  Function:
  Example :
  Returns : A GFF3-formatted string representation of the SeqFeature
  Args    : A Bio::SeqFeatureI implementing object to be GFF3-stringified

=cut

sub _gff3_string {
    my ($gff, $origfeat) = @_;
    my $feat;
    if ($origfeat->isa('Bio::SeqFeature::FeaturePair')){
        $feat = $origfeat->feature2;
    } else {
        $feat = $origfeat;
    }

    my $ID = $gff->_incrementGFF3ID();

    my ($score,$frame,$name,$strand);

    if( $feat->can('score') ) {
        $score = $feat->score();
    }
    $score = '.' unless defined $score;

    if( $feat->can('frame') ) {
        $frame = $feat->frame();
    }
    $frame = '1' unless defined $frame;

    $strand = $feat->strand();

    if(! $strand) {
        $strand = ".";
    } elsif( $strand == 1 ) {
        $strand = '+';
    } elsif ( $feat->strand == -1 ) {
        $strand = '-';
    }

    if( $feat->can('seqname') ) {
        $name = $feat->seq_id();
    } 
    $name = 'SEQ' if ! length($name);
    
    my @groups;

    # force leading ID and Parent tags
    my @all_tags =  grep { ! exists $GFF3_ID_Tags{$_} } $feat->all_tags;
    for my $t ( sort { $GFF3_ID_Tags{$b} <=> $GFF3_ID_Tags{$a} }
                keys %GFF3_ID_Tags ) {
        unshift @all_tags, $t if $feat->has_tag($t);
    }

    for my $tag ( @all_tags ) {
    next if exists $SKIPPED_TAGS{$tag};

        my @values = $feat->get_tag_values($tag);
        # next if $tag eq 'Target';
        if ($tag eq 'Target' && ! $origfeat->isa('Bio::SeqFeature::FeaturePair')){      
            if(scalar(@values) > 1){ # How is it possible that Target is has a value list ??
                # simple Target,start,stop
                my ($target_id, $b,$e,$strand) = @values;
                next unless(defined($e) && defined($b) && $target_id);
                ($b,$e)= ($e,$b) if(defined $strand && $strand<0);
                #if we have the strand we will print it
                if($strand){ push @groups, sprintf("Target=%s %d %d %s", $target_id,$b,$e,$strand); }
                else{ push @groups, sprintf("Target=%s %d %d", $target_id,$b,$e); }
                next;
            }
        }

        my $valuestr;
        # a string which will hold one or more values
        # for this tag, with quoted free text and
        # space-separated individual values.
        my @v;
        for my $value ( @values ) {
            if(  defined $value && length($value) ) {
                
                # Multiple tag=value pairs are separated by semicolons. 
                # URL escaping rules are used for tags or values containing the following characters: ",=;". 
                # Spaces are allowed in this field, but tabs must be replaced with the %09 URL escape. 
                # Attribute values do not need to be and should not be quoted. 
                #The quotes should be included as part of the value by parsers and not stripped.
                if ( ref $value eq 'Bio::Annotation::Comment') {
                    $value = $value->text;
                }

                #if ($value =~ /[^a-zA-Z0-9\,\;\=\.:\%\^\*\$\@\!\+\_\?\-]/) {
                #    $value =~ s/\t/\\t/g; # substitute tab and newline
                #    # characters
                #    $value =~ s/\n/\\n/g; # to their UNIX equivalents

                    # Unescaped quotes are not allowed in GFF3
                    #                    $value = '"' . $value . '"';
                #}
                $value = escape($value);
            } else {
                # if it is completely empty, then just make empty double quotes
                $value = '""';
            }
            push @v, $value;
        }
        # can we figure out how to improve this?
        $tag = lcfirst($tag) unless ( $tag =~
            /^(ID|Name|Alias|Parent|Gap|Target|Derives_from|Note|Dbxref|Ontology_term)$/);

        push @groups, "$tag=".join(",",@v);
    }
    # Add Target information for Feature Pairs
    if( $feat->has_tag('Target') &&
        ! $feat->has_tag('Group') &&
        $origfeat->isa('Bio::SeqFeature::FeaturePair') ) {

        my $target_id = $origfeat->feature1->seq_id;
        $target_id =~ s/([\t\n\r%&\=;,])/sprintf("%%%X",ord($1))/ge;

        push @groups, sprintf("Target=%s %d %d",
                              $target_id,
                              ( $origfeat->feature1->strand < 0 ?
                                ( $origfeat->feature1->end,
                                  $origfeat->feature1->start) :
                                ( $origfeat->feature1->start,
                                  $origfeat->feature1->end)
                                ));
    }

    # unshift @groups, "ID=autogenerated$ID" unless ($feat->has_tag('ID'));
    if ( $feat->can('name') && defined($feat->name) ) {
        # such as might be for Bio::DB::SeqFeature
        unshift @groups, 'Name=' . $feat->name;
    }

    my $gff_string = "";
    my $source = $feat->source_tag() || '.';
    my $primary = $feat->primary_tag();
    if ($feat->location->isa("Bio::Location::SplitLocationI")) {
        my @locs = $feat->location->each_Location;
        foreach my $loc (@locs) {
            $gff_string .= join("\t",
                                $name,
                                $source,
                                $primary,
                                $loc->start(),
                                $loc->end(),
                                $score,
                                $strand,
                                $frame,
                                join(';', @groups)) . "\n";
        }
        chop $gff_string;
        return $gff_string;
    } else {
        $gff_string = join("\t",
                           $name,
                           $source,
                           $primary,
                           $feat->start(),
                           $feat->end(),
                           $score,
                           $strand,
                           $frame,
                           join(';', @groups));
    }
    return $gff_string;
}


=head2 check_version

  Title   : check_version
  Usage   : $version = $gffio->check_version
  Function:
  Example :
  Returns : The GFF version this parser will accept and emit.
  Args    : none
  Comment : priority of VERSION over GFF_VERSION (GFF_VERSION is kept for retro-compatibility)


=cut

sub check_version {
    my ($self, $gff_version, $version) = @_;

    # GFF case
    if ($self->{'TYPE'} eq "GFF"){
        if (! $version){
            if ($gff_version){
                $version = $gff_version;
            }
            else{
                $version = 3;
            }
        }
        if( exists_keys(\%GFF_VERSION, ($version) ) ) {
            $self->{'GFF_VERSION'} = $version;
            $self->{'VERSION'} = $version;
        }
        else  {
            $self->throw("Can't build a GFF object with the unknown version ". $version  .". Accepted values are: 1, 2, 2.5, 3");
        }
    }
    elsif ($self->{'TYPE'} eq "GTF"){
        if (! $version){
            if ($gff_version){
                $version = $gff_version;
            }
            else{
                $version = "relax";
            }
        }
        if( exists_keys(\%GTF_VERSION, ($version) ) ) {
            $self->{'GFF_VERSION'} = $version;
            $self->{'VERSION'} = $version;
        }
        else  {
            $self->throw("Can't build a GTF object with the unknown version ". $version . ". Accepted values are: 1, 2, 2.1, 2.2, 2.5, 3, relax");
        }
    }
    else {
        $self->throw("Can't build GFF/GTF object of the unknown type ".
            $$self->{'TYPE'} . ". Accepted value is GFF and GTF. If none provided GFF is used as default.");
    }
    return $self->{'VERSION'};
}

=head2 check_type

  Title   : check_type
  Usage   : $type = $gffio->check_type
  Function:
  Example :
  Returns : The file type GFF or GTF this parser will accept and emit.
  Args    : none

=cut
sub check_type {
    my ($self, $type) = @_;
    if( $type ) {
        if (lc($type) eq "gff"){
            $self->{'TYPE'} = "GFF";
        } elsif (lc($type) eq "gtf") {
            $self->{'TYPE'} = "GTF";
        } else {
            $self->throw("Can't build GFF/GTF object of the unknown type ".
            $type . ". Accepted value is GFF and GTF. If none provided GFF is used as default.");
        }   
    } else { # default is GFF
        $self->{'TYPE'} = "GFF";
    }
    return $self->{'TYPE'};
}

# Make filehandles

=head2 newFh

 Title   : newFh
 Usage   : $fh = Bio::Tools::GFF->newFh(-file=>$filename,-format=>'Format')
 Function: does a new() followed by an fh()
 Example : $fh = Bio::Tools::GFF->newFh(-file=>$filename,-format=>'Format')
           $feature = <$fh>;            # read a feature object
           print $fh $feature;          # write a feature object
 Returns : filehandle tied to the Bio::Tools::GFF class
 Args    :

=cut

sub newFh {
    my $class = shift;
    return unless my $self = $class->new(@_);
    return $self->fh;
}


=head2 fh

 Title   : fh
 Usage   : $obj->fh
 Function:
 Example : $fh = $obj->fh;      # make a tied filehandle
           $feature = <$fh>;    # read a feature object
           print $fh $feature;  # write a feature object
 Returns : filehandle tied to Bio::Tools::GFF class
 Args    : none

=cut


sub fh {
    my $self = shift;
    my $class = ref($self) || $self;
    my $s = Symbol::gensym;
    tie $$s,$class,$self;
    return $s;
}

# This accessor is used for accessing the Bio::Seq objects from a GFF3
# file; if the file you are using has no sequence data you can ignore
# this accessor

# This accessor returns a hash reference containing Bio::Seq objects,
# indexed by Bio::Seq->primary_id

sub _seq_by_id_h {
    my $self = shift;

    return $self->{'_seq_by_id_h'} = shift if @_;
    $self->{'_seq_by_id_h'} = {}
    unless $self->{'_seq_by_id_h'};
    return $self->{'_seq_by_id_h'};
}


=head2 get_seqs

 Title   : get_seqs
 Usage   :
 Function: Returns all Bio::Seq objects populated by GFF3 file
 Example :
 Returns :
 Args    :

=cut

sub get_seqs {
    my ($self,@args) = @_;
    return values %{$self->_seq_by_id_h};
}


=head2 features_attached_to_seqs

 Title   : features_attached_to_seqs
 Usage   : $obj->features_attached_to_seqs(1);
 Function: For use with GFF3 containing sequence only

    Setting this B<before> parsing ensures that all Bio::Seq object
    created will have the appropriate features added to them

    defaults to false (off)

    Note that this mode will incur higher memory usage because features
    will have to be cached until the relevant feature comes along

 Example :
 Returns : value of features_attached_to_seqs (a boolean)
 Args    : on set, new value (a boolean, optional)


=cut

sub features_attached_to_seqs{
    my $self = shift;

    return $self->{'_features_attached_to_seqs'} = shift if @_;
    return $self->{'_features_attached_to_seqs'};
}


=head2 ignore_sequence

 Title   : ignore_sequence
 Usage   : $obj->ignore_sequence(1);
 Function: For use with GFF3 containing sequence only

    Setting this B<before> parsing means that all sequence data will be
    ignored

 Example :
 Returns : value of ignore_sequence (a boolean)
 Args    : on set, new value (a boolean, optional)

=cut

sub ignore_sequence{
    my $self = shift;

    return $self->{'_ignore_sequence'} = shift if @_;
    return $self->{'_ignore_sequence'};
}


sub DESTROY {
    my $self = shift;
    $self->close();
}

sub TIEHANDLE {
    my ($class,$val) = @_;
    return bless {'gffio' => $val}, $class;
}

sub READLINE {
    my $self = shift;
    return $self->{'gffio'}->next_feature() || undef unless wantarray;
    my (@list, $obj);
    push @list, $obj while $obj = $self->{'gffio'}->next_feature();
    return @list;
}

sub PRINT {
    my $self = shift;
    $self->{'gffio'}->write_feature(@_);
}

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

1;