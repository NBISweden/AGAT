<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 3.2//EN">
<html>
  Copy from <a href="https://web.archive.org/web/20010208224442/http://www.sanger.ac.uk/Software/formats/GFF/GFF_Spec.shtml">web.archive.org</a>.

<h1>The Sanger Centre : GFF2 </h1>
    
    
<table width="100%" cellspacing="0" cellpadding="0">
  <tr class="h2bg" valign="top">
    <td width="100%">
      <br>
      <h2 align="center">GFF (General Feature Format) Specifications Document</h2>
    </td>
  </tr>
</table>

<a name="TOC">
<!-- INDEX BEGIN -->
<ul>
	<li><a href="#introduction">Introduction</a>
	<li><a href="#fields">Definition</a>
        <ul>
	   <li><a href="#standard_feature_table">Standard Table of Features</a>
           <li><a href="#attribute_field">Attribute Field</a>
	   <li><a href="#comments">Comments</a>
           <ul>
              <li><a href="#meta_info">Comments for Meta-Information</a>
           </ul>
	   <li><a href="#file_names">File Naming</a>
        </ul>
	<li><a href="#semantics">Semantics</a>
	<li><a href="#GFF_use">Ways to use GFF</a>
	<ul>
	   <li><a href="#examples">Complex Examples</a>
	   <ul>
              <li><a href="#homology_feature">Similarities to Other Sequences</a>
           </ul>
	   <li><a href="#cum_score_array">Cumulative Score Arrays</a>
	</ul>
	<li><a href="#mailing_list"> Mailing list</a>
	<li><a href="#edit_history">Edit History</a>
	<li><a href="#authors">Authors</a>
	<li><a href="index.shtml">Back to the GFF Home Page</a>
</ul>
<!-- INDEX END -->
<p>
<font color="#a00000">
<b>2000-9-29</b> The default version for GFF files is now Version 2.  This document has
been changed to show version 2 as default, with version one
alternatives shown where appropriate.
The main change from Version 1 to Version 2 is the requirement for a tag-value
type structure (essentially semicolon-separated .ace format) for any additional material on the
line, following the mandatory fields.  Version 2 also allows
'.' as a score, for features for which there is no score.  Dumping in version
2 format is implemented in ACEDB.
</font>
<hr>
<a name="introduction"><h3>Introduction</h3></a>
<p>
Essentially all current approaches to feature finding in higher organisms
use a variety of recognition methods that give scores to likely
signals (starts, splice sites, stops, motifs, etc.) or to extended regions
(exons, introns, protein domains etc.), and then combine these to give complete gene,
RNA transcript or protein structures.  Normally the combination step is done in the 
same program as the feature detection, often using dynamic programming methods.  To enable 
these processes to be decoupled, a format called GFF ('Gene-Finding
Format' or 'General Feature Format') 
was proposed as a protocol for the transfer of feature information.  
It is now possible to take features from an outside source and add them in to 
an existing program, or in the extreme to write a dynamic programming system 
which only took external features.
</p><p>
GFF allows people to develop features
and have them tested without having to maintain a complete
feature-finding system.  Equally, it would help those developing and
applying integrated gene-finding programs to test new feature
detectors developed by others, or even by themselves.
</p><p>
We want the GFF format to be easy to parse and process by a variety of
programs in different languages.  e.g. it would be useful if Unix
tools like grep, sort and simple perl and awk scripts could easily
extract information out of the file.  For these reasons, for the
primary format, we propose a record-based structure, where each
feature is described on a single line, and line order is not relevant.
</p><p>
We do not intend GFF format to be used for complete data management of
the analysis and annotation of genomic sequence.  Systems such as
Acedb, Genotator etc. that have much richer data representation
semantics have been designed for that purpose.  The disadvantages in
using their formats for data exchange (or other richer formats such as
ASN.1) are (1) they require more complexity in parsing/processing, (2)
there is little hope on achieving consensus on how to capture all
information.  GFF is intentionally aiming for a low common
denominator.
</p><p>
With the changes taking place to version 2 of the format, we also 
allow for feature sets to be defined over RNA and Protein sequences,
as well as genomic DNA.  This is used for example by the <a href="https://web.archive.org/web/20010208224442/http://www.uk.embnet.org/Software/EMBOSS/">EMBOSS</a> project to
provide standard format output for all features as an option.
In this case the &lt;strand&gt; and &lt;frame&gt; 
fields should be set to '.'.  To assist this transition in specification,
a new <a href="#Type_Meta_Comment">#Type Meta-Comment</a> has been added.
</p><p>
Here are some example records:
</p>
<pre>
SEQ1	EMBL	atg	103	105	.	+	0
SEQ1	EMBL	exon	103	172	.	+	0
SEQ1	EMBL	splice5	172	173	.	+	.
SEQ1	netgene	splice5	172	173	0.94	+	.
SEQ1	genie	sp5-20	163	182	2.3	+	.
SEQ1	genie	sp5-10	168	177	2.1	+	.
SEQ2	grail	ATG	17	19	2.1	-	0
</pre>
<p>
Back to <a href="#TOC">Table of Contents</a>
<p>
<hr>
<a name="fields"><h3>Definition</h3></a>

Fields are:
&#060;seqname&#062; &#060;source&#062; &#060;feature&#062; &#060;start&#062; &#060;end&#062; &#060;score&#062; &#060;strand&#062; &#060;frame&#062; [attributes] [comments]<p>
 <dl>

 <dt>&#060;seqname&#062; 
 <dd>The name of the sequence.  Having an explicit sequence name
allows a feature file to be prepared for a data set of multiple
sequences.  Normally the seqname will be the identifier of the
sequence in an accompanying fasta format file.  An alternative is that
&lt;seqname&gt; is the identifier for a sequence in a public database, such
as an EMBL/Genbank/DDBJ accession number.  Which is the case, and
which file or database to use, should be explained in accompanying
information.<p>

 <dt>&#060;source&#062; 
 <dd> The source of this feature.  This field will normally be used to
indicate the program making the prediction, or if it comes from public
database annotation, or is experimentally verified, etc.<p>

 <dt>&#060;feature&#062; 
 <dd> The feature type name.  We hope to suggest a standard set of
features, to facilitate import/export, comparison etc..  Of course,
people are free to define new ones as needed.  For example, Genie
splice detectors account for a region of DNA, and multiple detectors
may be available for the same site, as shown above.<p>
<a name="standard_feature_table">
<p> 
We would like to enforce a standard nomenclature for
common GFF features. This does not forbid the use of other features,
rather, just that if the feature is obviously described in the standard
list, that the standard label should be used. For this standard table
we propose to fall back on the international public standards for genomic 
database feature annotation, specifically, the 
<a href="https://web.archive.org/web/20010208224442/http://www3.ebi.ac.uk/Services/WebFeat/">
DDBJ/EMBL/GenBank feature table documentation</a>).</font></b>
</p><p>
 <dt>&#060;start&#062;, &#060;end&#062;
 <dd> Integers.  &#060;start&#062; must be less than or equal to
&#060;end&#062;.  Sequence numbering starts at 1, so these numbers
should be between 1 and the length of the relevant sequence,
inclusive. (<b>Version 2 change</b>: version 2 condones values of
&#060;start&#062; and &#060;end&#062; that extend outside the
reference sequence.  This is often more natural when dumping from
acedb, rather than clipping.  It means that some software using the
files may need to clip for itself.)<p>

 <dt>&#060;score&#062; 
 <dd> A floating point value.  When there is no score (i.e. for a
sensor that just records the possible presence of a signal, as for the
EMBL features above) you should use '.'. (<b>Version 2 change</b>: in
version 1 of GFF you had to write 0 in such circumstances.)<p>

 <dt>&#060;strand&#062; 
 <dd> One of '+', '-' or '.'.  '.' should be used when
strand is not relevant, e.g. for dinucleotide repeats.
<b>Version 2 change</b>: This field is left empty '.' for RNA and protein features.<p>

 <dt>&#060;frame&#062;
 <dd> One of '0', '1', '2' or '.'.  '0' indicates that the specified
region is in frame, i.e. that its first base corresponds to the first
base of a codon.  '1' indicates that there is one extra base,
i.e. that the second base of the region corresponds to the first base
of a codon, and '2' means that the third base of the region is the
first base of a codon.  If the strand is '-', then the first base of
the region is value of &#060;end&#062;, because the corresponding
coding region will run from &#060;end&#062; to &#060;start&#062; on
the reverse strand.  As with &#060;strand&#062;, if the frame is not
relevant then set &#060;frame&#062; to '.'.  
It has been pointed out that "phase" might be a better descriptor than
"frame" for this field.
<b>Version 2 change</b>: This field is left empty '.' for RNA and protein features.<p>

 <dt><a name="attribute_field">[attribute] </a>
<dd> From version 2 onwards, the attribute field 
must have an tag value structure following the syntax used within
objects in a .ace file, flattened onto one line by semicolon
separators.  Tags must be standard identifiers
([A-Za-z][A-Za-z0-9_]*).  Free text values must be quoted with double
quotes. <em>Note: all non-printing characters in such free text value strings
(e.g. newlines, tabs, control characters, etc)
must be explicitly represented by their C (UNIX) style backslash-escaped
representation (e.g. newlines as '\n', tabs as '\t').</em>
As in ACEDB, multiple values can follow a specific tag.  The
aim is to establish consistent use of particular tags, corresponding
to an underlying implied ACEDB model if you want to think that way
(but acedb is not required).  Examples of these would be:
<font size="3"><pre>
seq1     BLASTX  similarity   101  235 87.1 + 0	Target "HBA_HUMAN" 11 55 ; E_value 0.0003
dJ102G20 GD_mRNA coding_exon 7105 7201   .  - 2 Sequence "dJ102G20.C1.1"
</pre></font>
The semantics of tags in attribute field tag-values pairs has
intentionally not been formalized.  Two useful guidelines are to use 
DDBJ/EMBL/GenBank feature 'qualifiers' (see 
<a href="https://web.archive.org/web/20010208224442/http://www3.ebi.ac.uk/Services/WebFeat/">DDBJ/EMBL/GenBank
feature table documentation</a>), or the features that ACEDB generates 
when it dumps GFF.
</p><p>
<b>Version 1 note</b> In version 1 the attribute field was called the
group field, with the following specification:<br>
An optional string-valued field that can be used as a name to
group together a set of records.  Typical uses might be to group the
introns and exons in one gene prediction (or experimentally verified
gene structure), or to group multiple regions of match to another
sequence, such as an EST or a protein.
</p>

</dl>
All of the above described fields should be separated by TAB characters ('\t').
All values of the mandatory fields should not include whitespace
(i.e. the strings for &#060;seqname&#062;,
&#060;source&#062; and &#060;feature&#062; fields).
<p>
<b>Version 1 note</b> In version 1 each string had to be under 256
characters long, and the whole line should under 32k long.  This was
to make things easier for guaranteed conforming parsers, but seemed
unnecessary given modern languages.
<p>
Back to <a href="#TOC">Table of Contents</a>
<p>
<hr>
<a name="comments"><h3> Comments </h3>

Comments are allowed, starting with "#" as in Perl, awk etc.
Everything following # until the end of the line is ignored.
Effectively this can be used in two ways.  Either it must be at the
beginning of the line (after any whitespace), to make the whole line a
comment, or the comment could come after all the required fields on
the line.

<a name="meta_info"><h4> ## comment lines for meta information </h4>

There is a set of standardised (i.e. parsable) ## line types that can
be used optionally at the top of a gff file.  The philosophy is a
little like the special set of %% lines at the top of postscript
files, used for example to give the BoundingBox for EPS files.<p>

Current proposed ## lines are:

<dl>

  <dt><pre> ##gff-version 2 </pre>
  <dd> GFF version - in case it is a real success and we want to
change it.  The current default version is 2, so if this line is not
present version 2 is assumed.

  <dt><pre> ##source-version &lt;source&gt; &lt;version text&gt; </pre>
  <dd> So that people can record what version of a program or package was
used to make the data in this file. I suggest the version is text
without whitespace.  That allows things like 1.3, 4a etc.  There
  should be at most one source-version line per source.

  <dt> <pre> ##date &lt;date&gt; </pre>
  <dd> The date the file was made, or perhaps that the prediction
programs were run.  We suggest to use astronomical format: 1997-11-08
for 8th November 1997, first because these sort properly, and second
to avoid any US/European bias.

  <dt><a name="Type_Meta_Comment">
   <pre> ##Type &lt;type&gt; [&lt;seqname&gt;] </pre>
  <dd> The type of host sequence described by the features. Standard types
are 'DNA', 'Protein' and 'RNA'. The optional &lt;seqname&gt; allows multiple
##Type definitions describing multiple GFF sets in one file, each of
which have a distinct type. If the name is not provided,
then all the features in the file are of the given type. Thus, with this
meta-comment, a single file could contain DNA, RNA and Protein features,
for example, representing a single genomic locus or 'gene', alongside type-specific
features of its transcribed mRNA and translated protein sequences.
If no ##Type meta-comment is provided for a given GFF file, then the type
is assumed to be DNA.

<dt> <pre> 
 ##DNA &lt;seqname&gt;
 ##acggctcggattggcgctggatgatagatcagacgac
 ##...
 ##end-DNA
</pre>

<dd> To give a DNA sequence.  Several people have pointed out that it may
be convenient to include the sequence in the file.  It should not
become mandatory to do so, and in our experience this has been very
little used.  Often the seqname will be a well-known
identifier, and the sequence can easily be retrieved from a
database, or an accompanying file.  

<dt> <pre> 
 ##RNA &lt;seqname&gt;
 ##acggcucggauuggcgcuggaugauagaucagacgac
 ##...
 ##end-RNA
</pre>

<dd> Similar to DNA. Creates an implicit ##Type RNA
&lt;seqname&gt; directive.

<dt> <pre> 
 ##Protein &lt;seqname&gt;
 ##MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSF
 ##...
 ##end-Protein
</pre>

<dd> Similar to DNA.  Creates an implicit ##Type Protein
&lt;seqname&gt; directive.

<dt> <pre> ##sequence-region &lt;seqname&gt; &lt;start&gt; &lt;end&gt; </pre>
<dd> To indicate that this file only contains entries for the
specified subregion of a sequence.

</dl>

Please feel free to propose new ## lines.

The ## line proposal came out of some discussions including Anders
Krogh, David Haussler, people at the Newton Institute on 1997-10-29
and some email from Suzanna Lewis.  Of course, naive programs can
ignore all of these...

<a name="file_names"><h3> File Naming </h3>

We propose that the format is called "GFF", with conventional file
name ending ".gff".
<p>
Back to <a href="#TOC">Table of Contents</a>
<p>
<hr>
<a name="semantics"><h3> Semantics </h3>

We have intentionally avoided overspecifying the semantics of the
format.  For example, we have not restricted the items expressible in
GFF to a specified set of feature types (splice sites, exons etc.)
with defined semantics.  Therefore, in order for the information in a
gff file to be useful to somebody else, the person producing the
features must describe the meaning of the features.  <p>

In the example given above the feature "splice5" indicates that there
is a candidate 5' splice site between positions 172 and 173.  The
"sp5-20" feature is a prediction based on a window of 20 bp for the
same splice site.  To use either of these, you must know the position
within the feature of the predicted splice site.  This only needs to
be given once, possibly in comments at the head of the file, or in a
separate document.  <p>

Another example is the scoring scheme; we ourselves would like the
score to be a log-odds likelihood score in bits to a defined null
model, but that is not required, because different methods take
different approaches.

Avoiding a prespecified feature set also leaves open the possibility
for GFF to be used for new feature types, such as CpG islands,
hypersensitive sites, promoter/enhancer elements, etc.
<p>
Back to <a href="#TOC">Table of Contents</a>
<p>
<hr>
<a name="GFF_use"><h3> Ways to use GFF </h3>

Here are a few suggestions on how the GFF format might be used.
 <ol>
 <li> Simple sharing of sensors. In this case, researcher A has a sensor,
such as a 3' splice site sensor, and researcher B wants to test that
sensor.  They agree on a set of sequences, researcher A runs the
sensor on these sequences and sends the resulting GFF file to
researher B, who then evaluates the result.<p>

 <li> Representing experimental results.  GFF feature records can also
be created for experimentally confirmed exons and other features.  In
these cases there will presumably be no score.  Such "confirmed" GFF
files will be useful for evaluating predictions, using the same
software as you would to compare predictions.<p>

 <li> Integrated gene parsing. Several GFF files from different
researchers can be combined to provide the features used by an
integrated genefinder.  As mentioned above, this has the advantage
that different combinations of sensors and dynamic programming methods
for assembling sensor scores into consistent gene parses can be easily
explored.<p>

 <li> Reporting final predictions. GFF format can also be used to
communicate finished gene predictions. One simply reports final
predicted exons and other predicted gene features, either with their
original scores. or with some sort of posterior scores, rather than,
or in addition to, reporting all candidate gene features with their
scores.  To show that a set of the components belong to a single
prediction, a "attribute" field can be added to all the accepted sites.
This is useful for comparing the outputs of several integrated
genefinders among themselves, and to "confirmed" GFF files.  A
particular advantage of having the same format for both raw sensor
feature score files and final gene parse files is that one can easily
explore the possibility of combining the final gene parses from
several different genefinders, using another round of dynamic
programming, into a single integrated predicted parse.<p>

 <li> Visualisation. GFF will also provide a simple standard format for
standardising input to visualisation programs, showing predicted and
experimentally determined features, gene structures etc.

</ol>

<p>
Back to <a href="#TOC">Table of Contents</a>
<p>
<hr>
<a name="examples"><h3> Complex Examples</h3>

<a name="homology_feature">
<h4> Similarities to Other Sequences </h4>

A major source of information about a sequence comes from similarities
to other sequences.  For example, BLAST hits to protein sequences help
identify potential coding regions.  We can represent these as a set of
"similarity features": <font size="3"><pre>
seq1 BLASTX similarity 101 235 87.1 + 0    Target "HBA_HUMAN" 11 54 ; E_value 0.0003
</pre></font>
The proposed tag-value structure for gapped alignments is
<pre>
	Align &lt;seq_start&gt; &lt;target_start&gt; [&lt;length&gt;] ;
</pre>
to define each ungapped block in the alignment, with multiple Align
tags to give a full gapped alignment.  The &lt;length&gt; field is
optional because in its absence a block is presumed to extend until it
reaches the next specified block, or the end of the complete
similarity.  This corresponds to the standard case with alignments
that they don't have simultaneous gaps on both strands.  For example,
for the above HBA_HUMAN similarity, the Align information could be
<pre>
	Align 101 11 ; Align 179 36 ;
</pre>
which leaves the DNA triplet from 176 to 178 aligned to a gap in the
protein sequence.
<p>
Back to <a href="#TOC">Table of Contents</a>
<p>
<hr>
<a name="cum_score_array"><h3> Cumulative Score Arrays </h3>

One issue that comes up with a record-based format such as the GFF
format is how to cope with large numbers of overlapping segments.  For
example, in a long sequence, if one tries to include a separate record
giving the score of every candidate exon, where a candidate exon is
defined as a segment of the sequence that begins and ends at candidate
splice sites and consists of an open reading frame in between, then
one can have an infeasibly large number of records.  The problem is
that there can be a huge number of highly overlapping exon
candidates. <p>

Let us assume that the score of an exon can be decomposed into three
parts: the score of the 5' splice site, the score of the 3' splice
site, and the sum of the scores of all the codons in between. In such
a case it can be much more efficient to use the GFF format to report
separate scores for the splice site sensors and for the individual
codons in all three (or six, including reverse strand) frames, and let
the program that interprets this file assemble the exon scores.  The
exon scores can be calculated efficiently by first creating three
arrays, each of which contains in its [i]th position a value A[i] that
is the partial sum of the codon scores in a particular frame for the
entire sequence from position 1 up to position i.  Then for any
positions i &lt; j, the sum of the scores of all codons from i to j can
be obtained as A[j] - A[i]. Using these arrays, along with the
candidate splice site scores, a very large number of scores for
overlapping exons are implicitly defined in a data structure that
takes only linear space with respect to the number of positions in the
sequence, and such that the score for each exon can be retrieved in
constant time. <p>

When the GFF format is used to transmit scores that can be summed for
efficient retrieval as in the case of the codon scores above, we ask
that the provider of the scores indicate that these scores are
summable in this manner, and provide a recipe for calculating the
scores that are to be derived from these summable scores, such as the
exon scores described above. We place no limit on the complexity of
this recipe, nor do we provide a standard protocol for such assembly,
other than providing examples.  It behooves the sensor score provider
to keep the recipe simple enough that others can easily implement it.
<p>
Back to <a href="#TOC">Table of Contents</a>
<p>
<hr>

<a name="mailing_list"><h3> Mailing list </h3>
<p>
There is a <a href="https://web.archive.org/web/20010208224442/mailto:gff-list@sanger.ac.uk"> mailing list </a>
to which you can send comments, enquiries, complaints etc. about GFF.
If you want to be added to the mailing list, please send
mail to <a href="https://web.archive.org/web/20010208224442/mailto:Majordomo@sanger.ac.uk">Majordomo@sanger.ac.uk</a> with the 
following command in the body of your email message:
<p>

<code>
    subscribe gff-list
</code>
<p>
<p>
Back to <a href="#TOC">Table of Contents</a>
<p>
<hr>
<a name="edit_history"><h3>Edit History</h3></a>
<p>
000929 rd: make version 2 default and propose Align tag-value syntax
<p>
0003022 rbsk: small clarification to #comment rules
<p>
991711 rbsk: (overdue changes as per September '99 gff-list commentaries)
<ul>
   <li>GFF acronym renamed to mean 'General Feature Format' rather than just
        'Gene-Finding Features', in order to conceptually accommodate RNA and 
        Protein as well as DNA features </li>
   <li>added ##Type metacomment field, </li>
   <li>changed name of [group] field to [attribute] field</li>
</ul>
<p>
990816 rbsk: standard list of features and group tags (first attempt at clarification)      
<p>
990317 rbsk:
<ul>
   <li>End of line comments following Version 2 [group] field tag-value structures must be 
       tab '\t' or hash '#' delimited.
</ul>
<p>
990226 rbsk: incorporated amendments to the version 2 specification as follows:
<ul>
     <li>Non-printing characters (e.g. newlines, tabs) in Version 2 double quoted
"free text values" must be explicitly represented by their C (UNIX) style 
backslash escaped character (i.e. '\t' for tabs, '\n' for newlines, etc.)<br>
     <li>Removed field (256) and line (32K) character size limitations for Version 2.
     <li>Removed arbitrary whitespace field delimiter permission from specification.
TAB ('\t') field delimiters now enforced again, as in Version 1.<br>
</ul>
<p>
981216 rd: introduced version 2 changes.
<p>
980909 ihh: fixed some small things and put this page on the Sanger
GFF site.
<p>
971113 rd: added section on mailing list.
<p>
971113 rd: added extra "source" field as discussed at Newton Institute
meeting 971029.  There are two main reasons.  First, to help prevent
name space clashes -- each program would have their own source
designation.  Second, to help reuse feature names, so one could have
"exon" for exon predictions from each prediction program.
<p>
971108 rd: added ## line proposals - moved them into main text 971113.
<p>
971028 rd: I added the section about name space.
<p>
971028 rd: I considered switching from start-end notation to
start-length notation, on the suggestion of Anders Krogh.  This seems
nicer in many cases, but is a debatable point.  I then switched back!
<p>
971028 rd: We also now allow extra text after &#060;group&#062;
without a comment character, because this immediately proved useful.
<p>
971028 rd: I changed the comment initiator to '#' from '//' because a 
single symbol is easier for simple parsers.<p>
<p>
Back to <a href="#TOC">Table of Contents</a>
<p>
<hr>
<a name="authors"><h3>Authors</h3></a>
<p>
GFF Protocol Specification initially proposed by: 
<a href="https://web.archive.org/web/20010208224442/mailto:rd@sanger.ac.uk">Richard Durbin</a> and 
<a href="https://web.archive.org/web/20010208224442/mailto:haussler@cse.ucsc.edu">David Haussler</a>
<p>with amendments proposed by: 
<a href="https://web.archive.org/web/20010208224442/mailto:lstein@cshl.org">Lincoln Stein</a>, Suzanna Lewis, Anders Krogh and others.
<p>
Back to <a href="#TOC">Table of Contents</a>
<p>
        <br></td>
      </tr>
    </table>
<!-- page content ends here -->
    <table border="0" cellpadding="0" cellspacing="0" width="100%" align="center">
      <tr valign="top">
        <td colspan="2" bgcolor="#aaaaaa"><img src="/web/20010208224442im_/http://www.sanger.ac.uk/icons/blank.gif" height="1" width="640" alt=""></td>
      </tr>
      <tr valign="top">
        <td colspan="2" bgcolor="#f5f5ff"><img src="/web/20010208224442im_/http://www.sanger.ac.uk/icons/blank.gif" height="2" width="640" alt=""></td>
      </tr>
      <tr valign="top" bgcolor="#f5f5ff">
        <td align="left" class="headerinactive">&nbsp;last modified 04-Dec-2000, 01:01 PM</td>
	<td align="right" class="headerinactive"><a href="https://web.archive.org/web/20010208224442/http://www.sanger.ac.uk/feedback/">webmaster@sanger.ac.uk</a>&nbsp;</td>
      </tr>
      <tr valign="top">
        <td colspan="2" bgcolor="#f5f5ff"><img src="/web/20010208224442im_/http://www.sanger.ac.uk/icons/blank.gif" height="2" width="640" alt=""></td>
      </tr>
      <tr valign="top">
        <td colspan="2" bgcolor="#aaaaaa"><img src="/web/20010208224442im_/http://www.sanger.ac.uk/icons/blank.gif" height="1" width="640" alt=""></td>
      </tr>
    </table>
  </body>
</html>


<!--
     FILE ARCHIVED ON 22:44:42 Feb 08, 2001 AND RETRIEVED FROM THE
     INTERNET ARCHIVE ON 12:48:32 Nov 12, 2018.
     JAVASCRIPT APPENDED BY WAYBACK MACHINE, COPYRIGHT INTERNET ARCHIVE.

     ALL OTHER CONTENT MAY ALSO BE PROTECTED BY COPYRIGHT (17 U.S.C.
     SECTION 108(a)(3)).
-->
<!--
playback timings (ms):
  LoadShardBlock: 69.805 (3)
  esindex: 0.009
  captures_list: 91.845
  CDXLines.iter: 12.965 (3)
  PetaboxLoader3.datanode: 81.123 (4)
  exclusion.robots: 0.531
  exclusion.robots.policy: 0.509
  RedisCDXSource: 2.189
  PetaboxLoader3.resolve: 153.316
  load_resource: 201.369
-->
