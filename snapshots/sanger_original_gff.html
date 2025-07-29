<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 3.2//EN">
<html>
Copy from <a href="https://web.archive.org/web/19990209163146/http://www.sanger.ac.uk:80/Software/GFF/">web.archive.org</a>.  
  
Other interesting information from <a href="https://web.archive.org/web/20010619132027/http://www.sanger.ac.uk:80/Software/formats/GFF/">web.archive.org</a>.

  <h1>
  The Sanger Centre : Gene-Finding Format
  </h1>

<!-- page content starts here -->

<h2 align="CENTER">GFF: an exchange format for gene-finding features</h2>


<a href="gff.shtml">GFF</a> (Gene-Finding Features) is a format specification
for describing genes and other features associated with genomic sequences.
This page is a starting-point for finding out about this format and its
use in bioinformatics.
In particular, since its proposal a considerable amount of software has
been developed for use with GFF and this page is intended as a focus
for the collation of this software, whether developed in the Sanger Centre
or elsewhere.
<p>

A GFF record is an extension of a basic <i>(name,start,end)</i> tuple (or "NSE")
 that can be used to identify a substring of a biological sequence.
(For example, the NSE <i>(ChromosomeI,2000,3000)</i> specifies the
 third kilobase of the sequence named "ChromosomeI".)
GFF allows for moderately verbose annotation of single NSEs.
It also provides limited support for NSE pairs in a rather asymmetrical way.
An alternative format for representing NSE pairs that is used by several of the programs listed below
 is EXBLX, as used by <a href="https://web.archive.org/web/19990209163146/http://bioweb.pasteur.fr/docs/softgen.html#MSPCRUNCH">MSPcrunch</a>
 (Sonnhammer and Durbin (1994), "An expert system for processing sequence homology data", Proceedings of ISMB 94, 363-368).
<p>

The most common operations that one tends to want to perform
 on sets of NSEs and NSE-pairs include intersection, exclusion, union, filtration, sorting, transformation
 (to a new co-ordinate system) and dereferencing (access to the described sequence).
With a suitably flexible definition of NSE "similarity", these operations form a basis for more
 sophisticated algorithms like clustering and joining-together by dynamic programming.
Programs to perform all of these tasks are described below, with links to local copies.
<p>

Criticism of and new links for this page are always welcome. Please contact the page
 administrator, whose email address appears at the foot of the page.
<p>

<ul>

<li> <a href="gff.shtml">Introduction to and full specification of GFF</a> - read this first, if you're new to GFF.
<p>

<li><a href="https://web.archive.org/web/19990209163146/http://www.cse.ucsc.edu/~haussler/genefindingpaper/paper.html">Review of gene-finding methods</a>
 by David Haussler.
The <a href="https://web.archive.org/web/19990209163146/http://www.cse.ucsc.edu/research/compbio/research.html">UCSC computational biology research
 projects page</a> has many interesting links.
<p>

<li><a href="mailinglist.shtml">The GFF mailing list</a>
- a mailing list for discussion related to the GFF feature file format.
To place a message on this list send it to
<a href="https://web.archive.org/web/19990209163146/mailto:gff-list@sanger.ac.uk">gff-list@sanger.ac.uk</a>.
To join, send email to <a href="https://web.archive.org/web/19990209163146/mailto:majordomo@sanger.ac.uk">majordomo@sanger.ac.uk</a>
with the following command in the body of your email message: 
<p>

<code>
subscribe gff-list
</code>
<p>

<li> GFF Perl Modules
- broad-functionality Perl5.0 modules developed by <a href="https://web.archive.org/web/19990209163146/mailto:th@sanger.ac.uk">Tim Hubbard</a>
at the Sanger Centre, including: <ul>
<li> <a href="https://web.archive.org/web/19990209163146/http://www.sanger.ac.uk/Software/PerlModule/GFF.html">GFF</a>
<li> <a href="https://web.archive.org/web/19990209163146/http://www.sanger.ac.uk/Software/PerlModule/GFFPair.html">GFFPair</a>
<li> <a href="https://web.archive.org/web/19990209163146/http://www.sanger.ac.uk/Software/PerlModule/GeneFeature.html">GeneFeature</a>
<li> <a href="https://web.archive.org/web/19990209163146/http://www.sanger.ac.uk/Software/PerlModule/GeneFeaturePair.html">GeneFeaturePair</a>
</ul>
<p>

<li> <a href="gffdp.pl">GFF dynamic programming: <code>gffdp.pl</code></a> - a Perl program for joining together GFF
segments using <b>Generalised Hidden Markov Models</b> with stacks,
written by <a href="https://web.archive.org/web/19990209163146/mailto:ihh@sanger.ac.uk">Ian Holmes</a>.
(Requires the <a href="BraceParser.pm"><code>BraceParser.pm</code></a> module.)
The architecture and scoring schemes of the underlying models are
entirely flexible and can be specified in a separate file.
Example model files include: <ul>
<li> <a href="gene.model"><code>gene.model</code></a> - a model for assembling exon predictions
<li> <a href="transposon.model"><code>transposon.model</code></a> - a model for finding DNA transposons (or indeed any proteins flanked by inverted repeats)
</ul>
More information about this program is available on request.
<p>

<li> <a href="bigdp.tar.gz">EXBLX dynamic programming: <code>bigdp</code></a> - a C++ program that assembles EXBLX segments using an affine gap penalty by doing linear-space divide-and-conquer dynamic programming, written by Ian Holmes. The program does not examine the sequences to which the EXBLX data refer, but finds optimal connections between the segments given their co-ordinates. GFF pair format can be converted to EXBLX using <a href="gff2exblx.pl"><code>gff2exblx.pl</code></a>.
<p>
EXBLX records are single lines comprising eight whitespace-delimited fields:
 (SCORE, PERCENT-ID, START#1, END#1, NAME#1, START#2, END#2, NAME#2).
<code>bigdp</code> requires that the two NSEs are the same length (i.e. END#1- START#1= END#2- START#2).
The output of <code>bigdp</code> is modified EXBLX.
Each line of the ouput describes a set of several input segments joined together;
 the percent-ID field is replaced by the number of input segments that were used
 and a ninth field, compactly describing the co-ordinates of the input segments, is added.
The algorithm used by the program is documented more fully in my <a href="https://web.archive.org/web/19990209163146/http://www.sanger.ac.uk/Users/ihh/thesis.ps.gz">PhD thesis</a>.
<p>

<li><a href="gffhitcount.tar.gz"><code>gffhitcount</code></a> - a C++ program that counts the number of times each base in a set of sequences is spanned by a GFF record and returns the results in GFF format.
<p>

<li>Miscellaneous Perl scripts:<ul>
 <li><a href="gffintersect.pl"><code>gffintersect.pl</code></a> - efficiently finds the intersection (or exclusion) of two GFF streams, reporting intersection information in the Group field. Definition of "intersection" allows for near-neighbours and minimum-overlap
 <li><a href="intersectlookup.pl"><code>intersectlookup.pl</code></a> - used with <a href="gffintersect.pl"><code>gffintersect.pl</code></a> to do reverse lookups and other manipulations on the results of an intersection test. Useful for e.g. pruning the lowest-scoring redundant entries from a GFF file
 <li><a href="gffmask.pl"><code>gffmask.pl</code></a> - uses a GFF file to mask out specified sections of a FASTA-format DNA database with "n"'s (or any other character)
 <li><a href="gfftransform.pl"><code>gfftransform.pl</code></a> - transforms a GFF stream from one co-ordinate system to another (e.g. from clone to chromosome co-ordinates), given another GFF file describing the transformation. Requires <a href="GFFTransform.pm"><code>GFFTransform.pm</code></a>
 <li><a href="gff2seq.pl"><code>gff2seq.pl</code></a> - given chromosome co-ordinates, a clone database and a physical map co-ordinate file, returns the specified section of chromosomal sequence, even if it spans multiple clones. Requires <a href="SeqFileIndex.pm"><code>SeqFileIndex.pm</code></a> and <a href="FileIndex.pm"><code>FileIndex.pm</code></a>
 <li><a href="gfffilter.pl"><code>gfffilter.pl</code></a> - filters lines out of a GFF stream according to user-specified criteria
 <li><a href="gffsort.pl"><code>gffsort.pl</code></a> - sorts GFF streams by sequence name and startpoint
 <li><a href="gffmerge.pl"><code>gffmerge.pl</code></a> - merges sorted GFF streams
 <li><a href="cluster2gff.pl"><code>cluster2gff.pl</code></a> - converts a list of whitespace-separated NSE clusters (in the format "name/start-end") into a GFF data set.
 <li><a href="exblxgffintersect.pl"><code>exblxgffintersect.pl</code></a> - similar to <a href="gffintersect.pl"><code>gffintersect.pl</code></a>, but finds NSE pairs in an EXBLX file that intersect with single NSEs in a GFF file. Useful for e.g. filtering out all hits between known genes from an all-vs-all BLAST comparison of genomic DNA
 <li><a href="GFFTransform.pm"><code>GFFTransform.pm</code></a> - module to convert between GFF co-ordinate systems. Used by <a href="gfftransform.pl"><code>gfftransform.pl</code></a>, <a href="blasttransform.pl"><code>blasttransform.pl</code></a> and <a href="exblxtransform.pl"><code>exblxtransform.pl</code></a>
 <li><a href="SeqFileIndex.pm"><code>SeqFileIndex.pm</code></a> - module to access a clone database using a map file. Requires <a href="FileIndex.pm"><code>FileIndex.pm</code></a>. Used by <a href="gff2seq.pl"><code>gff2seq.pl</code></a>
 <li><a href="FileIndex.pm"><code>FileIndex.pm</code></a> - module to build a quick lookup table for flatfiles. Used by <a href="exblxsym.pl"><code>exblxsym.pl</code></a>, <a href="gff2seq.pl"><code>gff2seq.pl</code></a> and <a href="SeqFileIndex.pm"><code>SeqFileIndex.pm</code></a>
 <li><a href="BraceParser.pm"><code>BraceParser.pm</code></a> - module to parse <a href="gffdp.pl"><code>gffdp.pl</code></a> model files, wherein fields are enclosed by braces {like this}
</ul>
<p>
Several of these scripts duplicate functionality provided
 by Tim Hubbard's perl modules (see above), but may be less algorithmically complex
 (a significant consideration for chromosome-sized GFF files!).
<p>
Please do email Ian Holmes if you require documentation for these programs.
<p>

<li>Programs that are only tangentially related to GFF, but complement the GFF tools well:<ul>
 <li><a href="exblxsym.pl"><code>exblxsym.pl</code></a> - symmetrises an EXBLX file (ensures that for every A:B pair there is a single corresponding pair B:A)
 <li><a href="exblxasym.pl"><code>exblxasym.pl</code></a> - asymmetrises an EXBLX file (filters through only those pairs A:B for which B&gt;A)
 <li><a href="exblxcluster.pl"><code>exblxcluster.pl</code></a> - builds optimal clusters from an EXBLX stream
 <li><a href="exblxfastcluster.pl"><code>exblxfastcluster.pl</code></a> - builds clusters from an EXBLX stream using a fast incremental heuristic
 <li><a href="seqcluster.pl"><code>seqcluster.pl</code></a> - builds optimal clusters from an EXBLX stream, ignoring sequence start and endpoint
 <li><a href="exblxindex.pl"><code>exblxindex.pl</code></a> - builds a quick lookup index for an EXBLX file
 <li><a href="exblxsingles.pl"><code>exblxsingles.pl</code></a> - filters through only non-overlapping entries from an EXBLX stream
 <li><a href="exblxsort.pl"><code>exblxsort.pl</code></a> - sorts an EXBLX stream
 <li><a href="exblxtidy.pl"><code>exblxtidy.pl</code></a> - tidies up an EXBLX stream (joins overlapping matches, prunes out lines corresponding to BLAST errors, etc.)
 <li><a href="exblxtransform.pl"><code>exblxtransform.pl</code></a> - transforms from one co-ordinate system to another (e.g. clones to chromosomes). Requires <a href="GFFTransform.pm"><code>GFFTransform.pm</code></a>
 <li><a href="cfilter.pl"><code>cfilter.pl</code></a> - flags low-complexity regions in a FASTA DNA database. The complexity is calculated as the entropy of variable-length oligomer composition in a variable-length sliding window
 <li><a href="blasttransform.pl"><code>blasttransform.pl</code></a> - BLASTs a clone database against itself then transforms, sorts and merges the results into chromosome co-ordinates according to a physical (sequence) map file, which is in GFF format. Requires <a href="GFFTransform.pm"><code>GFFTransform.pm</code></a>
 <li><a href="SequenceIterator.pm"><code>SequenceIterator.pm</code></a> - module to assist iterations on FASTA DNA databases; creates temporary files for each sequence
</ul>
<p>

<li>Output format conversion utilities:<ul>
 <li><a href="hmm2gff.pl">HMMER 1.7 to GFF</a>
 <li><a href="hmmsearch2gff.pl">HMMER 2.0 to GFF</a>
 <li><a href="exblx2gff.pl">EXBLX to GFF</a>
 <li><a href="gff2exblx.pl">GFF to EXBLX</a>
 <li><a href="spcwise2gff.pl">GeneWise to GFF</a>
 <li><a href="scan2gff.pl">GCG's <code>scan</code> to GFF</a>
 <li><a href="tandem2gff.pl">GCG's <code>tandem</code> to GFF</a>
</ul>
<p>

</ul>



<!-- page content ends here -->

</td></tr></table></center>  <!-- close table for page content -->

 <hr align="CENTER" width="90%">

<!-- open table for page footer -->
<table border="0" width="100%">
 <tr>
  <td align="LEFT">
   <i>
   last modified : 04-Jan-1999, 07:46 PM
   </i>
  </td>

  <td align="RIGHT">
   <a href="/web/19990209163146/http://www.sanger.ac.uk:80/Users/ihh/">Ian Holmes</a>
   <i>(<a href="https://web.archive.org/web/19990209163146/mailto:ihh@sanger.ac.uk">ihh@sanger.ac.uk</a>)</i>
  </td>
 </tr>
</table>  <!-- close table for page footer -->

</body>
</html>


<!--
     FILE ARCHIVED ON 16:31:46 Feb 09, 1999 AND RETRIEVED FROM THE
     INTERNET ARCHIVE ON 11:24:11 Oct 29, 2018.
     JAVASCRIPT APPENDED BY WAYBACK MACHINE, COPYRIGHT INTERNET ARCHIVE.

     ALL OTHER CONTENT MAY ALSO BE PROTECTED BY COPYRIGHT (17 U.S.C.
     SECTION 108(a)(3)).
-->
<!--
playback timings (ms):
  LoadShardBlock: 51.404 (3)
  esindex: 0.014
  captures_list: 354.797
  CDXLines.iter: 10.86 (3)
  PetaboxLoader3.datanode: 50.977 (4)
  exclusion.robots.fetch: 284.623 (4)
  exclusion.robots: 285.088
  exclusion.robots.policy: 0.184
  RedisCDXSource: 1.678
  PetaboxLoader3.resolve: 23.166
  load_resource: 31.258
-->
