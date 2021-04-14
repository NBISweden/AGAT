# NAME

agat\_sp\_manage\_functional\_annotation.pl

# DESCRIPTION

The script take a gff3 file as input and blast and/or interpro output in order
to attach functional annotation to corresponding features within the gff file.

\>The blast against Protein Database (outfmt 6) allows to fill the field/attribute
NAME for gene and PRODUCT for mRNA.

\>The Interpro result (.tsv) file allows to fill the DBXREF field/attribute with
pfam, tigr, interpro, GO, KEGG, etc... terms data.

With the &lt;id> option the script will change all the ID field by an Uniq ID
created from the given prefix, a letter to specify the kind of feature (G,T,C,E,U),
and the feature number.

The result is written to the specified output file, or to STDOUT.

About the TSV format from interproscan:
&#x3d;======================================

The TSV format presents the match data in columns as follows:

```
1.Protein Accession (e.g. P51587)
2.Sequence MD5 digest (e.g. 14086411a2cdf1c4cba63020e1622579)
3.Sequence Length (e.g. 3418)
4.Analysis (e.g. Pfam / PRINTS / Gene3D)
5.Signature Accession (e.g. PF09103 / G3DSA:2.40.50.140)
6.Signature Description (e.g. BRCA2 repeat profile)
7.Start location
8.Stop location
9.Score - is the e-value (or score) of the match reported by member database method (e.g. 3.1E-52)
10.Status - is the status of the match (T: true)
11.Date - is the date of the run
12.(InterPro annotations - accession (e.g. IPR002093) - optional column; only displayed if -iprlookup option is switched on)
13.(InterPro annotations - description (e.g. BRCA2 repeat) - optional column; only displayed if -iprlookup option is switched on)
14.(GO annotations (e.g. GO:0005515) - optional column; only displayed if --goterms option is switched on)
15.(Pathways annotations (e.g. REACT_71) - optional column; only displayed if --pathways option is switched on)
```

P.S: The 9th column contains most of time e-value, but can contain also score (e.g Prosite). To understand the difference: https://myhits.isb-sib.ch/cgi-bin/help?doc=scores.html

About the outfmt 6 from blast:
&#x3d;=============================

```perl
1.  qseqid  query (e.g., gene) sequence id
2.  sseqid  subject (e.g., reference genome) sequence id
3.  pident  percentage of identical matches
4.  length  alignment length
5.  mismatch  number of mismatches
6.  gapopen   number of gap openings
7.  qstart  start of alignment in query
8.  qend  end of alignment in query
9.  sstart  start of alignment in subject
10.   send  end of alignment in subject
11.   evalue  expect value
12.   bitscore  bit score
```

Currently the best e-value win... That means another hit with a lower e-value
(but still over the defined threshold anyway) even if it has a better PE value
will not be reported.

# SYNOPSIS

```
agat_sp_manage_functional_annotation.pl -f infile.gff [ -b blast_infile --db uniprot.fasta -i interpro_infile.tsv --id ABCDEF --output outfile ]
agat_sp_manage_functional_annotation.pl --help
```

# OPTIONS

- **-f**, **--reffile**,**-ref** , **--gff** or **--gff3**

    String - Input GTF/GFF file.

- **-b** or **--blast**

    String - Input blast ( outfmt 6 = tabular ) file that will be used to complement the features
    read from the first file (specified with --ref).

- **-d** or **--db**

    String - The fasta file that has been used as DB for the blast. Gene names and products/descriptions will be fished from this file.

- **--be** or **--blast\_evalue**

    Integer - Maximum e-value to keep the annotation from the blast file. By default 1e-6.

- **--pe**

    Integer - The PE (protein existence) in the uniprot header indicates the type of evidence that supports the existence of the protein.
    You can decide until which protein existence level you want to consider to lift the finctional information. Default 5.

    1\. Experimental evidence at protein level
    2\. Experimental evidence at transcript level
    3\. Protein inferred from homology
    4\. Protein predicted
    5\. Protein uncertain

- **-i** or **--interpro**

    String - Input interpro file (.tsv) that will be used to complement the features read from
    the first file (specified with **--ref**).

- **-id**

    String - This option will changed the id name. It will create from id prefix (usually 6 letters) given as input, uniq IDs like prefixE00000000001. Where E mean exon. Instead E we can have C for CDS, G for gene, T for mRNA, U for Utr.
    In the case of discontinuous features (i.e. a single feature that exists over multiple genomic locations) the same ID may appear on multiple lines. All lines that share an ID collectively represent a signle feature.

- **-idau**

    Boolean - This option (id all uniq) is similar to -id option but Id of features that share an ID collectively will be change by different and uniq ID.

- **-nb**

    Integer - Usefull only if -id is used.
    This option is used to define the number that will be used to begin the numbering. By default begin by 1.

- **-o** or **--output**

    String - Output GFF file.  If no output file is specified, the output will be
    written to STDOUT.

- **-v**

    Boolean - Verbose, for debug purpose.

- **-h** or **--help**

    Boolean - Display this helpful text.

