# NAME

agat\_sp\_prokka\_fragmented\_gene\_annotations.pl

# DESCRIPTION

The script aims to look at fragmented gene annotations (FRAGS) within prokka annotations.
The FRAGS represent two (or more) ORFs that are in close proximity and are annotated
with homology to the same gene. In such cases, Prokka ads an \_n suffix to the gene ID.
For example, a splitted genX can then be found as genX\_1 and genX\_2 in the GFF.
See here for a case: https://github.com/tseemann/prokka/issues/502

\* The script will inform you how many case there is in your annotation.
\* If you think the FRAGS is due to a sequencing error (frameshift due to short indel),
using the --frags parameter will fix the FRAGS if genX\_1 and genX\_2 are not in the same frame.
The gff and the fasta file will be modified. The gene are merged, an insertion of
one or two N will be added in between the genes to fix the frameshift.
\* If you think the FRAGS is not due to a sequencing error, use the --pseudo parameter,
the gff will be fix (gene merged) and the agat\_pseudo attribute (the value is the position of the codon stop)
will be added to the related features.
\* using --frags and --pseudo is similar to use only --frags, except when no frameshift
is found for a detected FRAGS (both gene are in the same frame), the agat\_pseudo
attribute is also added to the related features.

How the tool detecte the FRAGS?
\* Search for cases where contiguous genes have the same name (e.g. lpxA\_1 lpxA\_2).
\* If so we look at the size of the protein of each of those genes (lpxA\_1 AA=175 ; lpxA\_2 AA=116),
and compute the size when merged togeter (devoided of the overlap if any) => here 270 AA
\* Then we look at the size of the protein used to infer the name (lpxA\_1 inferred from Q9PIM1 = 263 AA ; lpxA\_2 inferred from P0A722 = 262 AA )
and compute the average length of the reference protein: here 262AA. We add 20% to the length to be sure to include border cases => 282AA.
\* Compare the length of the merged proteins (262 AA) against the reference protein length (282).
If the the expected protein length (282 AA) is longer we have a FRAGS.

# SYNOPSIS

```
agat_sp_prokka_fragmented_gene_annotations.pl -gff infile.gff --fasta genome.fa --db prokka/prokka_bacteria_sprot.fa  -o outfolder
agat_sp_prokka_fragmented_gene_annotations.pl --help
```

# OPTIONS

- **--gff**

    Input genome GTF/GFF file. Mandatory.

- **-f**, **--fa** or **--fasta**

    Input genome fasta file. Mandatory.

- **--db**

    Input Uniprot fasta file used by prokka. Mandatory.

- **--frags**

    Merge and fix detected FRAGS if not in the same frame

- **--pseudo**

    Merge detected FRAGS and add the agat\_pseudo attribute (value will be the location of the first stop codon met).

- **--hamap\_size**

    Some protein function are not infered by Uniprot but by Hamap. In such case the information
    is retrieved from the web. As hamap provide a family profile, the protein size if a range.
    "low" option will use the low value of the range,
    "middle" option will use the average of the range,
    "high" option will the the high value of the range.
    Default "high".

- **--ct**, **--codon** or **--table**

    Codon table to use. \[default 1\]

- **--skip\_hamap**

    For test purpose it could be useful to skip hamap, because it requires fetching information through internet.

- **-o** , **--output** or **--out**

    Output folder. Mandatory.

- **-v**

    verbose mode. Default off.

- **-h** or **--help**

    Display this helpful text.

