## NAME

agat\_convert\_sp\_gff2tsv.pl

## DESCRIPTION

The script aims to convert gtf/gff file into tabulated file.
Attribute's tags from the 9th column become column titles.

## SYNOPSIS

```
agat_convert_sp_gff2tsv.pl -gff file.gff [ -o outfile ]
agat_convert_sp_gff2tsv.pl --help
```

## OPTIONS

- **--gff** or **-f**

    Input GTF/GFF file.

- **-c** or **--ct**

    When the features doesn't have Parent/ID relationships, the parser will try to group
    features using a common/shared attribute (i.e. a locus tag.). By default locus\_tag and gene\_id.
    You can provide another specific common/shared attribute using this option.

- **--ml** or **--merge\_loci**

    Merge loci parameter, default deactivated. You turn on the parameter if you want to merge loci into one locus when they overlap.
    (at CDS level for mRNA, at exon level for other level2 features. Strand has to be the same). Prokaryote can have overlaping loci so it should not use it for prokaryote annotation.
    In eukaryote, loci rarely overlap. Overlaps could be due to error in the file, mRNA can be merged under the same parent gene if you acticate the option.

- **-o** , **--output** , **--out** or **--outfile**

    Output GFF file.  If no output file is specified, the output will be
    written to STDOUT.

- **-h** or **--help**

    Display this helpful text.

