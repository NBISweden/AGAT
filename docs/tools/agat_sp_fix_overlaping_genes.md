# agat\_sp\_fix\_overlaping\_genes.pl

## DESCRIPTION

Check a GTF/GFF annotation file to find cases where different gene features
have CDS that overlap. In this case the gene features will be merged in only one.
One gene is chosen as reference, and the mRNA from the other gene will be linked to it.
So, it creates isoforms.

## SYNOPSIS

```
agat_sp_fix_overlaping_genes.pl -f infile  [-o outfile]
agat_sp_fix_overlaping_genes.pl --help
```

## OPTIONS

- **-f**, **--file**, **--gff3** or **--gff**

    Input GTF/GFF file.

- **-m** or **--merge**

    Bolean: Merge/add the attributes of gene feature that are merged (except ID and Parent).

- **-o**, **--out**, **--output** or **--outfile**

    Output file. If none given, will be display in standard output.

- **-v** or **--verbose**

    BOLEAN: Add verbosity.

- **--help** or **-h**

    Display this helpful text.
