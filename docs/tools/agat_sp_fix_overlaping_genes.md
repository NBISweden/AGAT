# agat\_sp\_fix\_overlaping\_genes.pl

## DESCRIPTION

Check a gtf/gff annotation file to find cases where differents gene features
have CDS that overlap. In this case the gene features will be merged in only one.
One gene is choosen as reference, and the mRNA from the other gene will be linked to it.
So, it creates isoforms.

## SYNOPSIS

```
agat_sp_fix_overlaping_genes.pl -f infile  [-o outfile]
agat_sp_fix_overlaping_genes.pl --help
```

## OPTIONS

- **-f**, **--file**, **--gff3** or **--gff**

    Input GTF/GFF file.

- **-o**, **--out**, **--output** or **--outfile**

    Output file. If none given, will be display in standard output.

- **--help** or **-h**

    Display this helpful text.

