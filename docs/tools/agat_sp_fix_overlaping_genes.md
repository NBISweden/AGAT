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

- **-c** or **--config**

    String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any,
    otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
    The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

- **--help** or **-h**

    Display this helpful text.
