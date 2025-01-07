# agat_sp_filter_by_ORF_size.pl

## DESCRIPTION

The script reads a gff annotation file, and create two output files,
one contains the gene models with ORF passing the test, the other contains the rest.
By default the test is "> 100" that means all gene models that have ORF longer
than 100 Amino acids, will pass the test.
In the case of isoforms, the isoforms that do not pass the test are removed
(If all isoforms are removed, the gene is removed).
A gene with with any transcript having any CDS will be considered as non
coding gene and will not be removed.

## SYNOPSIS

```
agat_sp_filter_by_ORF_size.pl --gff infile.gff [ -o outfile ]
agat_sp_filter_by_ORF_size.pl -h
```

## OPTIONS

- **-g** or **--gff**

    Input GTF/GFF file.

- **-s** or **--size**

    ORF size to apply the test. Default 100.

- **-t** or **--test**

    Test to apply (> < = >= <=). If you us one of these two character >, <, please don't forget to quote you parameter like that "<=" otherwise your terminal will complain.
    By default it will be ">"


- **-v**

    Verbose. Useful for debugging purpose. Bolean

- **-o** or **--out** or **--output** or **--outfile**

    Output GFF file.  If no output file is specified, the output will be
    written to STDOUT.

- **-c** or **--config**

    String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any,
    otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
    The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

- **-h** or **--help**

    Display this helpful text.

