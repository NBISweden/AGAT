# NAME

agat\_sp\_filter\_by\_ORF\_size.pl

# DESCRIPTION

The script reads a gff annotation file, and create two output files,
one contains the gene models with ORF passing the test, the other contains the rest.
By default the test is "> 100" that means all gene models that have ORF longer
than 100 Amino acids, will pass the test.

# SYNOPSIS

```
agat_sp_filter_by_ORF_size.pl --gff infile.gff [ -o outfile ]
agat_sp_filter_by_ORF_size.pl -h
```

# OPTIONS

- **-g** or **--gff**

    Input GTF/GFF file.

- **-s** or **--size**

    ORF size to apply the test. Default 100.

- **-t** or **--test**
Test to apply (> < = >= <=). If you us one of these two character >, <, please don't forget to quote you parameter liket that "<=". Else your terminal will complain.
By default it will be ">"
- **-v**

    Verbose. Useful for debugging purpose. Bolean

- **-o** or **--out** or **--output** or **--outfile**

    Output GFF file.  If no output file is specified, the output will be
    written to STDOUT.

- **-h** or **--help**

    Display this helpful text.

