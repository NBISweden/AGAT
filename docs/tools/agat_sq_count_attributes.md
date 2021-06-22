# agat\_sq\_count\_attributes.pl

# DESCRIPTION

The script count the number of a choosen attribute and also count the number of
unique value of this attribute.

# SYNOPSIS

```
agat_sq_count_attributes.pl --gff file.gff  --att gene_id [ -o outfile ]
agat_sq_count_attributes.pl --help
```

# OPTIONS

- **--gff** or **-f**

    Input GTF/GFF file.

- **--tag**, **--att**

    The name of the attribute that will be investigated.

- **-o** , **--output** , **--out** or **--outfile**

    Output GFF file.  If no output file is specified, the output will be
    written to STDOUT.

- **-h** or **--help**

    Display this helpful text.
