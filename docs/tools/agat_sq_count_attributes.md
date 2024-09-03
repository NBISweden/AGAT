# agat_sq_count_attributes.pl

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

- **-c** or **--config**

    String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any,
    otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
    The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

- **-h** or **--help**

    Display this helpful text.
