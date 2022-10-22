# agat\_convert\_sp\_gff2tsv.pl

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

- **-o** , **--output** , **--out** or **--outfile**

    Output GFF file.  If no output file is specified, the output will be
    written to STDOUT.

- **-h** or **--help**

    Display this helpful text.
