# agat_convert_sp_gff2tsv.pl

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

- **-c** or **--config**

    String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any,
    otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
    The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

- **-h** or **--help**

    Display this helpful text.
