## NAME

gaas\_converter\_embl2gff.pl

## DESCRIPTION

The script takes an EMBL file as input, and will translate it in gff format.

## SYNOPSIS

```
gaas_converter_embl2gff.pl --embl infile.embl [ -o outfile ]
```

## OPTIONS

- **--embl**

    Input EMBL file that will be read

- **--primary\_tag**, **--pt**, **-t**

    List of "primary tag". Useful to discard or keep specific features.
    Multiple tags must be coma-separated.

- **-d**

    Means that primary tags provided by the option "primary\_tag" will be discarded.

- **-o**, **--output**, **--out**, **--outfile** or **--gff**

    Output GFF file. If no output file is specified, the output will be
    written to STDOUT.

- **-h** or **--help**

    Display this helpful text.

