# agat\_sp\_ensembl\_output\_style.pl

## DESCRIPTION

This script takes a normal gff3 annotation format file and convert it to gff3
like ensembl format.

## SYNOPSIS

```
agat_sp_ensembl_output_style.pl -g infile.gff [ -o outfile ]
agat_sp_ensembl_output_style.pl --help
```

## OPTIONS

- **-g**, **--gff** or **-ref**

    Input GTF/GFF file.

- **-c** or **--ct**

    When the gff file provided is not correcly formated and features are linked
    to each other by a comon tag (by default locus\_tag), this tag can be provided
    to parse the input file correctly.

- **-v**

    Verbose option to see the warning messages when parsing the gff file.

- **-o** or **--output**

    Output GFF file.  If no output file is specified, the output will be
    written to STDOUT.

- **-h** or **--help**

    Display this helpful text.

