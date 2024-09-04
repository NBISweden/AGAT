# agat_sp_ensembl_output_style.pl

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
    to each other by a comon tag (by default locus_tag), this tag can be provided
    to parse the input file correctly.

- **-v**

    Verbose option to see the warning messages when parsing the gff file.

- **-o** or **--output**

    Output GFF file.  If no output file is specified, the output will be
    written to STDOUT.

- **-c** or **--config**

    String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any,
    otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
    The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

- **-h** or **--help**

    Display this helpful text.

