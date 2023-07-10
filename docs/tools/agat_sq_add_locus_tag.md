# agat\_sq\_add\_locus\_tag.pl

## DESCRIPTION

Add a shared locus tag per record. A record is all features linked by each other
by parent/children relationship (e.g Gene,mRNA,exon, CDS).

## SYNOPSIS

```
agat_sq_add_locus_tag.pl --gff <input file> [-o <output file>]
agat_sq_add_locus_tag.pl --help
```

## OPTIONS

- **--gff**, **--file** or **--input**

    STRING: Input GTF/GFF file.

- **-p**,  **--type** or  **-l**

    Primary tag option, case insensitive, list. Allow to specied the Level1 feature types that will be handled.
    By default all feature Level1 are taken into account.

- **--lo** or **--to**

    Locus tag output, by defaut it will be called locus\_tag, but using this option you can specied the name of this attribute.

- **--li** or **--ti**

    Tag input, by default the value of the locus tag attribute will be locusX where X is an incremented number.
    You can use the values of an existing attribute instead e.g the ID value: --li ID.

- **--of**

    Output format, if no ouput format is given, the same as the input one detected will be used.
    Otherwise you can force to have a gff version 1 or 2 or 3 by giving the corresponding number.

- **-o** or **--output**

    STRING: Output file.  If no output file is specified, the output will be written to STDOUT. The result is in tabulate format.

- **-q** or **--quiet**

    To remove verbosity.

- **-c** or **--config**

    String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any,
    otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
    The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

- **--help** or **-h**

    Display this helpful text.

