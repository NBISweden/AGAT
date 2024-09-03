# agat_convert_embl2gff.pl

## DESCRIPTION

The script takes an EMBL file as input, and will translate it in gff format.

## SYNOPSIS

```
agat_converter_embl2gff.pl --embl infile.embl [ -o outfile ]
```

## OPTIONS

- **--embl**

    Input EMBL file that will be read

- **--emblmygff3**

    Bolean - Means that the EMBL flat file comes from the EMBLmyGFF3 software. 
    This is an EMBL format dedicated for submission and contains particularity to deal with.
    This parameter is needed to get a proper sequence id in the GFF3 from an embl made with EMBLmyGFF3.

- **--primary_tag**, **--pt**, **-t**

    List of "primary tag". Useful to discard or keep specific features.
    Multiple tags must be coma-separated.

- **-d**

    Bolean - Means that primary tags provided by the option "primary_tag" will be discarded.

- **-k**

    Bolean - Means that only primary tags provided by the option "primary_tag" will be kept.

- **-o**, **--output**, **--out**, **--outfile** or **--gff**

    Output GFF file. If no output file is specified, the output will be
    written to STDOUT.

- **-c** or **--config**

    String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any,
    otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
    The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

- **-h** or **--help**

    Display this helpful text.
