# agat_sq_list_attributes.pl

## DESCRIPTION

The script take a gff3 file as input. -
The script give information about attribute tags used within you file.

## SYNOPSIS

```
agat_sq_list_attributes.pl -gff file.gff -p level2,cds,exon [ -o outfile ]
agat_sq_list_attributes.pl --help
```

## OPTIONS

- **--gff** or **-f**

    Input GTF/GFF file.

- **-p**,  **-t** or  **-l**

    primary tag option, case insensitive, list. Allow to specied the feature types that will be handled.
    You can specified a specific feature by given its primary tag name (column 3) as: cds, Gene, MrNa
    You can specify directly all the feature of a particular level:
          level2=mRNA,ncRNA,tRNA,etc
          level3=CDS,exon,UTR,etc
    By default all feature are taking in account. fill the option by the value "all" will have the same behaviour.

- **-o** , **--output** , **--out** or **--outfile**

    Output GFF file.  If no output file is specified, the output will be
    written to STDOUT.

- **-c** or **--config**

    String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any,
    otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
    The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

- **-h** or **--help**

    Display this helpful text.

