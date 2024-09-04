# agat_sp_extract_attributes.pl

## DESCRIPTION

The script takes a gtf/gff file as input.
The script allows to extract choosen attributes of all or specific feature types.
The 9th column of a gff/gtf file contains a list of attributes.
An attribute (gff3) looks like that tag=value

## SYNOPSIS

```
agat_sp_extract_attributes.pl -gff file.gff  -att locus_tag,product,name -p level2,cds,exon [ -o outfile ]
agat_sp_extract_attributes.pl --help
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

- **--attribute**, **--att**, **-a**

    attribute tag. The value of the attribute tag specified will be extracted from the feature type specified by the option -p. List of attributes must be coma separated.

- **--merge** or **-m**

    By default the values of each attribute tag is writen in its dedicated file. To write the values of all tags in only one file use this option.

- **-d**

    By default when an attribute is not found for a feature, a dot (.) is reported. If you don't want anything to be printed in such case use this option.

- **-o** , **--output** , **--out** or **--outfile**

    Output GFF file.  If no output file is specified, the output will be
    written to STDOUT.

- **-c** or **--config**

    String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any,
    otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
    The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

- **-h** or **--help**

    Display this helpful text.

