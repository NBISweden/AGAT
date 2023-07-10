# agat\_sp\_split\_by\_level2\_feature.pl

## DESCRIPTION

The script will split the gff input file into different files according to
the different Level2 feature that it contains.

## SYNOPSIS

```
agat_sp_split_by_level2_feature.pl -g infile.gff [ -o outfolder ]
agat_sp_split_by_level2_feature.pl --help
```

## OPTIONS

- **-g**, **--gff** or **-ref**

    Input GTF/GFF file.

- **-o** or **--output**

    Output folder.  If no output folder provided, the default name will be &lt;split\_result>.

- **-c** or **--config**

    String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any,
    otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
    The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

- **-h** or **--help**

    Display this helpful text.

