# agat_sp_add_introns.pl

## DESCRIPTION

The script aims to add intron features to gtf/gff file without intron features.

## SYNOPSIS

```
agat_sp_add_introns.pl --gff infile --out outFile
agat_sp_add_introns.pl --help
```

## OPTIONS

- **--gff**, **-f**, **--ref** or **-reffile**

    Input GTF/GFF file.

- **--out**, **--output** or **-o**

    Output GFF3 file.

- **-c** or **--config**

    String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any,
    otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
    The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

- **--help** or **-h**

    Display this helpful text.

