# agat\_sp\_list\_short\_introns.pl

## DESCRIPTION

The script aims to list all the introns inferior to a certain size.
Introns are calculated on the fly from exons. (intron feature will not be used).

## SYNOPSIS

```
agat_sp_list_short_introns.pl --gff infile [ --out outFile ]
agat_sp_list_short_introns.pl --help
```

## OPTIONS

- **--gff**, **-f**, **--ref** or **-reffile**

    Input GTF/GFF file.

- **--size** or **-s**

    Minimum intron size accepted in nucleotide. All introns under this size will be reported.
    Default value = 10.

- **--out**, **--output** or **-o**

    Output gff3 file where the gene incriminated will be write.

- **-c** or **--config**

    String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any,
    otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
    The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

- **--help** or **-h**

    Display this helpful text.

