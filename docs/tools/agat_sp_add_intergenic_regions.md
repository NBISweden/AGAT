# agat\_sp\_add\_intergenic\_regions.pl

## DESCRIPTION

The script aims to add intergenic features (intergenic\_region) to gtf/gff file.
The intergenic regions are deduced from gene features (feature type gene from the 3rd column).

## SYNOPSIS

```
agat_sp_add_intergenic_regions.pl --gff infile --out outFile
agat_sp_add_intergenic_regions.pl --help
```

## OPTIONS

- **--gff**, **-f** or **--ref**

    Input GTF/GFF file.

- **--out**, **--output** or **-o**

    Output file (default GFF3 - see config to modify output format).

- **-c** or **--config**

    String - Input agat config file. By default AGAT takes as input agat\_config.yaml file from the working directory if any, 
    otherwise it takes the orignal agat\_config.yaml shipped with AGAT. To get the agat\_config.yaml locally type: "agat config --expose".
    The --config option gives yo the possibility to use your own AGAT config file (located elsewhere or named differently).

- **-v** or **--verbose**

    Add verbosity

- **-h** or **--help**

    Display this helpful text.

