# agat_sp_add_splice_sites.pl

## DESCRIPTION

The script aims to add splice sites features (five_prime_cis_splice_site and three_prime_cis_splice_site) to gtf/gff file.
The splice sites are deduced from CDS features.

## SYNOPSIS

```
agat_sp_add_splice_sites.pl --gff infile --out outFile
agat_sp_add_splice_sites.pl --help
```

## OPTIONS

- **--gff**, **-f**, **--ref** or **-reffile**

    Input GTF/GFF file.

- **--out**, **--output** or **-o**

    Output file (default GFF3 - see config to modify output format).

- **-c** or **--config**

    String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any, 
    otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
    The --config option gives yo the possibility to use your own AGAT config file (located elsewhere or named differently).

- **--help** or **-h**

    Display this helpful text.

