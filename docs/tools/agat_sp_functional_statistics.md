# agat_sp_functional_statistics.pl

## DESCRIPTION

The script aims to summerize functional information stored in the file.

## SYNOPSIS

```
agat_sp_functional_statistics.pl --gff file.gff  [ -o outfile ]
agat_sp_functional_statistics.pl --help
```

## OPTIONS

- **--gff** or **-f**

    Input GTF/GFF file.

- **--gs** or **-g**

    This option inform about the genome size in oder to compute more statistics.
    You can give the size in Nucleotide or directly the fasta file.

- **--output** or **-o**

    Folder where will be written the results. [Default output_functional_statistics]

- **-c** or **--config**

    String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any,
    otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
    The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

- **-h** or **--help**

    Display this helpful text.

