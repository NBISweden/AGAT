# agat\_sp\_statistics.pl

## DESCRIPTION

The script provides exhaustive statistics of a gft/gff file.
/!\\ If you have isoforms in your file, even if correct, some values calculated
might sounds incoherent: e.g. total length mRNA can be superior than the genome size.
Because all isoforms length is added... It is why by default
we always compute the statistics twice when there are isoforms, once with the
isoforms, once without (In that case we keep the longest isoform per locus).


## SYNOPSIS

```
agat_sp_statistics.pl --gff file.gff  [ -o outfile ]
agat_sp_statistics.pl --help
```

## OPTIONS

- **--gff** or **-i**

    Input GTF/GFF file.

- **--gs**, **-f** or **-g**

    This option inform about the genome size in oder to compute more statistics. You can give the size in Nucleotide or directly the fasta file.

- **-d** or **-p**

    When this option is used, an histogram of distribution of the features will be printed in pdf files. (d means distribution, p means plot).

- **-v** or **--verbose**

    Verbose option. To modify verbosity. Default is 1. 0 is quiet, 2 and 3 are increasing verbosity.

- **--output** or **-o**

    File where will be written the result. If no output file is specified, the output will be written to STDOUT.

- **-c** or **--config**

    String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any,
    otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
    The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

- **-h** or **--help**

    Display this helpful text.

