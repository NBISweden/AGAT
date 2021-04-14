# NAME

agat\_sp\_functional\_statistics.pl

# DESCRIPTION

The script aims to summerize functional information stored in the file.

# SYNOPSIS

```
agat_sp_functional_statistics.pl --gff file.gff  [ -o outfile ]
agat_sp_functional_statistics.pl --help
```

# OPTIONS

- **--gff** or **-f**

    Input GTF/GFF file.

- **--gs** or **-g**

    This option inform about the genome size in oder to compute more statistics.
    You can give the size in Nucleotide or directly the fasta file.

- **--output** or **-o**

    File where will be written the result. If no output file is specified,
    the output will be written to STDOUT.

- **-h** or **--help**

    Display this helpful text.

