# NAME

agat\_sp\_filter\_by\_mrnaBlastValue.pl

# DESCRIPTION

The script aims to remove from a gff file all the sequence that have a similarity
over THRESHOLD with another sequence (will keep only one).
This is typically useful when creating a list of mRNA to use to train abinitio gene finder.
A reciprocal blast of the sequences need to have been performed prior
to the use of this script in order to get the blastp input file.

# SYNOPSIS

```
agat_sp_filter_by_mrnaBlastValue.pl --gff infile.gff --blast blastfile --outfile outFile
agat_sp_filter_by_mrnaBlastValue.pl --help
```

# OPTIONS

- **--gff**

    Input GTF/GFF file.

- **--blast**

    The list of the all-vs-all blast file (outfmt 6, blastp)

- **--outfile**

    The name of the output file. By default the output is the standard output.

- **--help** or **-h**

    Display this helpful text.

