# agat\_sq\_reverse\_complement.pl

## DESCRIPTION

This script will reverse complement the annotation of all annotation from the gff that are hold by sequences described in the fasta file.
The match between sequence name in the fasta file and the 1st column of the gff3 file is case sensitive.

## SYNOPSIS

```
agat_sq_reverse_complement.pl --gff <gff_file.gff> --fasta <fasta_file.fa> [-o <output file>]
agat_sq_reverse_complement.pl --help
```

## OPTIONS

- **--gff**, **--file** or **--input**

    STRING: Input GTF/GFF file.

- **-f**, **--fasta**

    STRING: fasta file.

- **-v**, **--verbose**

    BOOLEAN: For verbosity.

- **-o** or **--output**

    STRING: Output file.  If no output file is specified, the output will be written to STDOUT.

- **--help** or **-h**

    BOOLEAN: Display this helpful text.
