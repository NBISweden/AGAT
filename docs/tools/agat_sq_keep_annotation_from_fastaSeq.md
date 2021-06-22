# agat\_sq\_filter\_feature\_from\_fasta.pl

## DESCRIPTION

This script is a kind of annotation filter by sequence name.
It goes through the gff annotation features and remove those that are not linked to a sequence from the fasta file provided.
The match between sequence name in the fasta file and the 1st column of the gff3 file is case sensitive.

## SYNOPSIS

```
agat_sq_filter_feature_from_fasta.pl --gff <gff_file.gff> --fasta <fasta_file.fa> [-o <output file>]
agat_sq_filter_feature_from_fasta.pl --help
```

## OPTIONS

- **--gff**, **--file** or **--input**

    STRING: Input GTF/GFF file.

- **-f** or **--fasta**

    STRING: fasta file.

- **-v** or **--verbose**

    For verbosity

- **-o** or **--output**

    STRING: Output file.  If no output file is specified, the output will be written to STDOUT. The result is in tabulate format.

- **--help** or **-h**

    Display this helpful text.
