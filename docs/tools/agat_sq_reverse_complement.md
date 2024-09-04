# agat_sq_reverse_complement.pl

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

- **-c** or **--config**

    String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any,
    otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
    The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

- **--help** or **-h**

    BOOLEAN: Display this helpful text.
