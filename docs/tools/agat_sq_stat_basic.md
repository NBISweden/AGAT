# agat\_sq\_stat\_basic.pl

## DESCRIPTION

The script aims to provide basic statistics of a gtf/gff file.

## SYNOPSIS

```
agat_sq_stat_basic.pl -i <input file> [-g <integer or fasta> -o <output file>]
agat_sq_stat_basic.pl --help
```

## OPTIONS

- **-i**, **--gff**, **--file** or **--input**

    STRING: Input GTF/GFF file. Several files can be processed at once: -i file1 -i file2

- **-g**, **--genome**

    That input is design to know the genome size in order to calculate the percentage of the genome represented by each kind of feature type.
    You can provide an INTEGER or the genome in fasta format. If you provide the fasta, the genome size will be calculated on the fly.

- **--inflate**

    Inflate the statistics taking into account feature with multi-parents.
    Indeed to avoid redundant information, some gff factorize identical features.
    e.g: one exon used in two different isoform will be defined only once, and will have multiple parent.
    By default the script count such feature only once. Using the inflate option allows
    to count the feature and its size as many time there are parents.

- **-o** or **--output**

    STRING: Output file. If no output file is specified, the output will be written to STDOUT. The result is in tabulate format.

- **--help** or **-h**

    Display this helpful text.

