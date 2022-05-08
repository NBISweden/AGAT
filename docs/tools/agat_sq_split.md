# agat\_sq\_split.pl

## DESCRIPTION

split gff3 file into several files.
By default we create files containing 1000 genes and all sub-features associated.
GFF3 input file must be sequential.

## SYNOPSIS

```
agat_sq_split.pl --input <input file> -o <output file>
agat_sq_split.pl --help
```

## OPTIONS

- **--gff**, **--file** or **--input**

    STRING: Input GTF/GFF file.

- **-i** or **--interval**
Integer.  Number of group of feature to include in each file. 1000 by default.
- **--ft** or **--feature\_type**
The top feature of the feature group. By default "gene".
- **-o** or **--output**

    STRING: Output file.  If no output file is specified, the output will be written to STDOUT. The result is in tabulate format.

- **--help** or **-h**

    Display this helpful text.

