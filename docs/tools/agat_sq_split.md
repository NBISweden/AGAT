# agat_sq_split.pl

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
- **--ft** or **--feature_type**
The top feature of the feature group. By default "gene".
- **-o** or **--output**

    STRING: Output file.  If no output file is specified, the output will be written to STDOUT. The result is in tabulate format.

- **-c** or **--config**

    String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any,
    otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
    The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

- **--help** or **-h**

    Display this helpful text.

