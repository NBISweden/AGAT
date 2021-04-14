## NAME

agat\_sq\_manage\_ID.pl

## DESCRIPTION

The script changes IDs to give uniq one and reflect the change in Parent attribute
of impacted features.

## SYNOPSIS

```
agat_sq_manage_ID.pl --gff <input file> [-o <output file>]
agat_sq_manage_ID.pl --help
```

## OPTIONS

- **--gff**, **--file** or **--input**

    STRING: Input GTF/GFF file.

- **--of**

    Output format, if no ouput format is given, the same as the input one detected will be used. Otherwise you can force to have a gff version 1 or 2 or 3 by giving the corresponding number.

- **-o** or **--output**

    STRING: Output file.  If no output file is specified, the output will be written to STDOUT. The result is in tabulate format.

- **--help** or **-h**

    Display this helpful text.

