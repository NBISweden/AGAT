# agat\_sq\_remove\_redundant\_entries.pl

## DESCRIPTION

The script remove redundant entries: same seq\_id,primary\_tag,start,stop,ID,Parent.
If ID and Parent attribute is not present, we do no remove the feature. If one of them
do not exists we use "" instead.

## SYNOPSIS

```
agat_sq_remove_redundant_entries.pl -i <input file> [-o <output file>]
agat_sq_remove_redundant_entries.pl --help
```

## OPTIONS

- **-i**, **--gff**, **--file** or **--input**

    STRING: Input GTF/GFF file.

- **-o** or **--output**

    STRING: Output file.  If no output file is specified, the output will be written to STDOUT. The result is in tabulate format.

- **--help** or **-h**

    Display this helpful text.

