# agat\_sq\_add\_hash\_tag.pl

## DESCRIPTION

The script aims to introduce hash tag (####) into the file. It allows for some tools
using gff3 to handle independantly file chucks separated by the #### signal. Can make
them more efficient.

## SYNOPSIS

```
agat_sq_add_hash_tag.pl -i <input file> [-o <output file>]
agat_sq_add_hash_tag.pl --help
```

## OPTIONS

- **--gff**, **--file** or **--input**

    STRING: Input GTF/GFF file.

- **-i** or **--interval**

    Integer: 1 or 2. 1 will add #### after each new sequence (column1 of the gff), while 2 will add the ### after each group of feature (gene).
    By default the value is 1.

- **-o** or **--output**

    STRING: Output file.  If no output file is specified, the output will be written to STDOUT. The result is in tabulate format.

- **--help** or **-h**

    Display this helpful text.

