# agat\_sp\_fix\_features\_locations\_duplicated.pl

## DESCRIPTION

The script aims to fix/remove feature with duplicated locations. Even if it
not an error by itself in a gtf/gff file, it becomes problematic when submitting
the file to ena (after convertion).

## SYNOPSIS

```
agat_sp_fix_features_locations_duplicated.pl --gff infile  [-o outfile]
agat_sp_fix_features_locations_duplicated.pl --help
```

## OPTIONS

- **-f**, **--file**, **--gff3** or **--gff**

    Input GTF/GFF file.

- **-o**, **--out**, **--output** or **--outfile**

    Output file. If none given, will be display in standard output.

- **--help** or **-h**

    Display this helpful text.

