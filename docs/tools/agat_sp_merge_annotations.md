# agat\_sp\_merge\_annotations.pl

## DESCRIPTION

This script merge different gff annotation files in one.
It uses the Omniscient parser that takes care of duplicated names and fixes other oddities met in those files.

## SYNOPSIS

```
agat_sp_merge_annotations.pl --gff infile1 --gff infile2 --out outFile
agat_sp_merge_annotations.pl --help
```

## OPTIONS

- **--gff** or **-f**

    Input GTF/GFF file(s). You can specify as much file you want like so: -f file1 -f file2 -f file3

- **--out**, **--output** or **-o**

    Output gff3 file where the gene incriminated will be write.

- **--help** or **-h**

    Display this helpful text.

