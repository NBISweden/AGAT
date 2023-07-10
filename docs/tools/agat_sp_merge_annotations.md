# agat\_sp\_merge\_annotations.pl

## DESCRIPTION

This script merge different gff annotation files in one.
It uses the AGAT parser that takes care of duplicated names and fixes other oddities met in those files.

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

- **-c** or **--config**

    String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any,
    otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
    The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

- **--help** or **-h**

    Display this helpful text.
