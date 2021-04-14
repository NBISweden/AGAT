# NAME

agat\_sp\_flag\_short\_introns.pl

# DESCRIPTION

The script flags the short introns with the attribute &lt;pseudo>.
Is is usefull to avoid ERROR when submiting the data to EBI.
(Typical EBI error message: \*\*\*\*\*\*\*\*ERROR: Intron usually expected to be at least 10 nt long. Please check the accuracy)

# SYNOPSIS

```
agat_sp_flag_short_introns.pl --gff infile --out outfile
agat_sp_flag_short_introns.pl --help
```

# OPTIONS

- **--gff**, **-f**, **--ref** or **-reffile**

    Input GTF/GFF file.

- **--intron\_size** or **-i**

    Minimum intron size, default 10. All genes with an intron < of this size will be
    flagged with the pseudo attribute (the value will be the size of the smallest
    intron found within the incriminated gene)

- **--out**, **--output** or **-o**

    Output gff3 file where the result will be printed.

- **-v**

    Bolean. Verbose for debugging purpose.

- **--help** or **-h**

    Display this helpful text.

