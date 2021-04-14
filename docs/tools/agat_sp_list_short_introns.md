# agat\_sp\_list\_short\_introns.pl

## DESCRIPTION

The script aims to list all the introns inferior to a certain size.
Introns are calculated on the fly from exons. (intron feature will not be used).

## SYNOPSIS

```
agat_sp_list_short_introns.pl --gff infile [ --out outFile ]
agat_sp_list_short_introns.pl --help
```

## OPTIONS

- **--gff**, **-f**, **--ref** or **-reffile**

    Input GTF/GFF file.

- **--size** or **-s**

    Minimum intron size accepted in nucleotide. All introns under this size will be reported.
    Default value = 10.

- **--out**, **--output** or **-o**

    Output gff3 file where the gene incriminated will be write.

- **--help** or **-h**

    Display this helpful text.

