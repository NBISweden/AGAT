# agat\_convert\_bed2gff.pl

## DESCRIPTION

The script takes a bed file as input, and will translate it in gff format.
The BED format is described [here](https://genome.ucsc.edu/FAQ/FAQformat.html##format1)
The script converts 0-based, half-open \[start-1, end) bed file to
1-based, closed \[start, end\] General Feature Format v3 (GFF3).

## SYNOPSIS

```
agat_convert_bed2gff.pl --bed infile.bed [ -o outfile ]
agat_convert_bed2gff.pl -h
```

## OPTIONS

- **--bed**

    Input bed file that will be converted.

- **--source**

    The source informs about the tool used to produce the data and is stored in 2nd field of a gff file.
    Example: Stringtie,Maker,Augustus,etc. \[default: data\]

- **--primary\_tag**

    The primary\_tag corresponds to the data type and is stored in 3rd field of a gff file.
    Example: gene,mRNA,CDS,etc.  \[default: gene\]

- **--inflate\_off**

    By default we inflate the block fields (blockCount, blockSizes, blockStarts) to create subfeatures
    of the main feature (primary\_tag). The type of subfeature created is based on the
    inflate\_type parameter. If you do not want this inflating behaviour you can deactivate it
    by using the --inflate\_off option.

- **--inflate\_type**

    Feature type (3rd column in gff) created when inflate parameter activated \[default: exon\].

- **--verbose**

    add verbosity

- **-o** , **--output** , **--out** , **--outfile** or **--gff**

    Output GFF file. If no output file is specified, the output will be
    written to STDOUT.

- **-h** or **--help**

    Display this helpful text.
