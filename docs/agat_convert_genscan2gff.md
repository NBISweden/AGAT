# NAME

agat\_convert\_genscan2gff.pl

# DESCRIPTION

The script takes a genscan file as input, and will translate it in gff format.
The genscan format is described here: http://genome.crg.es/courses/Bioinformatics2003\_genefinding/results/genscan.html
/!\\ vvv Known problem vvv /!\\
You must have submited only DNA sequence, wihtout any header!!
Indeed the tool expects only DNA sequences and does not crash/warn if an header
is submited along the sequence.
e.g If you have an header ">seq" s-e-q are seen as the 3 first nucleotides of the sequence.
Then all prediction location are shifted accordingly.
(checked only on the online version http://argonaute.mit.edu/GENSCAN.html. I don't
know if there is the same pronlem elsewhere.)
/!\\ ^^^ Known problem ^^^^ /!\\

# SYNOPSIS

```
agat_convert_genscan2gff.pl --genscan infile.bed [ -o outfile ]
agat_convert_genscan2gff.pl -h
```

# OPTIONS

- **--genscan** or **-g**

    Input bed file that will be convert.

- **--source**

    The source informs about the tool used to produce the data and is stored in 2nd field of a gff file.
    Example: Stringtie,Maker,Augustus,etc. \[default: data\]

- **--primary\_tag**

    The primary\_tag corresponf to the data type and is stored in 3rd field of a gff file.
    Example: gene,mRNA,CDS,etc.  \[default: gene\]

- **--inflate\_off**

    By default we inflate the block fields (blockCount, blockSizes, blockStarts) to create subfeatures
    of the main feature (primary\_tag). Type of subfeature created based on the
    inflate\_type parameter. If you don't want this inflating behaviour you can deactivate it
    by using the option --inflate\_off.

- **--inflate\_type**

    Feature type (3rd column in gff) created when inflate parameter activated \[default: exon\].

- **--verbose**

    add verbosity

- **-o** , **--output** , **--out** , **--outfile** or **--gff**

    Output GFF file. If no output file is specified, the output will be
    written to STDOUT.

- **-h** or **--help**

    Display this helpful text.

