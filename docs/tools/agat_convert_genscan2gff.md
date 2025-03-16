# agat_convert_genscan2gff.pl

## DESCRIPTION

The script takes a genscan file as input, and will translate it in gff format.
The genscan format is described here: http://genome.crg.es/courses/Bioinformatics2003_genefinding/results/genscan.html
/! vvv Known problem vvv /!
You must have submited only DNA sequence, wihtout any header!!
Indeed the tool expects only DNA sequences and does not crash/warn if an header
is submited along the sequence.
e.g If you have an header ">seq" s-e-q are seen as the 3 first nucleotides of the sequence.
Then all prediction location are shifted accordingly.
(checked only on the online version http://argonaute.mit.edu/GENSCAN.html. I don't
know if there is the same pronlem elsewhere.)
/! ^^^ Known problem ^^^^ /!

## SYNOPSIS

```
agat_convert_genscan2gff.pl --genscan infile.bed [ -o outfile ]
agat_convert_genscan2gff.pl -h
```

## OPTIONS

- **--genscan** or **-g**

    Input bed file that will be convert.

- **--seqid**

    String - Sequence ID. [default: unknown]

- **--primary_tag**

    The primary_tag corresponf to the data type and is stored in 3rd field of a gff file.
    Example: gene,mRNA,CDS,etc.  [default: gene]

- **--verbose** or **-v**

    Add verbosity

- **-o**, **--output**, **--out**, **--outfile** or **--gff**

    Output GFF file. If no output file is specified, the output will be
    written to STDOUT.

- **-c** or **--config**

    String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any,
    otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
    The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

- **-h** or **--help**

    Display this helpful text.

