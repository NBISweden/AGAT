# agat\_sp\_sensitivity\_specificity.pl

## DESCRIPTION

The script aims to compute the Sensitivity and Specificity in order to assess the quality
of an annotation according to a reference (that is supposed to be true high-quality annotation).
The Sensitivity (Sn) is the proportion of true predictions compared to the total number of correct genes (including missed predictions)
Sn = TP / TP+FN
The Specificity (Sp) is the proportion of true predictions among all predicted genes (including incorrectly predicted ones)
Sp = TP / TP+FP

reference annotation:     -------------
prediction          :           ------------
                            FN     TP    FP    TN

Sensitivity and Specificity will be computed for each feature types.
(and computed independentaly if part of different Level2 type. i.e. exons Sn Sp
for tRNA will not be mixed up with the exon Sn Sp of mRNA exons)

## SYNOPSIS

```
agat_sp_sensitivity_specificity.pl --gff1 infile1.gff --gff2 infile2.gff  [ -o outfile ]
agat_sp_sensitivity_specificity.pl --help
```

## OPTIONS

- **-gff1**

    Input GTF/GFF file 1.

- **-gff2**

    Input GTF/GFF file 2.

- **-o** , **--output** , **--out** or **--outfile**

    Output GFF file.  If no output file is specified, the output will be
    written to STDOUT.

- **-v**

    Verbose option for debug purposes.

- **-h** or **--help**

    Display this helpful text.

