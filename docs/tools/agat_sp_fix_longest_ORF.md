# agat\_sp\_fix\_longest\_ORF.pl

## DESCRIPTION

The script aims to fix the ORFs of gene models described in the gff file.
By fixing it means replacing the original ORF (defined by the cds)
when the longest predicted one within the mRNA is different. See the --model parameter
for more details about the different cases. Currently the tool does not perform
incomplete prediction (It always look for a start codon). It is consequently advised
to not use the model5 except if you understand what you do.
Several ouput files will be written if you specify an output.
One will contain the gene not modified (intact), one with the gene models fixed (modified),
one will both together (all).

## SYNOPSIS

```
agat_sp_fix_longest_ORF.pl -gff infile.gff --fasta genome.fa [ -o outfile ]
agat_sp_fix_longest_ORF.pl --help
```

## OPTIONS

- **--gff**

    Input GTF/GFF file.

- **-f**, **--fa** or **--fasta**

    Imput fasta file.

- **--ct**, **--codon** or **--table**

    Codon table to use. \[default 1\]

- **-m** or **--model**

    Kind of ORF Model you want to fix. By default all are used. To select specific models writte e.g --model 1,4

    Model1 = The original ORF is part of the new ORF; the new ORF is longer

    Model2 = The original ORF and the new one are different; the new one is longer, they do not overlap each other.

    Model3 = The original ORF and the new one are different; the new one is longer, they overlap each other.

    Model4 = The new ORF is shorter due to the presence of stop codon in the original ORF.

    Model5 = The new ORF is shorter but the original ORF has not premature stop codon.
             The shorter predicted ORF can be due to the fact that the original ORF does not start by a start codon,
    				 while we force here the prediction to have a start codon.
    				 A ORF wihtout start can be the fact of an incomplete or fragmented ORF:
    				 annotation tool didn't predict the start because:
    				 \* the start region is NNNN
    				 \* the start region is XXXX
    				 \* correct nucleotides but prediction tool did not annotate this part (e.g incomplete evidence in evidence-based prediction)

    Model6 = The ORF is same size but not correct frame (+1 or +2 bp gives a frame shift).

- **-s** or **--split**

    This option is usefull for Model2. Indeed when the prediction is non overlapping the original cds, it is possible to split the gene into two different genes. By default we don't split it.
    We keep the longest. If you want to split it type: -s

- **-o** , **--output** , **--out** or **--outfile**

    Output GFF file.  If no output file is specified, the output will be
    written to STDOUT.

- **-v**

    verbose mode. Default off. -v 1 minimum verbosity, -v 3 maximum verbosity

- **-h** or **--help**

    Display this helpful text.

