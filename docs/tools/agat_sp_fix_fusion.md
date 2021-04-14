# agat\_sp\_fix\_fusion.pl

## DESCRIPTION

The script looks for other ORF in UTRs (UTR3 and UTR5) of each gene model described in the gff file.
Several ouput files will be written if you specify an output.
One will contain the gene not modified (intact), one the gene models fixed.

## SYNOPSIS

```
agat_sp_fix_fusion.pl --gff infile.gff --fasta genome.fa [ -o outfile ]
agat_sp_fix_fusion.pl --help
```

## OPTIONS

- **-gff**

    Input GTF/GFF file.

- **-fa** or **--fasta**

    Input fasta file.

- **--ct**, **--codon** or **--table**

    Codon table to use. \[default 1\]

- **-t** or **--threshold**

    This is the minimum length of new protein predicted that will be taken in account.
    By default this value is 100 AA.

- **-s** or **--stranded**

    By default we predict protein in UTR3 and UTR5 and in both direction. The fusion assumed can be between gene in same direction and in opposite direction.
    If RNAseq data used during the annotation was stranded, only fusion of close genes oriented in same direction are expected. In that case this option should be activated.
    When activated, we will try to predict protein in UTR3 and UTR5 only in the same orientation than the gene investigated.

- **-v** or **--verbose**

    Output verbose information.

- **-o** , **--output** , **--out** or **--outfile**

    Output GFF file.  If no output file is specified, the output will be
    written to STDOUT.

- **-h** or **--help**

    Display this helpful text.

