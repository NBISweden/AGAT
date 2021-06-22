# agat\_fix\_small\_exon\_from\_extremities.pl

## DESCRIPTION

The script aims to extend the small exons to make them longer.
When submitting annotation to ENA they expect exon size of 15 nt minimum.
Currently we extend only the exon from extremities, otherwise we risk to break the predicted ORF.
/!\\ When we extend an exon and the CDS has to be extended too (because is was a partial CDS), we exit;

## SYNOPSIS

```
agat_fix_small_exon_from_extremities.pl -gff infile.gff --fasta genome.fa [ -o outfile ]
agat_fix_small_exon_from_extremities.pl --help
```

## OPTIONS

- **-gff**

    Input GTF/GFF file.

- **-fa** or **--fasta**

    Genome fasta file
    The name of the fasta file containing the genome to work with.

- **--ct** or **--table** or **--codon**

    This option allows specifying the codon table to use - It expects an integer (1 by default = standard)

- **--size** or **-s**

    Minimum exon size accepted in nucleotide. All exon below this size will be extended to this size. Default value = 15.

- **-o** , **--output** , **--out** or **--outfile**

    Output GFF file.  If no output file is specified, the output will be
    written to STDOUT.

- **-v**

    Verbose option, make it easier to follow what is going on for debugging purpose.

- **-h** or **--help**

    Display this helpful text.

