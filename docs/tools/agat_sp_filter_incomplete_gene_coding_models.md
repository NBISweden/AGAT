# agat\_sp\_filter\_incomplete\_gene\_coding\_models.pl

## DESCRIPTION

The script aims to remove incomplete gene models. An incomplete gene coding model
is a gene coding with start and/or stop codon missing in its cds.
You can modify the behavior using the skip\_start\_check or skip\_stop\_check options.

## SYNOPSIS

```
agat_sp_filter_incomplete_gene_coding_models.pl --gff infile.gff --fasta genome.fa [ -o outfile ]
agat_sp_filter_incomplete_gene_coding_models.pl --help
```

## OPTIONS

- **-gff**

    Input GTF/GFF file.

- **-fa** or **--fasta**

    Genome fasta file.
    The name of the fasta file containing the genome to work with.

- **--ct** or **--table** or **--codon**

    This option allows specifying the codon table to use.
    It expects an integer \[default 1\]

- **--ad** or **--add\_flag**

    Instead of filter the result into two output files, write only one and add the flag &lt;incomplete> in the gff.(tag = inclomplete, value = 1, 2, 3.  1=start missing; 2=stop missing; 3=both)

- **--skip\_start\_check** or **--sstartc**

    Gene model must have a start codon. Activated by default.

- **--skip\_stop\_check** or **--sstopc**

    Gene model must have a stop codon. Activated by default.

- **-o** , **--output** , **--out** or **--outfile**

    Output GFF file.  If no output file is specified, the output will be
    written to STDOUT.

- **-v**

    Verbose option, make it easier to follow what is going on for debugging purpose.

- **-h** or **--help**

    Display this helpful text.

