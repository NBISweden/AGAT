# agat_sp_filter_incomplete_gene_coding_models.pl

## DESCRIPTION

The script aims to remove incomplete gene models. An incomplete gene coding model
is a gene coding with start and/or stop codon missing in its cds.
You can modify the behavior using the skip_start_check or skip_stop_check options.

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
    It expects an integer [default 1]

- **--af** or **--add_flag**

    Instead of filter the result into two output files, write only one and add the flag &lt;incomplete> in the gff.(tag = inclomplete, value = 1, 2, 3.  1=start missing; 2=stop missing; 3=both)

- **--skip_start_check** or **--sstartc**

    Gene model must have a start codon. Activated by default.

- **--skip_stop_check** or **--sstopc**

    Gene model must have a stop codon. Activated by default.

- **-o** , **--output** , **--out** or **--outfile**

    Output GFF file.  If no output file is specified, the output will be
    written to STDOUT.

- **-v**

    Verbose option, make it easier to follow what is going on for debugging purpose.

- **-c** or **--config**

    String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any,
    otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
    The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

- **-h** or **--help**

    Display this helpful text.

