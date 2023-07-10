# agat\_sp\_flag\_premature\_stop\_codons.pl

## DESCRIPTION

The script flags the mRNAs containing premature stop codons.
It will add the attribute "pseudo" and the value will be the positions of all premature stop codons.
Gene are flagged as pseudogene only if all the isoforms are pseudogenes. The attribute
will also be "pseudo" but will not contains any location.

## SYNOPSIS

```
agat_sp_flag_premature_stop_codons.pl --gff infile.gff --fasta infile.fa --out outfile
agat_sp_flag_premature_stop_codons.pl --help
```

## OPTIONS

- **--gff**, **--ref** or **-reffile**

    Input GTF/GFF file.

- **-f**, **--fa** or **--fasta**

    Imput fasta file.

- **--ct**, **--codon** or **--table**

    Codon table to use. \[default 1\]

- **--out**, **--output** or **-o**

    Output gff3 file where the result will be printed.

- **-c** or **--config**

    String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any,
    otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
    The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

- **--help** or **-h**

    Display this helpful text.

