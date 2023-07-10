# agat\_sp\_add\_start\_and\_stop.pl.pl

## DESCRIPTION

The script adds start and stop codons when a CDS feature exists.
The script looks at the nucleotide sequence and checks the presence of start and stop codons.
The script works even if the start or stop codon are split over several CDS features.

## SYNOPSIS

```
agat_sp_add_start_and_stop.pl.pl --gff infile.gff --fasta genome.fa --out outfile.gff
agat_sp_add_start_and_stop.pl.pl --help
```

## OPTIONS

- **--gff**, **-i** or **-g**

    Input GTF/GFF file.

- **--fasta**, **--fa** or **-f**

    Input fasta file. Needed to check that CDS sequences start by start codon and stop by stop codon.

- **--ct**, **--codon** or **--table**

    Codon table to use. \[default 1\]

- **--out**, **--output** or **-o**

    Output gff file updated

- **-v**

    Verbose for debugging purpose.

- **-c** or **--config**

    String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any,
    otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
    The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

- **--help** or **-h**

    Display this helpful text.

