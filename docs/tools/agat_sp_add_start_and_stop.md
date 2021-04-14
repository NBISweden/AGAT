## NAME

agat\_sp\_add\_start\_and\_stop.pl.pl

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

- **--help** or **-h**

    Display this helpful text.

