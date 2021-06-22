# agat\_sp\_clipN\_seqExtremities\_and\_fixCoordinates.pl

## DESCRIPTION

The script aims to clip the N's extremities of the sequences.
The annotation from the sequence clipped are modified accrodingly to stay consistent

## SYNOPSIS

```
agat_sp_clipN_seqExtremities_and_fixCoordinates.pl -g infile.gff -f infile.fasta  [ -o outfile ]
agat_sp_clipN_seqExtremities_and_fixCoordinates.pl --help
```

## OPTIONS

- **-g**, **--gff** or **-ref**

    Input GTF/GFF file.

- **-f** or **--fasta**

    Input fasta file.

- **--of**

    Output fixed fasta file.  If no output file is specified, the output will be
    written to STDOUT.

- **--og**

    Output fixed GFF file.  If no output file is specified, the output will be
    written to STDOUT

- **-h** or **--help**

    Display this helpful text.

