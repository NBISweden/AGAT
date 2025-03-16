# agat_sp_clipN_seqExtremities_and_fixCoordinates.pl

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

- **-f**, **--fa** or **--fasta**

    Input fasta file.

- **--of**

    Output fixed fasta file.  If no output file is specified, the output will be
    written to STDOUT.

- **--og**

    Output fixed GFF file.  If no output file is specified, the output will be
    written to STDOUT

- **-c** or **--config**

    String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any,
    otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
    The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

- **-h** or **--help**

    Display this helpful text.

