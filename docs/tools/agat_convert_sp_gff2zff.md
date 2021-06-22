# agat\_convert\_sp\_gff2zff.pl

## DESCRIPTION

The script converts GTF/GFF file into zff file a format used by the ab initio
tool SNAP. The script produces a .ann file containing the annotation and .dna
file containing the fasta file. The .ann and .dna are identicaly sorted by
sequence identifier (This is mandatory for usage with fathom).

## SYNOPSIS

```
agat_convert_sp_gff2zff.pl --gff file.gff  --fasta file.fasta [ -o outfile ]
agat_convert_sp_gff2zff.pl --help
```

## OPTIONS

- **--gff**

    Input GTF/GFF file

- **--fasta**

    Input fasta file

- **--outfile**, **--out**, **--output**, or **-o**

    File prefix where will be written the results (e.g. outfile.ann and outfile.dna).
    If no output file is specified, the output will be written to STDOUT.

- **-h** or **--help**

    Display this helpful text.

