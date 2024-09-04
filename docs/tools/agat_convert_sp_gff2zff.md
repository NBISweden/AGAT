# agat_convert_sp_gff2zff.pl

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

- **-c** or **--config**

    String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any,
    otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
    The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

- **-h** or **--help**

    Display this helpful text.

