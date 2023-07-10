# agat\_sp\_fix\_cds\_frame.pl

## DESCRIPTION

This script aims to fix the cds phases.

## SYNOPSIS

```
agat_sp_fix_cds_frame.pl --gff infile.gff -f fasta [ -o outfile ]
agat_sp_fix_cds_frame.pl --help
```

## OPTIONS

- **-g**, **--gff** or **-ref**

    Input GTF/GFF file.

- **-fa** or **--fasta**

    Input fasta file.

- **-v** or **--verbose**

    Add verbosity.

- **-o** or **--output**

    Output GFF file.  If no output file is specified, the output will be
    written to STDOUT.

- **-c** or **--config**

    String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any,
    otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
    The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

- **-h** or **--help**

    Display this helpful text.

