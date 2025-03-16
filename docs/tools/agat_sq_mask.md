# agat_sq_mask.pl

## DESCRIPTION

The script masks (hard or soft) GFF-denoted segments out of a FASTA format file.
It needs 3 input parameters: a gff3 file, a fasta file, and a Mask method.
The result is written to the specified output file, or to STDOUT.

## SYNOPSIS

```
agat_sq_mask.pl -g infile.gff -f infile.fasta  (-hm or -sm) [ -o outfile ]
agat_sq_mask.pl --help
```

## OPTIONS

- **-g**, **--gff** or **-ref**

    Input GTF/GFF file.

- **-f**, **--fa**  or **--fasta**

    Input fasta file that will be masked

- **-sm**

    SoftMask option =>Sequences masked will be in lowercase

- **-hm**

    HardMask option => Sequences masked will be replaced by a character. By default the character used is 'n'. But you are allowed to speceify any character of your choice. To use 'z' instead of 'n' type: -hm z

- **-o** or **--output**

    Output GFF file.  If no output file is specified, the output will be
    written to STDOUT.

- **-c** or **--config**

    String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any,
    otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
    The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

- **-h** or **--help**

    Display this helpful text.

