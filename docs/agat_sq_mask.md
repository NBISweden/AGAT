# NAME

agat\_sq\_mask.pl

# DESCRIPTION

The script masks GFF-denoted segments out of a FASTA format file.
This script masks (hard or soft) repeats among sequences.
It needs 3 input parameters: a gff3 file, a fasta file, and a Mask method.
The result is written to the specified output file, or to STDOUT.

# SYNOPSIS

```
agat_sq_mask.pl -g infile.gff -f infile.fasta  (-hm or -sm) [ -o outfile ]
agat_sq_mask.pl --help
```

# OPTIONS

- **-g**, **--gff** or **-ref**

    Input GTF/GFF file.

- **-f** or **--fasta**

    Input fasta file that will be masked

- **-sm**

    SoftMask option =>Sequences masked will be in lowercase

- **-hm**

    HardMask option => Sequences masked will be replaced by a character. By default the character used is 'n'. But you are allowed to speceify any character of your choice. To use 'z' instead of 'n' type: -hm z

- **-o** or **--output**

    Output GFF file.  If no output file is specified, the output will be
    written to STDOUT.

- **-h** or **--help**

    Display this helpful text.

