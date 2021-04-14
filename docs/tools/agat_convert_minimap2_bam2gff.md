# NAME

agat\_convert\_sp\_minimap2\_bam2gff.pl

# DESCRIPTION

The script converts output from minimap2 (bam or sam) into gff file.
To get bam from minimap2 use the following command:
minimap2 -ax splice:hq genome.fa Asecodes\_parviclava.nucest.fa | samtools sort -O BAM -o output.bam
To use bam with this script you will need samtools in your path.

# SYNOPSIS

```
agat_convert_sp_minimap2_bam2gff.pl -i infile.bam [ -o outfile ]
agat_convert_sp_minimap2_bam2gff.pl -i infile.sam [ -o outfile ]
agat_convert_sp_minimap2_bam2gff.pl --help
```

# OPTIONS

if ( !GetOptions( 'i|input=s' => \\$opt\_in,

- **-i** or **--input**

    Input file in sam (.sam extension) or bam (.bam extension) format.

- **-b** or **--bam**

    To force to use the input file as sam file.

- **-s** or **--sam**

    To force to use the input file as sam file.

- **-o**, **--out** or **--output**

    Output GFF file.  If no output file is specified, the output will be
    written to STDOUT.

- **-h** or **--help**

    Display this helpful text.

