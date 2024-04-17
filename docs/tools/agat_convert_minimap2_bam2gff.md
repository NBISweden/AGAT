# agat\_convert\_minimap2\_bam2gff.pl

## DESCRIPTION

The script converts output from minimap2 (bam or sam) into GFF file.
To get bam from minimap2 use the following command:  

minimap2 -ax splice:hq genome.fa Asecodes\_parviclava.nucest.fa | samtools sort -O BAM -o output.bam  

To use bam with this script you will need samtools in your path.

## SYNOPSIS

```
agat_convert_minimap2_bam2gff.pl -i infile.bam [ -o outfile ]
agat_convert_minimap2_bam2gff.pl -i infile.sam [ -o outfile ]
agat_convert_minimap2_bam2gff.pl --help
```

## OPTIONS

- **-i** or **--input**

    Input file in sam (.sam extension) or bam (.bam extension) format.

- **-b** or **--bam**

    To force to use the input file as sam file.

- **-s** or **--sam**

    To force to use the input file as sam file.

- **-o**, **--out** or **--output**

    Output GFF file.  If no output file is specified, the output will be
    written to STDOUT.

- **-c** or **--config**

    String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any,
    otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
    The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

- **-h** or **--help**

    Display this helpful text.

