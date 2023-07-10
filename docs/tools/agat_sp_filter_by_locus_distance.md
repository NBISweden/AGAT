# agat\_sp\_filter\_by\_locus\_distance.pl

## DESCRIPTION

The script aims to remove or flag loci that are too close to each other.
Close loci are important to remove when training abinitio tools in order
to train intergenic region properly. Indeed if intergenic region
(surrouneded part of a locus) contain part of another locus,
the training on intergenic part will be biased.

## SYNOPSIS

```
agat_sp_filter_by_locus_distance.pl -gff infile.gff [ -o outfile ]
agat_sp_filter_by_locus_distance.pl --help
```

## OPTIONS

- **-gff**

    Input GTF/GFF file.

- **--dist** or **-d**

    The minimum inter-loci distance to allow.  No default (will not apply
    filter by default).

- **--add** or **--add\_flag**

    Instead of filter the result into two output files, write only one and add the flag &lt;low\_dist> in the gff.(tag = Lvalue or tag = Rvalue  where L is left and R right and the value is the distance with accordingle the left or right locus)

- **-o** , **--output** , **--out** or **--outfile**

    Output GFF file.  If no output file is specified, the output will be
    written to STDOUT.

- **-v**

    Verbose option, make it easier to follow what is going on for debugging purpose.

- **-c** or **--config**

    String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any,
    otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
    The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

- **-h** or **--help**

    Display this helpful text.

