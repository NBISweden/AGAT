# agat\_sp\_filter\_gene\_by\_length.pl

## DESCRIPTION

The script aims to filter level1 feature (e.g. gene, match, etc) by length.
It will create two files. one with the feature passing the length filter,
the other one with the remaining features.
If the level1 feature has exon features, the size is computed by concatenating
the exon together. If the level1 feature has several level2 features (e.g. mRNA)
we apply the test on the longest one (the longest concatenated exon set).

Some examples:
Select L1 feature shorter than 1000bp:
agat\_sp\_filter\_gene\_by\_length.pl --gff infile.gff  --size 1000 --test "<" -o result.gff
Select genes longer than 200bp:
agat\_sp\_filter\_gene\_by\_length.pl --gff infile.gff --size 200 --test ">" -o result.gff

## SYNOPSIS

```
agat_sp_filter_gene_by_length.pl --gff infile.gff --test ">=" --nb 10 [ --output outfile ]
agat_sp_filter_gene_by_length.pl --help
```

## OPTIONS

- **-f**, **--reffile**, **--gff**  or **-ref**

    Input GFF3 file that will be read

- **-s** or **--size**

    Integer - Gene size in pb \[Default 100\]

- **-t** or **--test**

    Test to apply (>, <, =, >= or <=). If you use one of these two characters >, <,
    please do not forget to quote your parameter like that "<=". Else your terminal will complain.
    \[Default "="\]

- **-o** or **--output**

    Output GFF file.  If no output file is specified, the output will be
    written to STDOUT.

- **-v**

    Verbose option for debugging purpose.

- **-c** or **--config**

    String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any,
    otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
    The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

- **-h** or **--help**

    Display this helpful text.

