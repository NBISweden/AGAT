# agat\_sp\_add\_attribute\_shortest\_intron\_size.pl

## DESCRIPTION

The script add the attribute \<shortest_intron\> to each gene and rna, which will hold the size of the shortest intron in bp.

## SYNOPSIS

```
agat_sp_add_attribute_shortest_intron_size.pl --gff infile --out outfile
agat_sp_add_attribute_shortest_intron_size.pl --help
```

## OPTIONS

- **--gff**, **-f** or **--ref** 

    STRING: Input GTF/GFF file.

- **--out**, **--output** or **-o**

    STRING: Output gff3 file where the result will be printed.

- **-v**

    BOLEAN: Verbose for debugging purpose.

- **-c** or **--config**

    String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any,
    otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
    The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

- **--help** or **-h**

    Display this helpful text.
