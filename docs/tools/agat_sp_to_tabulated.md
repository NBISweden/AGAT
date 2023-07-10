# agat\_sp\_to\_tabulated.pl

## DESCRIPTION

The script aims to convert gtf/gff file into tabulated file.
Attribute's tags from the 9th column become column titles.

## SYNOPSIS

```
agat_sp_to_tabulated.pl -gff file.gff [ -o outfile ]
agat_sp_to_tabulated.pl --help
```

## OPTIONS

- **--gff** or **-f**

    Input GTF/GFF file.

- **-c** or **--ct**

    When the features doesn't have Parent/ID relationships, the parser will try to group
    features using a common/shared attribute (i.e. a locus tag.). By default locus\_tag and gene\_id.
    You can provide another specific common/shared attribute using this option.

- **--ml** or **--merge\_loci**

    Merge loci parameter, default deactivated. You turn on the parameter if you want to merge loci into one locus when they overlap.
    (at CDS level for mRNA, at exon level for other level2 features. Strand has to be the same). Prokaryote can have overlaping loci so it should not use it for prokaryote annotation.
    In eukaryote, loci rarely overlap. Overlaps could be due to error in the file, mRNA can be merged under the same parent gene if you acticate the option.

- **-o** , **--output** , **--out** or **--outfile**

    Output GFF file.  If no output file is specified, the output will be
    written to STDOUT.

- **-c** or **--config**

    String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any,
    otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
    The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

- **-h** or **--help**

    Display this helpful text.

