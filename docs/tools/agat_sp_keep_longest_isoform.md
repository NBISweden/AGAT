# agat\_sp\_keep\_longest\_isoform.pl

## DESCRIPTION

The script aims to filter isoforms when present. For a locus:
\- when all isoforms have CDS we keep the one with the longest CDS.
\- when some isoforms have CDS some others not, we keep the one with the longest CDS.
\- when none of the isoforms have CDS, we keep the one with the longest concatenated exons. 

## SYNOPSIS

```
agat_sp_keep_longest_isoform.pl -gff file.gff  [ -o outfile ]
agat_sp_keep_longest_isoform.pl --help
```

## OPTIONS

- **--gff** or **-f**

    GTF/GFF file.

- **--output** or **-o**

    File where will be written the result. If no output file is specified, the output will be written to STDOUT.

- **-c** or **--config**

    String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any,
    otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
    The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

- **-h** or **--help**

    Display this helpful text.

