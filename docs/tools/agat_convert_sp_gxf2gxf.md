# agat\_convert\_sp\_gxf2gxf.pl

## DESCRIPTION

This script fixes and/or standardizes any GTF/GFF file into full sorted GTF/GFF file.
It AGAT parser removes duplicate features, fixes duplicated IDs, adds missing ID and/or Parent attributes,
deflates factorized attributes (attributes with several parents are duplicated with uniq ID),
add missing features when possible (e.g. add exon if only CDS described, add UTR if CDS and exon described),
fix feature locations (e.g. check exon is embedded in the parent features mRNA, gene), etc...

All AGAT's scripts with the _sp_ prefix use the AGAT parser, before to perform any supplementary task.
So, it is not necessary to run this script prior the use of any other _sp_ script.

## SYNOPSIS

```
agat_convert_sp_gxf2gxf.pl -g infile.gff [ -o outfile ]
agat_convert_sp_gxf2gxf.pl --help
```

## OPTIONS

- **-g**, **--gff** or **-ref**

    String - Input GTF/GFF file. Compressed file with .gz extension is accepted.

- **-v**

    Integer - Verbose option. To modify verbosity. Default is 1. 0 is quiet, 2 and 3 are increasing verbosity.

- **-o** or **--output**

    String - Output GFF file.  If no output file is specified, the output will be
    written to STDOUT.

- **-h** or **--help**

    Boolean - Display this helpful text.
