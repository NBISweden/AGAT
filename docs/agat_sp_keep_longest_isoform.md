# NAME

agat\_sp\_keep\_longest\_isoform.pl

# DESCRIPTION

The script aims to filter isoforms when present. For a locus:
\- when all isoforms have CDS we keep the one with the longest CDS.
\- when some isoforms have CDS some others not, we keep the one with the longest CDS.
\- when none of the isoforms have CDS, we keep the one with the longest concatenated exons. 

# SYNOPSIS

```
agat_sp_keep_longest_isoform.pl -gff file.gff  [ -o outfile ]
agat_sp_keep_longest_isoform.pl --help
```

# OPTIONS

- **--gff** or **-f**

    GTF/GFF file.

- **--output** or **-o**

    File where will be written the result. If no output file is specified, the output will be written to STDOUT.

- **-h** or **--help**

    Display this helpful text.

