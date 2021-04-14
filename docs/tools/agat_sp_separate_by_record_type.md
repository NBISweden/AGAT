## NAME

agat\_sp\_separate\_by\_record\_type.pl

## DESCRIPTION

The script will separate the features from the gff input file into different files according to
the record type. A record represent all features linked collectively by Parent/ID relationships.
(e.g gene + mrna + exon + cds + utr of a locus).

a) When the record contains Level2 feature, the record type is the Level2 feature type (e.g tRNA,mRNA,ncRNA etc...)
b) Some features do not have children (top and standalone level1 features) e.g. location,region,chromosome.
In such case the record type is the level1 feature type.

## SYNOPSIS

```
agat_sp_separate_by_record_type.pl -g infile.gff [ -o outfolder ]
agat_sp_separate_by_record_type.pl --help
```

## OPTIONS

- **-g**, **--gff** or **-ref**

    Input GTF/GFF file.

- **-o** or **--output**

    Output folder.  If no output folder provided, the default name will be &lt;split\_result>.

- **-h** or **--help**

    Display this helpful text.

