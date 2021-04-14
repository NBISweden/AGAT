## NAME

agat\_sp\_filter\_record\_by\_coordinates.pl

## DESCRIPTION

The script aims to filter the records to keep only those contained within coordinates
defined in an input csv file.
A record can be a feature or a set of features with part-of relationships.
By default we keep records overlapping the coordinates. The --exclude parameter
allows to keep only record fully contained within the coordinates.

! With default paramater, an exon out of the coordinates can be kept if the gene
it is part of is overlaping the coordinates.

## SYNOPSIS

```
agat_sp_filter_record_by_coordinates.pl --gff infile.gff --tsv coordinates.tsv [ --output outfile ]
agat_sp_filter_record_by_coordinates.pl --help
```

## OPTIONS

- **-i**, **--input**, **--gtf**  or **--gff**

    Input GTF/GFF file

- **-c**, **--coordinates**, **--tsv**, **-r** or **--ranges**

    String - tsv file containing the coordinates.
    Coordinates must be one per line.
    Each line must contain 3 fields separated by a tabulation.
    Field1 is the sequence id
    Field2 is the start coordinate (included)
    Field3 is the end coordinate (included)

- **-e** or **--exclude**

    Select only the features fully containined within the coordinates, exclude the overlapping
    ones.

- **-o** or **--output**

    Output GFF file.  If no output file is specified, the output will be
    written to STDOUT.

- **-v** or **--verbose**

    Verbosity.

- **-h** or **--help**

    Display this helpful text.

