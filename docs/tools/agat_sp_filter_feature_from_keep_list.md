# agat\_sp\_filter\_feature\_from\_keep\_list.pl

## DESCRIPTION

The script aims to keep records based on a keeplist.
The default behaviour is to look at the features's ID. If the feature has an ID
(case insensitive) listed among the keeplist it will be kept along with all
related features (the whole record is kept. A record repsent all features linked
 by relationship e.g. gene+transcript+exon+cds of a same locus).

## SYNOPSIS

```
agat_sp_filter_feature_from_keep_list.pl --gff infile.gff --keep_list file.txt  [ --output outfile ]
agat_sp_filter_feature_from_keep_list.pl --help
```

## OPTIONS

- **-f**, **--reffile**, **--gff**  or **-ref**

    Input GFF3 file that will be read

- **-p**,  **--type** or  **-l**

    primary tag option, case insensitive, list. Allow to specied the feature types that will be handled.
    You can specified a specific feature by given its primary tag name (column 3) as: cds, Gene, MrNa
    You can specify directly all the feature of a particular level:
          level2=mRNA,ncRNA,tRNA,etc
          level3=CDS,exon,UTR,etc
    By default all feature are taking into account. fill the option by the value "all" will have the same behaviour.

- **--kl** or **--keep\_list**

    Keep list. One value per line.

- **-a** or **--attribute**

    Attribute tag to specify the attribute to analyse. Case sensitive. Default: ID

- **-o** or **--output**

    Output GFF file.  If no output file is specified, the output will be
    written to STDOUT.

- **-v**

    Verbose option for debugging purpose.

- **-h** or **--help**

    Display this helpful text.

