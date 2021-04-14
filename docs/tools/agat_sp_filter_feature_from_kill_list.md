## NAME

agat\_sp\_filter\_feature\_from\_kill\_list.pl

## DESCRIPTION

The script aims to remove features based on a kill list.
The default behaviour is to look at the features's ID. If the feature has an ID
(case insensitive) listed among the kill list it will be removed.
/!\\ Removing a level1 or level2 feature will automatically remove all linked subfeatures, and
removing all children of a feature will automatically remove this feature too.

## SYNOPSIS

```
agat_sp_filter_feature_from_kill_list.pl --gff infile.gff --kill_list file.txt  [ --output outfile ]
agat_sp_filter_feature_from_kill_list.pl --help
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

- **--kl** or **--kill\_list**

    Kill list. One value per line.

- **-a** or **--attribute**

    Attribute tag to specify the attribute to analyse. Case sensitive. Default: ID

- **-o** or **--output**

    Output GFF file.  If no output file is specified, the output will be
    written to STDOUT.

- **-v**

    Verbose option for debugging purpose.

- **-h** or **--help**

    Display this helpful text.

