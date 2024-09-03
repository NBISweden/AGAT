# agat_sp_filter_feature_from_kill_list.pl

## DESCRIPTION

The script aims to remove features based on a kill list.
The default behaviour is to look at the features's ID. If the feature has an ID
(case insensitive) listed among the kill list it will be removed.
/! Removing a level1 or level2 feature will automatically remove all linked subfeatures, and
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

- **--kl** or **--kill_list**

    Kill list. One value per line.

- **-a** or **--attribute**

    Attribute tag to specify the attribute to analyse. Case sensitive. Default: ID

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

