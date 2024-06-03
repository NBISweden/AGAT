# agat\_sp\filter\_feature\_by\_attribute\_presence.pl

## DESCRIPTION

The script aims to filter features according to attribute presence (9th column).
If the attribute exists, the feature is discarded.
Attribute are stored in the 9th column and have this shape: tag=value
/!\\ Removing a level1 or level2 feature will automatically remove all linked subfeatures, and
removing all children of a feature will automatically remove this feature too.

## SYNOPSIS

```
agat_sp_filter_feature_by_attribute_presence.pl --gff infile.gff -a <tag> [ --output outfile ]
agat_sp_filter_feature_by_attribute_presence.pl --help
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

- **--attribute**, **--att**, **-a**

    String - Attributes tag specified will be used to filter the feature type (feature type can also be specified by the option -p). List of attribute tags must be coma separated.

- **--flip**

    BOLEAN - In order to flip the test and keep features that do have the attribute and filter those without

- **-o** or **--output**

    Output GFF file.  If no output file is specified, the output will be
    written to STDOUT.

- **-c** or **--config**

    String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any,
    otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
    The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

- **-h** or **--help**

    Display this helpful text.

