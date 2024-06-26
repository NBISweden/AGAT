# agat\_sp\_filter\_feature\_by\_attribute\_value.pl

## DESCRIPTION

The script aims to filter features according to attribute value (9th column).
If the attribute tag is missing the feature will not be discarded.
If the attribute exists and the value pass the test, the feature is discarded.
Attribute are stored in the 9th column and have this shape: tag=value
/!\\ Removing a level1 or level2 feature will automatically remove all linked subfeatures, and
removing all children of a feature will automatically remove this feature too.

## SYNOPSIS

```
agat_sp_filter_feature_by_attribute_value.pl --gff infile.gff --value 1 -t "=" [ --output outfile ]
agat_sp_filter_feature_by_attribute_value.pl --help
```

## OPTIONS

- **-f**, **--reffile**, **--gff**  or **-ref**

    Input GFF3 file that will be read

- **-a** or **--attribute**

    Attribute tag to specify the attribute to analyse (attribute example: tag=value).

- **-p**,  **--type** or  **-l**

    primary tag option, case insensitive, list. Allow to specied the feature types that will be handled.
    You can specified a specific feature by given its primary tag name (column 3) as: cds, Gene, MrNa
    You can specify directly all the feature of a particular level:
          level2=mRNA,ncRNA,tRNA,etc
          level3=CDS,exon,UTR,etc
    By default all feature are taking into account. fill the option by the value "all" will have the same behaviour.

- **--value**

    Value(s) to check in the attribute. Case sensitive. List of values must be coma separated. 

- **--value\_insensitive**

    Bolean. Deactivated by default. When activated the values provided by the --value parameter are handled case insensitive.


- **-t** or **--test**

    Test to apply (> < = ! >= <=). default value "=". 
    If you use one of these two character >, <, please don't forget to quote the
    parameter like that "<=" otherwise your terminal will complain.
    Only = and ! tests can be used to compare string values.

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

