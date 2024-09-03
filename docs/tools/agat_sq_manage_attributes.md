# agat_sp_manage_attributes.pl

## DESCRIPTION

The script removes choosen attributes of selected features. It can also create new
attribute with 'empty' value, or copy paste an existing attribute using a new specified tag.
Attribute in a gff file have this shape (2 attributes here): tag=value;tag=value and
are stored within the 9th column.

## SYNOPSIS

```
agat_sq_manage_attributes.pl --gff file.gff  --att locus_tag,product,name/NewName -p level2,cds,exon [ -o outfile ]
agat_sq_manage_attributes.pl --help
```

## OPTIONS

- **--gff** or **-f**

    Input GTF/GFF file.

- **-p**,  **--type** or  **-l**

    primary tag option, case insensitive, list. Allow to specied the feature types that will be handled.
    You can specified a specific feature by given its primary tag name (column 3) as: cds, Gene, MrNa
    You can specify directly all the feature of a particular level:
          level2=mRNA,ncRNA,tRNA,etc
          level3=CDS,exon,UTR,etc
    By default all feature are taking in account.

- **--tag**, **--att**

    Attributes with the tag specified will be removed from the feature type specified by the option p (primary tag). List of tag must be coma separated.
    /! You must use "" if name contains spaces.
    Instead to remove an attribute, you can replace its Tag by a new Tag using this formulation tagName/newTagName.
    To remove all attributes non mandatory (only ID and Parent are mandatory) you can use the option with &lt;all_attributes> parameter.

- **--add**

    Attribute with the tag specified will be added if doesn't exist. The value will be 'empty'.

- **--cp**

    When tags specied are with this form: tagName/newTagName.
    By using this &lt;cp> parameter, the attribute with the tag tagName will be duplicated
    with the new tag newTagName if no attribute with the tag newTagName already exits.

- **--overwrite**

    When using --add parameter, if an attribute with the specificed tag already exists, it will not be modified.
    When using --cp parameter, if an attribute with the specificed newTagName already exists, it will not be modified.
    So using the --overwrite parameter allows to overwrite the value of the existing attribute.

- **--value**

    String. When a value is provided the attribute is taken into account only if
    the attribute contains (or match) a specific value

- **--strategy**

    String. [Default equal]. Strategy to use when --value parameter is in use. Can be equal or match.
    Equal => the attribute value must be identical. Match => the attribute value must match

- **-o** , **--output** , **--out** or **--outfile**

    Output GFF file.  If no output file is specified, the output will be
    written to STDOUT.

- **-c** or **--config**

    String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any,
    otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
    The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

- **-h** or **--help**

    Display this helpful text.

