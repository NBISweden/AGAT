# agat\_sp\_move\_attributes\_within\_records.pl

# DESCRIPTION

The script aims to keep move attributes within a record e.g. from Level1 to Level2 and/or Level3 features; and / or from Level2 to Level2 or Level3 features; and / or from Level3 to Level3 features.
Example of L1 feature: gene
Example of L2 featrue

# SYNOPSIS

```
agat_sp_move_attributes_within_records.pl --gff infile.gff --feature_copy mRNA  --feature_paste CDS --attribute Dbxref,Ontology [ --output outfile ]
agat_sp_move_attributes_within_records.pl --help
```

# OPTIONS

- **-f**, **--reffile**, **--gff**  or **-ref**

    Input GFF3 file that will be read

- **--feature\_copy** or **--fc**

    primary tag (feature type) option to list from which feature we will copy the attributes, case insensitive. 
    You can specified a feature (or a coma separated list) by giving its primary tag / feature type (column 3) value as: cds, Gene, MrNa, etc
    You can specify directly all the feature of a particular level:
          level2=mRNA,ncRNA,tRNA,etc
          level3=CDS,exon,UTR,etc
    By default all level2 feature are used. 

- **--feature\_paste** or **--fp**

    primary tag (feature type) option to list to which feature we will paste the attributes, case sensitive. 
    You can specified a feature (or a coma separated list) by giving its primary tag / feature type (column 3) value as: cds, Gene, MrNa, etc
    You can specify directly all the feature of a particular level:
          level2=mRNA,ncRNA,tRNA,etc
          level3=CDS,exon,UTR,etc
    By default all feature level3 are used. 

- **-a** or **--attribute**

    Attribute that will be copied and pasted. Case sensitive.
    You can specified an attribute (or a coma separated list) by giving its attribute tag value (column9) as: Ontology, Dbxref, etc
    Default: all\_attributes
    /!\\ &lt;all\_attributes> is a specific parameter meaning all the attributes will be use.

- **-o** or **--output**

    Output GFF file.  If no output file is specified, the output will be
    written to STDOUT.

- **-v**

    Verbose option for debugging purpose.

- **-c** or **--config**

    String - Input agat config file. By default AGAT takes as input agat\_config.yaml file from the working directory if any, 
    otherwise it takes the orignal agat\_config.yaml shipped with AGAT. To get the agat\_config.yaml locally type: "agat config --expose".
    The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

- **-h** or **--help**

    Display this helpful text.

