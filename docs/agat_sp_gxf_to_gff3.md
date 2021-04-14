# /!\ replaced by [agat_convert_sp_gxf2gxf](https://github.com/NBISweden/AGAT/wiki/agat_convert_sp_gxf2gxf) since AGAT v0.3.0

# NAME

agat\_sp\_gxf\_to\_gff3.pl

# DESCRIPTION

This script fixes and/or standardizes a gtf/gff file into full sorted gff3.
It will be read by the Omniscient parser, that will first detect automtically
which GFF parser to use from bioperl (GFF1,GFF2,GFF3) to read the lines.
Then the Omniscient parser removes duplicate features, fixes duplicated IDs,
add missing ID and/or Parent attributes, deflates factorized attributes
(attributes with several parents are duplicated with uniq ID), add missing features
when possible (e.g add exon if only CDS described, add UTR if CDS and exon described),
fix feature locations (e.g check exon is embeded in th parent features mRNA, gene), etc...
All AGAT's scripts with the \_sp\_ prefix use the same parser, before to perform suplement tasks.
whith the script you can tuned the Omniscient parser behaviour. I.e you can decide
to merge loci that have an overlap at their CDS features ( Only one top feature
is kept (gene), and the mRNA features become isoforms). This is not activated by
default in case you are working on a prokaryote annotation that often have overlaping
loci.
The Omniscient parser defines relationship between features using 3 levels.
e.g Level1=gene; Level2=mRNA,tRNA; Level3=exon,cds,utr.
The feature type information is stored within the 3rd column of a GTF/GFF file.
Which level a feature is part of is crucial for the parser. This information
is stored by default in a json file coming with the tool. We have implemented the
most common feature types met in gff. If a feature type in your file is not imprelement
the parser will not handle it and inform you. You could easily inform the parser how
to handle it (level1,level2 or level3) by adding it in the corresponding file. How to
access the json files? Easy just use the --expose option and the json files will appear in
the workoing folder. If they are present, the Omniscient parser use the json files
from the working direcrtory by default.

Omniscient parser phylosophy:

```
Parse by Parent/child relationship
  ELSE Parse by a comomn tag  (an attribute value shared by feature that must be grouped together.
  By default we are using locus_tag and gene_id as locus tag, but you can specify the one of your choice
    ELSE Parse sequentially (features are grouped in a bucket, and the bucket change at each level2 feature met, and bucket(s) are linked to the first l1 top feature met)
```

# SYNOPSIS

```
agat_sp_gxf_to_gff3.pl -g infile.gff [ -o outfile ]
agat_sp_gxf_to_gff3.pl --help
```

# OPTIONS

- **-g**, **--gff** or **-ref**

    Input GTF/GFF file.

- **-c** or **--ct**

    When the features doesn't have Parent/ID relationships, the parser will try to group
    features using a common/shared attribute (i.e. a locus tag.). By default locus\_tag and gene\_id.
    You can provide another specific common/shared attribute using this option.

- **--efl** or **--expose**

    If you want to see, add or modified the feature relationships you will have to use this option.
    It will copy past in you working directory the json files used to define the relation between feature types and their level organisation.
    Typical level organisation: Level1 => gene; Level2 => mRNA; level3 => exon,cds,utrs
    If you get warning from the Omniscient parser that a feature relationship is not defined, you can provide information about it within the exposed json files.
    Indeed, if the json files exists in your working directory, they will be used by default.

- **--ml** or **--merge\_loci**

    Merge loci parameter, default deactivated. You turn on the parameter if you want to merge loci into one locus when they overlap.
    (at CDS level for mRNA, at exon level for other level2 features. Strand has to be the same). Prokaryote can have overlaping loci so it should not use it for prokaryote annotation.
    In eukaryote, loci rarely overlap. Overlaps could be due to error in the file, mRNA can be merged under the same parent gene if you acticate the option.

- **-v**

    Verbose option. To modify verbosity. Default is 1. 0 is quiet, 2 and 3 are increasing verbosity.

- **--nc** or **--no\_check**

    To deacticate all check that can be performed by the parser (e.g fixing UTR, exon, coordinates etc...)

- **-o** or **--output**

    Output GFF file.  If no output file is specified, the output will be
    written to STDOUT.

- **--gvi** or **--gff\_version\_input**

    If you don't want to use the autodection of the gff/gft version you give as input, you can force the tool to use the parser of the gff version you decide to use: 1,2,2.5 or 3. Remind: 2.5 is suposed to be gtf.

- **--gvo** or **--gff\_version\_output**

    If you don't want to use the autodection of the gff/gft version you give as input, you can force the tool to use the parser of the gff version you decide to use: 1,2,2.5 or 3. Remind: 2.5 is suposed to be gtf.

- **-h** or **--help**

    Display this helpful text.

