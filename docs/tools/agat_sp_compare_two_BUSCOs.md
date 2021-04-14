# agat\_sp\_compare\_two\_BUSCOs.pl

## DESCRIPTION

The tool compares the results from two BUSCO runs (genome and proteome mode) in order to pinpoint the differences.
It compares the BUSCOs classification (complete,fragmented, duplicated) of the 1st run (genome mode)
against the classification found in the second run. It will report the results in txt files, and
extracts the complete,fragmented and duplicated annotated BUSCOs from the 1st run in gff files.
We add in the gff an attribute specifying the cases e.g. description=EOG090W00UK-complete2duplicated.
Where EOG090W00UK is the BUSCO name/label/group investigated, and complete2duplicated the case we found
(was complete in run1 and duplicated in run2).
By loading these gff tracks in a web browser and helped by other tracks (e.g the genome annotation/prediction)
can help to understand why the BUSCO have been classified differently from run1 to run2.
In other term it allows to catch potential problems in an annotation.
agat\_sp\_compare\_two\_BUSCOs.pl has been tested with results from BUSCO version 3 and 4.
/!\\ The tool expects a BUSCO run in genome mode as input folder 1 and a BUSCO run in proteins mode
as input folder 2. You can also decide to provide twice (--f1 --f2) the same BUSCO run in genome mode,
the tool will only extract the annotation of the complete,fragmented and duplicated annotated BUSCOs from the 1st run in gff.

## SYNOPSIS

```
agat_sp_compare_two_BUSCOs.pl --f1 <input busco folder1> --f2 <input busco folder2> [-o <output folder>]
agat_sp_compare_two_BUSCOs.pl --help
```

## OPTIONS

- **--f1**

    STRING: Input busco folder1

- **--f2**

    STRING: Input busco folder2

- **-v** or **--verbose**

    Integer: For displaying extra information use -v 1.
    For activating the verbosity in the omniscient parser use -v 66. (not recommended)

- **-o** or **--output**

    STRING: Output folder.

- **--help** or **-h**

    Display this helpful text.

