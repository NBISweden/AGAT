# agat_convert_mfannot2gff.pl

## DESCRIPTION

Conversion utility for MFannot "masterfile" annotation produced by the MFannot
pipeline (http://megasun.bch.umontreal.ca/RNAweasel/). Reports GFF3 format.

## SYNOPSIS

```
agat_convert_mfannot2gff.pl -m <mfannot> -o <gff>
agat_convert_mfannot2gff.pl --help
```

## COPYRIGHT AND LICENSE

Copyright (C) 2015, Brandon Seah (kbseah@mpi-bremen.de)
... GPL-3 ...
modified by jacques dainat 2017-11

## OPTIONS

- **-m** or **-i** or **--mfannot**

    The mfannot input file

- **-g** or **-o** or **--gff**

    the gff output file

- **-c** or **--config**

    String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any,
    otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
    The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

- **--verbose** or **-v**

    Add verbosity
        
- **-h** or **--help**

    Display this helpful text.
