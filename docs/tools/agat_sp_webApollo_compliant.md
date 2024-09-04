# agat_sp_webApollo_compliant.pl

## DESCRIPTION

This script aim to remove useless/problematic information for webapollo,
change some featuree type to avoid problem whem loading them into webapollo,
and optimize some attribute for a nice displaying.

## SYNOPSIS

```
agat_sp_webApollo_compliant.pl -g infile.gff [ -o outfile ]
agat_sp_webApollo_compliant.pl --help
```

## OPTIONS

- **-g**, **--gff** or **-ref**

    Input GTF/GFF file.

- **-o** or **--output**

    Output GFF file.  If no output file is specified, the output will be
    written to STDOUT.

- **-c** or **--config**

    String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any,
    otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
    The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

- **-h** or **--help**

    Display this helpful text.

