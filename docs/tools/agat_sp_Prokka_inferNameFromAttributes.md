# agat\_sp\_Prokka\_inferNameFromAttributes.pl

## DESCRIPTION

The script aims to fill a Name attribute based on &lt;gene> attribute in a prokka gff
annotation file. If no gene attribute is present it take if from the &lt;inference>
attribute.

## SYNOPSIS

```
agat_sp_Prokka_inferNameFromAttributes.pl -gff file.gff  [ -o outfile ]
agat_sp_Prokka_inferNameFromAttributes.pl --help
```

## OPTIONS

- **--gff** or **-f**

    Input GTF/GFF file.

- **--force**

    If Name attribute already exists, they will be replaced if a new one is found

- **-o** , **--output** , **--out** or **--outfile**

    Output GFF file.  If no output file is specified, the output will be
    written to STDOUT.

- **-c** or **--config**

    String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any,
    otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
    The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

- **-h** or **--help**

    Display this helpful text.

