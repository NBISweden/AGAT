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

- **-h** or **--help**

    Display this helpful text.

