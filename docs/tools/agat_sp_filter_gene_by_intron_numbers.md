## NAME

agat\_sp\_filter\_gene\_by\_intron\_numbers.pl

## DESCRIPTION

The script aims to filter genes by intron numbers.
It will create two files. one with the genes passing the intron number filter,
the other one with the remaining genes.

Some examples:
Select intronless genes:
agat\_sp\_filter\_gene\_by\_intron\_numbers.pl --gff infile.gff -o result.gff
Select genes with more or equal 10 introns:
agat\_sp\_filter\_gene\_by\_intron\_numbers.pl --gff infile.gff --test ">=" --nb 10 \[ --output outfile \]

## SYNOPSIS

```
agat_sp_filter_gene_by_intron_numbers.pl --gff infile.gff --test ">=" --nb 10 [ --output outfile ]
agat_sp_filter_gene_by_intron_numbers.pl --help
```

## OPTIONS

- **-f**, **--reffile**, **--gff**  or **-ref**

    Input GFF3 file that will be read

- **-n**,  **--nb** or **--number**

    Integer - Number of introns \[Default 0\]

- **-t** or **--test**
Test to apply (>, <, =, >= or <=). If you use one of these two characters >, <,
please do not forget to quote your parameter like that "<=". Else your terminal will complain.
\[Default "="\]
- **-o** or **--output**

    Output GFF file.  If no output file is specified, the output will be
    written to STDOUT.

- **-v**

    Verbose option for debugging purpose.

- **-h** or **--help**

    Display this helpful text.

