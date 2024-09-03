# agat_sp_filter_gene_by_intron_numbers.pl

## DESCRIPTION

The script aims to filter genes by intron numbers.
It will create two files. one with the genes passing the intron number filter,
the other one with the remaining genes.

Some examples:
Select intronless genes:
agat_sp_filter_gene_by_intron_numbers.pl --gff infile.gff -o result.gff
Select genes with more or equal 10 introns:
agat_sp_filter_gene_by_intron_numbers.pl --gff infile.gff --test ">=" --nb 10 [ --output outfile ]

## SYNOPSIS

```
agat_sp_filter_gene_by_intron_numbers.pl --gff infile.gff --test ">=" --nb 10 [ --output outfile ]
agat_sp_filter_gene_by_intron_numbers.pl --help
```

## OPTIONS

- **-f**, **--reffile**, **--gff**  or **-ref**

    Input GFF3 file that will be read

- **-n**,  **--nb** or **--number**

    Integer - Number of introns [Default 0]

- **-t** or **--test**
Test to apply (>, <, =, >= or <=). If you use one of these two characters >, <,
please do not forget to quote your parameter like that "<=". Else your terminal will complain.
[Default "="]
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

