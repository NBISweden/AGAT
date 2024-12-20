# agat_sp_compare_two_annotations.pl

## DESCRIPTION

The script aims to compare two annotation of the same assembly. It provided
information about split/fusion of genes between the two annotations.  
The most common case are:  
1 => 0 ( gene uniq to file1)  
0 => 1 ( gene uniq to file2)  
1 => 1 ( 1 gene from file 1 overlaps only 1 gene from file2)  
1 => &lt;many> ( 1 gene from file 1 overlaps &lt;many> genes from file2) => split case (with file 1 as reference)  
&lt;many> => 1 ( &lt;many> genes from file 1 overlap only 1 gene from file2) => fusion case (with file 1 as reference)

Then you can get more complex cases:  
&lt;many> => &lt;many>  (&lt;many> genes from file 1 overlap &lt;many> genes from file2)

The script output a folder containing a report of number of different cases as well as a file
per case type listing per line the gene feature's ID involved in each case.

## SYNOPSIS

```
agat_sp_compare_two_annotations.pl --gff1 infile1.gff -gff2 infile2.gff [ -o outfile ]
agat_sp_compare_two_annotations.pl --help
```

## OPTIONS

- **--gff1**

    Input GTF/GFF file1.

- **--gff2**

    Input GTF/GFF file2.

- **-o** , **--output** , **--out** or **--outfile**

    Output folder.  It contains a report that resume the type and number of cases, as well as a file per case type 
    containing one case per line with the list of gene feature's ID (or other type of feature level1) from file1 then file2 separated by a |.

- **-v**

    Verbose option, make it easier to follow what is going on for debugging purpose.

- **-c** or **--config**

    String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any,
    otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
    The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

- **-h** or **--help**

    Display this helpful text.

