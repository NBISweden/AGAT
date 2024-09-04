# agat_sp_load_function_from_protein_align.pl

## DESCRIPTION

The script takes an annotation in gff format, a protein alignment in gff format and a protein fasta file as input. It checks if protein alignement overlap gene models, and will load the gene name and/or the function to the gene model according to the user requirements.
The script applies the following steps:
For each gene model structure it take the proteins aligned against, and sort them by an overlaping score. The best coming first.
Then it filters them by applying the overlaping score threshold.
1) If you activated the PE and the species filtering, we will try to find the best protein that follows the defined requirement.
2.1) If you activated the PE filtering or the precedent filtering (1) didn't succeed, we take the best protein according to the PE requirement.
2.2) If you activated the species filtering or the precedent filtering (1) didn't succeed, we take the best protein according to the list of prioritized species defined.
3) If no option or the precedent filtering (1,2.1,2.2)didn't succeed, the best protein will be selected.
You can flip the 2.1 and 2.2 test using the priority option.

## SYNOPSIS

```
agat_sp_load_function_from_protein_align.pl -a annotation.gff --pgff protein.gff --pfasta protein.fasta [ -o outfile ]
agat_sp_load_function_from_protein_align.pl --help
```

## OPTIONS

- **-a** or **--annotation**

    Input gtf/gff file of an annotation.

- **-pgff**

    Input gff file of aligned proteins.

- **-pfasta**

    Input protein fasta file where the extra information will be retrieved for each aligned protein.

- **-m** or **--method**

    Rule to apply to lift function when a protein map properly.
    1) replace  => replace or add the product and Name attribute's values.
    2) complete => add the product and Name attribute's values only if doesn't exist.
    3) add      => add the lfp_product and lfp_name attributes with the corresponding values

- **--value**, **--threshold** or **-t**

    Gene mapping percentage over which a gene must be reported. By default the value is 50.

- **-w**

    Compute the overlap score based on the whole annotation sequence. By default we use only the coding sequence part.

- **--pe**

    Protein existence value. We will take the best overlap score protein according to the PE expected
    1. Experimental evidence at protein level
    2. Experimental evidence at transcript level
    3. Protein inferred from homology
    4. Protein predicted
    5. Protein uncertain

- **--test**

    Test to apply (> < = >= <=). If you us one of these two character >, <, please don't forget to quote you parameter liket that "<=". Else your terminal will complain.

- **--sp**

    Species, between the set of the best protein aligned we try first to take the one that follow the species prioritization defined. There is a default one, but you can define you own (quoted and coma separated value)like that: "mus Musculus, Homo Sapiens" from the most important to the less important. In that case Mus will be taken first even if a better overlaping one exist for human.
    If none of them is found we take anyway the best overlapping one.

- **-p** or **--priority**

    By default the priority is PE test before species test when both are applied. You can flip these two test by activating this option like this: -p species

- **-v**

    Be verbose.

- **-o** , **--output** or **--out**

    Output GFF file.  If no output file is specified, the output will be
    written to STDOUT.

- **-c** or **--config**

    String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any,
    otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
    The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

- **-h** or **--help**

    Display this helpful text.

