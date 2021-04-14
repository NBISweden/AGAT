# NAME

agat\_sp\_kraken\_assess\_lift\_coverage.pl

# DESCRIPTION

The script takes as input gtf produced by Kraken (lift-over tool).
It will analyse the kraken\_mapped attributes to calculate the mapped percentage of each mRNA.
According to a threshold (0 by default), gene with a mapping percentage over that value will be reported.
A plot nammed geneMapped\_plot.pdf is performed to visualize the result.
/!\\ The script handles chimeric files (i.e containg gene part mapped on the template genome and others on the de-novo one)
/!\\/!\\ If the file is complete (containing kraken\_mapped="TRUE" and kraken\_mapped="FALSE" attributes),
the script calcul the real percentage lentgh that has been mapped.
Else the calcul is only based on feature with kraken\_mapped="TRUE" attributes.
So in this case the result most of time will be 100%.
/!\\/!\\/!\\ We met rare cases where Kraken mapped a feature to several locations of the de-novo genome.
As result we could end up with mapping over > 100%. We report them as 100% mapped in the plot
and a warning is raised to allow to check thoses cases.

# SYNOPSIS

```
agat_sp_kraken_assess_lift_coverage --gtf infile.gtf [ -o outfile ]
agat_sp_kraken_assess_lift_coverage --help
```

# OPTIONS

- **-gtf**

    Input gtf file produced by Kraken.

- **--threshold** or **-t**

    Gene mapping percentage over which a gene must be reported. By default the value is 0.

- **--verbose** or **-v**

    Verbose information.

- **-o** , **--output** , **--out** or **--outfile**

    Output GFF file.  If no output file is specified, the output will be
    written to STDOUT.

- **-h** or **--help**

    Display this helpful text.

