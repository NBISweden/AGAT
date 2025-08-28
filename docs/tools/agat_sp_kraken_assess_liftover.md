# agat_sp_kraken_assess_lift_coverage.pl

## DESCRIPTION

The script takes as input gtf produced by Kraken (lift-over tool).
It will analyse the kraken_mapped attributes to calculate the mapped percentage of each mRNA.
According to a threshold (0 by default), gene with a mapping percentage over that value will be reported.
A plot nammed geneMapped_plot.pdf is performed to visualize the result.
/! The script handles chimeric files (i.e containg gene part mapped on the template genome and others on the de-novo one)
/!/! If the file is complete (containing kraken_mapped="TRUE" and kraken_mapped="FALSE" attributes),
the script calcul the real percentage lentgh that has been mapped.
Else the calcul is only based on feature with kraken_mapped="TRUE" attributes.
So in this case the result most of time will be 100%.
/!/!/! We met rare cases where Kraken mapped a feature to several locations of the de-novo genome.
As result we could end up with mapping over > 100%. We report them as 100% mapped in the plot
and a warning is raised to allow to check thoses cases.

## SYNOPSIS

```
agat_sp_kraken_assess_lift_coverage --gtf infile.gtf [ -o outfile ]
agat_sp_kraken_assess_lift_coverage --help
```

## OPTIONS

- **-gtf**

    Input gtf file produced by Kraken.

- **--threshold** or **-t**

    Gene mapping percentage over which a gene must be reported. By default the value is 0.

- **--verbose** or **-v**

    Verbose information.

- **--quiet** or **-q**

    Run quietly. Sets verbosity to 0 and disables the progress bar.

- **-o** , **--output** , **--out** or **--outfile**

    Output GFF file.  If no output file is specified, the output will be
    written to STDOUT.

- **-c** or **--config**

    String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any,
    otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
    The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

- **-h** or **--help**

    Display this helpful text.

