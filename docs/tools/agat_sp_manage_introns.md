# agat\_sp\_manage\_introns.pl

## DESCRIPTION

The script provides information about introns (longest, shortest size mean ...) using the statistic method,
then plot all the intron size values to get an overview of the introns size distribution.
It gives you as well the value of the longest intron after removing X percent(s) of the longest (removing potential biais / false positive).

## SYNOPSIS

```
agat_sp_manage_introns.pl --gff infile [--out outFile]
agat_sp_manage_introns.pl --help
```

## OPTIONS

- **--gff**, **-f**, **--ref** or **-reffile**

    Input GTF/GFF file. You can use several input files by doing: -f file1 -f file2 -f file3

- **-w**, **--window**, **--break**, **--breaks** or **-b**

    It the number of break used within the histogram plot. By default it's 1000. You can modify the value to get something more or less precise.

- **-x**, **--p**

    Allows to modify the X values to calculate the percentage of the longest introns to remove. By default the value is 1 (We remove 1 percent of the longest).

- **--plot**

    Allows to create an histogram in pdf of intron sizes distribution.

- **--out**, **--output** or **-o**

    Output gff3 file where the gene incriminated will be write.

- **--help** or **-h**

    Display this helpful text.

