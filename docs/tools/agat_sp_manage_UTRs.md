# agat_sp_manage_UTRs.pl

## DESCRIPTION

Detect the genes containing too much UTR's exon according to a choosen threshold.
If no UTR option (3, 5, 3 and 5, both) is given the threshold will be not used.
option 3 and 5 together is different of "both". In the first case the gene is discarded if either the 3' or the 5' UTR contains more exon than the threshold given.
In the second case, will be discarded only the genes where the addition of UTR's exon of both side is over the threshold given.

## SYNOPSIS

```
agat_sp_manage_UTRs.pl --ref infile --three --five -p --out outFile
agat_sp_manage_UTRs.pl --help
```

## OPTIONS

- **--gff**, **--ref**, **--reffile** or **-f**

    Input GTF/GFF file.

- **-n**, **-t**, **--nb** or **--number**

    Threshold of exon's number of the UTR. Over or equal to this threshold, the UTR will be discarded. Default value is 5.

- **-3**, **--three** or **--tree_prime_utr**

    The threshold of the option &lt;n> will be applied on the 3'UTR.

- **-5**, **--five** or **--five_prime_utr**

    The threshold of the option &lt;n> will be applied on the 5'UTR.

- **-b**, **--both** or **--bs**

    The threshold of the option &lt;n> will be applied on genes where the number of UTR exon (3' and 5' additioned) is over it.

- **--p** or **--plot**

    Allows to create an histogram in pdf of UTR sizes distribution.

- **--out**, **--output** or **-o**

    Output gff3 file where the gene incriminated will be write.

- **-c** or **--config**

    String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any,
    otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
    The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

- **--help** or **-h**

    Display this helpful text.

