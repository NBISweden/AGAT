# agat_sq_rfam_analyzer.pl

## DESCRIPTION

The script allows to generate a tabulated format report of rfam-id annotated from a gff file
containing rfam results (type of the 3rd column must be ncRNA or nc_RNA - not case sensitive. And the 9th column must contain the rfam-id attribute).
    e.g:
ScG6Pog_82  Rfam  ncRNA 737595  737663  20.7  + 0 ID=RF00134_ScG6Pog_82_737595;Name=RF00134_ScG6Pog_82_737595;evalue=0.45;gc-content=0.28;model_end=1;model_start=1;rfam-acc=RF00134;rfam-id=snoZ196
ScG6Pog_82  Rfam  ncRNA 305023  305103  20.8  + 0 ID=RF00227_ScG6Pog_82_305023;Name=RF00227_ScG6Pog_82_305023;evalue=0.35;gc-content=0.31;model_end=1;model_start=1;rfam-acc=RF00227;rfam-id=FIE3

## SYNOPSIS

```
agat_sq_rfam_analyzer.pl -i <input file> [-g <integer or fasta> -o <output file>]
agat_sq_rfam_analyzer.pl --help
```

## OPTIONS

- **-i**, **--gff**, **--file** or **--input**

    STRING: Input GTF/GFF file(s). Several files can be processed at once: -i file1 -i file2

- **-g**, **--genome**

    That input is design to know the genome size in order to calculate the percentage of the genome represented by each kind of rfam-id.
    You can provide an INTEGER or the genome in fasta format. If you provide the fasta, the genome size will be calculated on the fly.

- **-o** or **--output**

    STRING: Output file.  If no output file is specified, the output will be written to STDOUT. The result is in tabulate format.

- **-c** or **--config**

    String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any,
    otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
    The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

- **--help** or **-h**

    Display this helpful text.

