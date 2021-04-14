# NAME

agat\_sq\_repeats\_analyzer.pl

# DESCRIPTION

The script allows to generate a tabulated format report of repeats annotated
from a gff file containing repeats (feature type must be match or protein\_match).

# SYNOPSIS

```
agat_sq_repeats_analyzer.pl -i <input file> [-g <integer or fasta> -o <output file>]
agat_sq_repeats_analyzer.pl --help
```

# OPTIONS

- **-i**, **--gff**, **--file** or **--input**

    STRING: Input GTF/GFF file(s). Several files can be processed at once: -i file1 -i file2

- **-g**, **--genome**

    That input is design to know the genome size in order to calculate the percentage of the genome represented by each kind of repeats.
    You can provide an INTEGER or the genome in fasta format. If you provide the fasta, the genome size will be calculated on the fly.

- **-o** or **--output**

    STRING: Output file.  If no output file is specified, the output will be written to STDOUT. The result is in tabulate format.

- **--help** or **-h**

    Display this helpful text.

