# NAME

agat\_sq\_add\_attributes\_from\_tsv.pl

# DESCRIPTION

The script aims to add info from a tsv/csv file to the attributes of a gff file.
An attribute looks like that: tag=value1,value2
The first line of the tsv/csv must contains the headers, the other lines contain the values.
The header becomes the tag of the new attribute. If the tag already exists, the value will be added only
if the value does not already exists.
The first column does not become an attribute, indeed it must contain the feature ID
that will be used to know to which feature we will add the attributes.

\--- example ---

\* input.tsv:
ID	annot\_type1
gene1	anot\_x
cds1	anot\_y

\* gff:
chr1	irgsp	gene	1000	2000	.	+	.	ID=gene1
chr1	irgsp	CDS	2983	3268	.	+	.	ID=cds1

\* output.gff:
chr1	irgsp	gene	1000	2000	.	+	.	ID=gene1;annot\_type1=anot\_x
chr1	irgsp	CDS	2983	3268	.	+	.	ID=cds1;annot\_type1=anot\_y

# SYNOPSIS

```
agat_sq_add_attributes_from_tsv.pl --gff input.gff --tsv input.tsv [ -o output.gff3 ]
agat_sq_add_attributes_from_tsv.pl --help
```

# OPTIONS

- **--gff**

    STRING: Input GTF/GFF file.

- **--tsv**

    STRING: Input tsv file

- **--csv**

    BOLEAN: Inform the script that the tsv input file is actually a csv (coma-separated).

- **-v** or **--verbose**

    BOLEAN: Add verbosity

- **-o** or **--output**

    STRING: Output file. If no output file is specified, the output will be written
    to STDOUT. The result is in tabulate format.

- **--help** or **-h**

    Display this helpful text.

