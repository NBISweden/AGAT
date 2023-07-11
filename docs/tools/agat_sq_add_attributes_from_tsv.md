# agat\_sq\_add\_attributes\_from\_tsv.pl

## DESCRIPTION

The purpose of the script is to add information from a tsv/csv file to the attributes of a gff file (9th column).
e.g. an attribute looks like this in a GFF3 file: tag=value1,value2 
The first line of the tsv/csv file must contain the headers (corresponding to an attribute tag in the GFF/GTF file),
while the other lines contain the values (corresponding to an attribute value in the GFF/GTF file).
The first column is used to synchronize information between the tsv file and the GFF/GTF file. In other words, 
it's used to determine which feature we're going to add attributes to.
The other columns will be added as attribute in the GFF/GTF file. The header becomes the tag for the new attribute, 
and the value is that defined for the corresponding feature line. 
(If the tag already exists, we append the value only if the value doesn't already exist).

\--- example ---

\* input.tsv:  
```
ID	annot_type1  
gene1	annot_x  
cds1	annot_y  
```

\* input.gff:  
```
chr1	irgsp	gene	1000	2000	.	+	.	ID=gene1  
chr1	irgsp	CDS	2983	3268	.	+	.	ID=cds1  
```

\* output.gff: 
```
chr1	irgsp	gene	1000	2000	.	+	.	ID=gene1;annot_type1=annot_x  
chr1	irgsp	CDS	2983	3268	.	+	.	ID=cds1;annot_type1=annot_y  
```

\--- example2 ---

\* input.tsv:
```
gene_id	annot_type1
gene1	anot_x
cds1	anot_y
```

\* input gtf:
```
chr1	irgsp	gene	1000	2000	.	+	.	gene_id gene1;
chr1	irgsp	CDS	2983	3268	.	+	.	gene_id cds1;
```

\* output.gff:
```
chr1	irgsp	gene	1000	2000	.	+	.	gene_id gene1;annot_type1 anot_x
chr1	irgsp	CDS	2983	3268	.	+	.	gene_id=cds1;annot_type1=anot_y
```

## SYNOPSIS

```
agat_sq_add_attributes_from_tsv.pl --gff input.gff --tsv input.tsv [ -o output.gff3 ]
agat_sq_add_attributes_from_tsv.pl --help
```

## OPTIONS

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

- **-c** or **--config**

    String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any,
    otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
    The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

- **--help** or **-h**

    Display this helpful text.
