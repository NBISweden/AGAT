# agat\_sp\_fix\_features\_locations\_duplicated.pl

## DESCRIPTION

The script aims to modify/remove features with duplicated locations. Even if it
not an error by itself in a gtf/gff file, it becomes problematic when submitting
the file to ENA (after convertion).

* Case1: When isoforms have identical exon structures AGAT removes duplicates by keeping the one with longest CDS;
* Case2: When l2 (e.g. mRNA) from different gene identifier have identical exon but no CDS at all (AGAT removes one duplicate);
* Case3: When l2 (e.g. mRNA) from different gene identifier have identical exon and CDS structures (AGAT removes duplicates by keeping the one with longest CDS);
* Case4: When l2 (e.g. mRNA) from different gene identifier have identical exon structures and different CDS structures (AGAT reshapes UTRs to modify mRNA and gene locations);
* Case5: When 2 genes have same locations while their exon/CDS locations are differents. In that case AGAT modifies the gene locations by clipping UTRs;


## SYNOPSIS

```
agat_sp_fix_features_locations_duplicated.pl --gff infile  [-o outfile]
agat_sp_fix_features_locations_duplicated.pl --help
```

## OPTIONS

- **-f**, **--file**, **--gff3** or **--gff**

    Input GTF/GFF file.

- **-o**, **--out**, **--output** or **--outfile**

    Output file. If none given, will be display in standard output.

- **--help** or **-h**

    Display this helpful text.

