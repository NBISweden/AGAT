# agat\_sp\_fix\_features\_locations\_duplicated.pl

## DESCRIPTION

The script aims to modify/remove feature with duplicated locations. Even if it
not an error by itself in a gtf/gff file, it becomes problematic when submitting
the file to ENA (after convertion).
To modify locations, AGAT modify the UTRs (when available) by shortening them by 1 bp (and consequently the Parent features and the exons accordingly)

* Case1: When isoforms have identical exon structures, AGAT removes duplicates by keeping the one with longest CDS;
* Case2: When l2 (e.g. mRNA) from different gene identifier have identical exon but no CDS at all, AGAT removes one duplicate);
* Case3: When l2 (e.g. mRNA) from different gene identifier have identical exon and CDS structures, AGAT removes duplicates by keeping the one with longest CDS);
* Case4: When l2 (e.g. mRNA) from different gene identifier have identical exon structures and different CDS structures, AGAT reshapes UTRs to modify mRNA and gene locations);
* Case5: When l2 (e.g. mRNA) from different gene identifier overlap but have different exon structure. In that case AGAT modified the gene locations by clipping UTRs;


## SYNOPSIS

```
agat_sp_fix_features_locations_duplicated.pl --gff infile  [-o outfile]
agat_sp_fix_features_locations_duplicated.pl --help
```

## OPTIONS

- **-f**, **--file**, **--gff3** or **--gff**

    Input GTF/GFF file.

- **-m** or **--model**

    To select cases you want to fix. By default all are used.
    To select specific cases write e.g. --model 1,4,5

    Case1: When isoforms have identical exon structures AGAT removes duplicates by keeping the one with longest CDS;
    Case2: When l2 (e.g. mRNA) from different gene identifier have identical exon but no CDS at all (AGAT removes one duplicate);
    Case3: When l2 (e.g. mRNA) from different gene identifier have identical exon and CDS structures (AGAT removes duplicates by keeping the one with longest CDS);
    Case4: When l2 (e.g. mRNA) from different gene identifier have identical exon structures and different CDS structures (AGAT reshapes UTRs to modify mRNA and gene locations);
    Case5: When l2 (e.g. mRNA) from different gene identifier overlap but have different exon structure. In that case AGAT modified the gene locations by clipping UTRs;

- **-v** or **verbose**

    Add verbosity.

- **-o**, **--out**, **--output** or **--outfile**

    Output file. If none given, will be display in standard output.

- **--help** or **-h**

    Display this helpful text.
