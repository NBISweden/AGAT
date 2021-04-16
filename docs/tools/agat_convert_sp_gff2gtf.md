# agat\_convert\_sp\_gff2gtf.pl

## DESCRIPTION

The script aims to convert any GTF/GFF file into a proper GTF file.
Full information about the format can be found here: https://github.com/NBISweden/GAAS/blob/master/annotation/knowledge/gxf.md
You can choose among 6 different GTF types (1, 2, 2.1, 2.2, 2.5, 3).
Depending the version selected the script will filter out the features that are not accepted.
For GTF2.5 and 3, every level1 feature (e.g nc\_gene pseudogene) will be converted into
gene feature and every level2 feature (e.g mRNA ncRNA) will be converted into
transcript feature.
You can even produce a GFF-like GTF using the --relax option. It allows to keep all
original feature types (3rd column).

To be fully GTF compliant all feature have a gene\_id and a transcript\_id attribute.
The gene\_id	is unique identifier for the genomic source of the transcript, which is
used to group transcripts into genes.
The transcript\_id	is a unique identifier for the predicted transcript,
which is used to group features into transcripts.

## SYNOPSIS

```
agat_convert_sp_gff2gtf.pl --gff infile.gtf [ -o outfile ]
agat_convert_sp_gff2gtf -h
```

## OPTIONS

- **--gff** or **--in**

    Input GFF file that will be read

- **--gtf\_version**
version of the GTF output. Default 3 (for GTF3)

    GTF3 (9 feature types accepted): gene, transcript, exon, CDS, Selenocysteine, start\_codon, stop\_codon, three\_prime\_utr and five\_prime\_utr

    GTF2.5 (8 feature types accepted): gene, transcript, exon, CDS, UTR, start\_codon, stop\_codon, Selenocysteine

    GTF2.2 (9 feature types accepted): CDS, start\_codon, stop\_codon, 5UTR, 3UTR, inter, inter\_CNS, intron\_CNS and exon

    GTF2.1 (6 feature types accepted): CDS, start\_codon, stop\_codon, exon, 5UTR, 3UTR

    GTF2 (4 feature types accepted): CDS, start\_codon, stop\_codon, exon

    GTF1 (5 feature types accepted): 	CDS, start\_codon, stop\_codon, exon, intron

- **--relax**

    Relax option avoid to apply strict GTF format specification. All feature type will be kept.
    No modification e.g. mRNA to transcript.
    No filtering i.e. feature type not accepted by GTF format are kept.
    gene\_id and transcript\_id attributes will be added, and the attributes will follow the
    GTF formating.

- **-o** , **--output** , **--out** , **--outfile** or **--gtf**

    Output GTF file. If no output file is specified, the output will be
    written to STDOUT.

- **-h** or **--help**

    Display this helpful text.
