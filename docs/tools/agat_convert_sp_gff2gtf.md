# agat_convert_sp_gff2gtf.pl

## DESCRIPTION

The script aims to convert any GTF/GFF file into a proper GTF file.
Full information about the format can be found here: [https://agat.readthedocs.io/en/latest/gxf.html](https://agat.readthedocs.io/en/latest/gxf.html)
You can choose among 7 different GTF types (1, 2, 2.1, 2.2, 2.5, 3 or relax).
Depending the version selected the script will filter out the features that are not accepted.
For GTF2.5 and 3, every level1 feature (e.g nc_gene pseudogene) will be converted into
gene feature and every level2 feature (e.g mRNA ncRNA) will be converted into
transcript feature.
You can even produce a GFF-like GTF using the relax option. It allows to keep all
original feature types (3rd column). No modification will occur e.g. mRNA to transcript. 

To be fully GTF compliant all feature have a gene_id and a transcript_id attribute.
The gene_id is unique identifier for the genomic source of the transcript, which is
used to group transcripts into genes.
The transcript_id	is a unique identifier for the predicted transcript,
which is used to group features into transcripts.

## SYNOPSIS

```
agat_convert_sp_gff2gtf.pl --gff infile.gtf [ -o outfile ]
agat_convert_sp_gff2gtf -h
```

## OPTIONS

- **--gff** or **--in**

    Input GFF file that will be read

- **--gtf_version**
    version of the GTF output (1,2,2.1,2.2,2.5,3 or relax). Default 3.

    relax: all feature types are accepted.

    3: GTF3 (9 feature types accepted): gene, transcript, exon, CDS, Selenocysteine, start_codon, stop_codon, three_prime_utr and five_prime_utr

    2.5: GTF2.5 (8 feature types accepted): gene, transcript, exon, CDS, UTR, start_codon, stop_codon, Selenocysteine

    2.2: GTF2.2 (9 feature types accepted): CDS, start_codon, stop_codon, 5UTR, 3UTR, inter, inter_CNS, intron_CNS and exon

    2.1: GTF2.1 (6 feature types accepted): CDS, start_codon, stop_codon, exon, 5UTR, 3UTR

    2: GTF2 (4 feature types accepted): CDS, start_codon, stop_codon, exon

    1: GTF1 (5 feature types accepted): CDS, start_codon, stop_codon, exon, intron

- **-o**, **--output**, **--out** or  **--outfile**

    Output GTF file. If no output file is specified, the output will be
    written to STDOUT.

- **-c** or **--config**

    String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any,
    otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
    The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

- **--verbose** or **-v**

    Add verbosity
        
- **-h** or **--help**

    Display this helpful text.
