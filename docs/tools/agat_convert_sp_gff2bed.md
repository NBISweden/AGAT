# agat\_convert\_sp\_gff2bed.pl

## DESCRIPTION

The script aims to convert GTF/GXF file into bed file.
It will convert level2 features from gff (mRNA, transcripts) into bed features.
If  the selected level2 subfeatures (defaut: exon) exist, they will be reported
in the block fields (9-12th colum in bed).

Definintion of the bed format:: 

\## 1 chrom - The name of the chromosome (e.g. chr3, chrY, chr2\_random) or scaffold (e.g. scaffold10671).  
\## 2 chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.  
\## 3 chromEnd - The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99.  
\########### OPTIONAL fields ##########  
\## 4 name - Defines the name of the BED line. This label is displayed to the left of the BED line in the Genome Browser window when the track is open to full display mode or directly to the left of the item in pack mode.  
\## 5 score - A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray).  
\## 6 strand - Defines the strand - either '+' or '-'.  
\## 7 thickStart - The starting position at which the feature is drawn thickly  
\## 8 thickEnd - The ending position at which the feature is drawn thickly  
\## 9 itemRgb - An RGB value of the form R,G,B (e.g. 255,0,0). If the track line itemRgb attribute is set to "On", this RBG value will determine the display color of the data contained in this BED line. NOTE: It is recommended that a simple color scheme (eight colors or less) be used with this attribute to avoid overwhelming the color resources of the Genome Browser and your Internet browser.  
\## 10 blockCount - The number of blocks (exons) in the BED line.  
\## 11 blockSizes - A comma-separated list of the block sizes. The number of items in this list should correspond to blockCount.  
\## 12 blockStarts - A comma-separated list of block starts. All of the blockStart positions should be calculated relative to chromStart. The number of items in this list should correspond to blockCount.  

## SYNOPSIS

```
agat_convert_sp_gff2bed.pl --gff file.gff  [ -o outfile ]
agat_convert_sp_gff2bed.pl --help
```

## OPTIONS

- **--gff**

    Input GFF3 file that will be read

- **--sub**

    Define the subfeature (level3, e.g exon,cds,utr,etc...) to report as blocks in the bed output.
    Defaut: exon.

- **--outfile**, **--out**, **--output**, or **-o**

    File where will be written the result. If no output file is specified, the output will be written to STDOUT.

- **-h** or **--help**

    Display this helpful text.

