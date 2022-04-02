# What can AGAT do for you?

AGAT has the power to check, fix, pad missing information (features/attributes) of any kind of GTF and GFF to create complete, sorted and standardised gff3 format.
The toolkit comes with an exhaustive list of tools allowing to perform almost everything you might want to achieve ^^

## How does AGAT work?

All tools taking GFF/GTF as input can be divided in two groups: \_sp\_ and \_sq\_.

### Tools with \_sp\_ prefix

\_sp\_ stands for SLURP. Those tools will charge the file in memory Omniscient data structure. It has a memory cost but makes life smoother. Indeed, it allows to perform complicated tasks in a more time efficient way ( Any features can be accessed at any time by AGAT).
Moreover, it allows to fix all potential errors in the limit of the possibilities given by the format itself.
See the Omniscient section for more information about it.  

### with \_sq\_ prefix

 \_sq\_ stands for SEQUENTIAL. Those tools will read and process GFF/GTF files from the top to the bottom, line by line, performing tasks on the fly. This is memory efficient but the sanity check of the file is minimum. Those tools are not intended to perform complex tasks.

## Omniscient / parsing performed by \_sp\_ prefix tools / Standardisation for a full GFF3 compliant to any tool

The first step of AGAT' tools with the \_sp\_ prefix of is to fix the file to standardize it. (e.g. a file containing only exon will be modified to create mRNA and gene features). To perform this task AGAT parses and slurps the entire data into a data structure called Omniscient.
Below you will find more information about peculiarity of the Omniscient structure, and the parsing approach used.

### What performs the Omniscient parser

* It creates missing parental features. (e.g if a level2 or level3 feature do not have parental feature(s) we create the missing level2 and/or level1 feature(s)).    
* It creates missing mandatory attributes (ID and/or Parent).  
* It fixes identifier to be uniq.  
* It removes duplicated features (same position, same ID, same Parent).  
* It expands level3 features sharing multiple parents (e.g  if one exon has list of multiple parent mRNA in its Parent attribute, one exon per parent with uniq ID will be created.  
* It fixes feature location errors (e.g an mRNA spanning over its gene location, we fix the gene location).  
* It adds UTR if possible (CDS and exon present).  
* It adds exon if possible (CDS has to be present).  
* It groups features together (if related features are spread at different places in the file).  


### Omniscient data structure

The method create a hash structure containing all the data in memory. We call it OMNISCIENT. The OMNISCIENT structure is a three levels structure:
```
$omniscient{level1}{tag_l1}{level1_id} = feature <= tag could be gene, match  
$omniscient{level2}{tag_l2}{idY} = @featureListL2 <= tag could be mRNA,rRNA,tRNA,etc. idY is a level1_id (know as Parent attribute within the level2 feature). The @featureList is a list to be able to manage isoform cases.  
$omniscient{level3}{tag_l3}{idZ} =  @featureListL3 <= tag could be exon,cds,utr3,utr5,etc. idZ is the ID of a level2 feature (know as Parent attribute within the level3 feature). The @featureList is a list to be able to put all the feature of a same tag together.  
```

### How does the Omniscient parser work

To resume by priority of way to parse: **Parent/child relationship > common attribute/tag > sequential.**  
The parser may used only one or a mix of these approaches according of the peculiarity of the gtf/gff file you provide.
If you need to use the `--ct` option you will have to process the file `agat_convert_sp_gxf2gxf.pl` first  before running any other tool.

  1) Parsing approach 1: by Parent/child relationship  

Example of Parent/ID relationship used by the GFF format:

    chr12	HAVANA	gene	100	500	.	+	.	ID=gene1
    chr12	HAVANA	transcript	100	500	.	+	.	ID=transcript1;Parent=gene1
    chr12	HAVANA	exon	100	500	.	+	.	ID=exon1;Parent=transcript1  
    chr12	HAVANA	CDS	100	500	.	+	0	ID=cds-1;Parent=transcript1

Example of gene_id/transcript_id relationship used by the GTF format:

    chr12	HAVANA	gene	100	500	.	+	.	gene_id "gene1";
    chr12	HAVANA	transcript	100	500	.	+	.	gene_id "gene1"; transcript_id "transcript1";
    chr12	HAVANA	exon	100	500	.	+	.	gene_id "gene1"; transcript_id "transcript1"; exon_id=exon1;
    chr12	HAVANA	CDS	100	500	.	+	0	gene_id "gene1"; transcript_id "transcript1"; cds_id=cds-1;

  2) ELSE Parsing approach 2: by a common attribute/tag  

  a common attribute (or common tag) is an attribute value shared by feature that must be grouped together. AGAT uses default attributes (`gene_id` and `locus_tag`) displayed in the log but can be set by the user using the `--ct` parameter).  

Example of relationship made using a commong tag (here locus_tag):

    chr12	HAVANA	gene	100	500	.	+	.	locus_tag="gene1"
    chr12	HAVANA	transcript	100	500	.	+	.	locus_tag="gene1";ID="transcript1"
    chr12	HAVANA	exon	100	500	.	+	.	locus_tag="gene1";ID=exon1;
    chr12	HAVANA	CDS	100	500	.	+	0	locus_tag="gene1";ID=cds-1;

  3) ELSE Parsing approach 3: sequentially.

  Reading from top to th ebotoom of the file, level3 features (e.g. exon, CDS, UTR) are attached to the last level2 feature (e.g. mRNA) met, and level2 feature are attached to the last L1 feature (e.g. gene) met. To see the list of features of each level see the corresponding json file (In the share folder in the github repo or using `agat_convert_sp_gxf2gxf.pl --expose`).

Example of relationship made sequentially:

    chr12	HAVANA	gene	100	500	.	+	.	ID="aaa"
    chr12	HAVANA	transcript	100	500	.	+	.	ID="bbb"
    chr12	HAVANA	exon	100	500	.	+	.	ID="ccc"
    chr12	HAVANA	CDS	100	500	.	+	0	ID="ddd"
    chr12	HAVANA	gene	1000	5000	.	+	.	ID="xxx"
    chr12	HAVANA	transcript	1000	5000	.	+	.	ID="yyy"
    chr12	HAVANA	exon	1000	5000	.	+	.	ID="zzz"
    chr12	HAVANA	CDS	1000	5000	.	+	0	ID="www"


### Particular case

Below you will find more information about peculiarity of the Omniscient structure, and the parsing approach used.

#### Level1 feature type missing and no Parent/gene_id

If you have isoforms (for Eukaryote organism) in your files and the `common attribute` used is not set properly you can end up with isoforms having independent parent gene features. See below for more details.

Here an example of three transcripts from two different genes (isoforms exist):

    chr12   HAVANA  transcript  100 500 .   +   .   ID="bbb";common_tag="gene1";transcript_id="transcript1";gene_info="gene1"
    chr12   HAVANA  exon    100 500 .   +   .   ID="ccc";common_tag="gene1"
    chr12   HAVANA  CDS 100 500 .   +   0   ID="ddd";common_tag="gene1"
    chr12   HAVANA  transcript  100 600 .   +   .   ID="bbb2";common_tag="gene1";transcript_id="transcript2";gene_info="gene1"
    chr12   HAVANA  exon    100 600 .   +   .   ID="ccc2";common_tag="gene1"
    chr12   HAVANA  CDS 100 600 .   +   0   ID="ddd2";common_tag="gene1"
    chr12   HAVANA  transcript  1000    5000    .   +   .   ID="yyy";common_tag="gene2";transcript_id="transcript3";gene_info="gene2"
    chr12   HAVANA  exon    1000    5000    .   +   .   ID="zzz";common_tag="gene2"
    chr12   HAVANA  CDS 1000    5000    .   +   0   ID="www";common_tag="gene2"


/!\ AGAT will create a gene feature by transcript. This is WRONG here because we have ISOFORMS (transcript1 and transcript2 are from the same gene) that must be attached to the same gene:  

`agat_convert_sp_gxf2gxf.pl --gff input.gxf`

    chr12   HAVANA  gene    100 500 .   +   .   ID=nbisL1-transcript-1;common_tag="gene1";gene_info="gene1";transcript_id="transcript1"
    chr12   HAVANA  transcript  100 500 .   +   .   ID="bbb";Parent=nbisL1-transcript-1;common_tag="gene1";gene_info="gene1";transcript_id="transcript1"
    chr12   HAVANA  exon    100 500 .   +   .   ID="ccc";Parent="bbb";common_tag="gene1"
    chr12   HAVANA  CDS 100 500 .   +   0   ID="ddd";Parent="bbb";common_tag="gene1"
    chr12   HAVANA  gene    100 600 .   +   .   ID=nbisL1-transcript-2;common_tag="gene1";gene_info="gene1";transcript_id="transcript2"
    chr12   HAVANA  transcript  100 600 .   +   .   ID="bbb2";Parent=nbisL1-transcript-2;common_tag="gene1";gene_info="gene1";transcript_id="transcript2"
    chr12   HAVANA  exon    100 600 .   +   .   ID="ccc2";Parent="bbb2";common_tag="gene1"
    chr12   HAVANA  CDS 100 600 .   +   0   ID="ddd2";Parent="bbb2";common_tag="gene1"
    chr12   HAVANA  gene    1000    5000    .   +   .   ID=nbisL1-transcript-3;common_tag="gene2";gene_info="gene2";transcript_id="transcript3"
    chr12   HAVANA  transcript  1000    5000    .   +   .   ID="yyy";Parent=nbisL1-transcript-3;common_tag="gene2";gene_info="gene2";transcript_id="transcript3"
    chr12   HAVANA  exon    1000    5000    .   +   .   ID="zzz";Parent="yyy";common_tag="gene2"
    chr12   HAVANA  CDS 1000    5000    .   +   0   ID="www";Parent="yyy";common_tag="gene2"

A way to fix that is to use a common attribute. Here you could use `common_tag`, `transcript_id`, `gene_info`.  
Example:

  * `agat_convert_sp_gxf2gxf.pl --gff input.gxf --ct common_tag`  

  This will use the parsing approach 2 (only using common attribute). This will work even if transcript isoform exists.

Correct output:

    chr12   HAVANA  gene    100 600 .   +   .   ID=nbisL1-transcript-1;common_tag="gene1";gene_info="gene1";transcript_id="transcript1"
    chr12   HAVANA  transcript  100 500 .   +   .   ID="bbb";Parent=nbisL1-transcript-1;common_tag="gene1";gene_info="gene1";transcript_id="transcript1"
    chr12   HAVANA  exon    100 500 .   +   .   ID="ccc";Parent="bbb";common_tag="gene1"
    chr12   HAVANA  CDS 100 500 .   +   0   ID="ddd";Parent="bbb";common_tag="gene1"
    chr12   HAVANA  transcript  100 600 .   +   .   ID="bbb2";Parent=nbisL1-transcript-1;common_tag="gene1";gene_info="gene1";transcript_id="transcript2"
    chr12   HAVANA  exon    100 600 .   +   .   ID="ccc2";Parent="bbb2";common_tag="gene1"
    chr12   HAVANA  CDS 100 600 .   +   0   ID="ddd2";Parent="bbb2";common_tag="gene1"
    chr12   HAVANA  gene    1000    5000    .   +   .   ID=nbisL1-transcript-2;common_tag="gene2";gene_info="gene2";transcript_id="transcript3"
    chr12   HAVANA  transcript  1000    5000    .   +   .   ID="yyy";Parent=nbisL1-transcript-2;common_tag="gene2";gene_info="gene2";transcript_id="transcript3"
    chr12   HAVANA  exon    1000    5000    .   +   .   ID="zzz";Parent="yyy";common_tag="gene2"
    chr12   HAVANA  CDS 1000    5000    .   +   0   ID="www";Parent="yyy";common_tag="gene2"

  * `agat_convert_sp_gxf2gxf.pl --gff input.gxf --ct transcript_id`  

  This will use the parsing approach 2 (common attribute transcript_id) for transcrtipt features and approach 3 (sequential) for subfeatures, which do not have the transcript_id attribute.
  /!\ This will not work properly if you have isoforms. Indeed in the case you have isoforms, each transcript will have its own gene feature.

    chr12    HAVANA  gene    100 500 .   +   .   ID="transcript1";common_tag="gene1";gene_info="gene1";transcript_id="transcript1"
    chr12   HAVANA  transcript  100 500 .   +   .   ID="bbb";Parent="transcript1";common_tag="gene1";gene_info="gene1";transcript_id="transcript1"
    chr12   HAVANA  exon    100 500 .   +   .   ID="ccc";Parent="bbb";common_tag="gene1"
    chr12   HAVANA  CDS 100 500 .   +   0   ID="ddd";Parent="bbb";common_tag="gene1"
    chr12   HAVANA  gene    100 600 .   +   .   ID="transcript2";common_tag="gene1";gene_info="gene1";transcript_id="transcript2"
    chr12   HAVANA  transcript  100 600 .   +   .   ID="bbb2";Parent="transcript2";common_tag="gene1";gene_info="gene1";transcript_id="transcript2"
    chr12   HAVANA  exon    100 600 .   +   .   ID="ccc2";Parent="bbb2";common_tag="gene1"
    chr12   HAVANA  CDS 100 600 .   +   0   ID="ddd2";Parent="bbb2";common_tag="gene1"
    chr12   HAVANA  gene    1000    5000    .   +   .   ID="transcript3";common_tag="gene2";gene_info="gene2";transcript_id="transcript3"
    chr12   HAVANA  transcript  1000    5000    .   +   .   ID="yyy";Parent="transcript3";common_tag="gene2";gene_info="gene2";transcript_id="transcript3"
    chr12   HAVANA  exon    1000    5000    .   +   .   ID="zzz";Parent="yyy";common_tag="gene2"
    chr12   HAVANA  CDS 1000    5000    .   +   0   ID="www";Parent="yyy";common_tag="gene2"

  * `agat_convert_sp_gxf2gxf.pl --gff input.gxf --ct gene_info`

  This will use the parsing approach 2 (common attribute gene_info) for transcrtipt features and approach 3 (sequential) for subfeatures, which do not have the transcript_id attribute. 

    chr12   HAVANA  gene    100 600 .   +   .   ID="gene1";common_tag="gene1";gene_info="gene1";transcript_id="transcript1"
    chr12   HAVANA  transcript  100 500 .   +   .   ID="bbb";Parent="gene1";common_tag="gene1";gene_info="gene1";transcript_id="transcript1"
    chr12   HAVANA  exon    100 500 .   +   .   ID="ccc";Parent="bbb";common_tag="gene1"
    chr12   HAVANA  CDS 100 500 .   +   0   ID="ddd";Parent="bbb";common_tag="gene1"
    chr12   HAVANA  transcript  100 600 .   +   .   ID="bbb2";Parent="gene1";common_tag="gene1";gene_info="gene1";transcript_id="transcript2"
    chr12   HAVANA  exon    100 600 .   +   .   ID="ccc2";Parent="bbb2";common_tag="gene1"
    chr12   HAVANA  CDS 100 600 .   +   0   ID="ddd2";Parent="bbb2";common_tag="gene1"
    chr12   HAVANA  gene    1000    5000    .   +   .   ID="gene2";common_tag="gene2";gene_info="gene2";transcript_id="transcript3"
    chr12   HAVANA  transcript  1000    5000    .   +   .   ID="yyy";Parent="gene2";common_tag="gene2";gene_info="gene2";transcript_id="transcript3"
    chr12   HAVANA  exon    1000    5000    .   +   .   ID="zzz";Parent="yyy";common_tag="gene2"
    chr12   HAVANA  CDS 1000    5000    .   +   0   ID="www";Parent="yyy";common_tag="gene2"


#### Level1 and Level2 feature type missing (Only Level3 features!)

In such case the sequential approach cannot be used (Indeed no level1 (e.g. gene) and no lelve2 (e.g. mrna) feature is present in the file)

 * Case with Parent/ID transcript_id/gene_id relationships.

If you have isoforms (for Eukaryote organism) in your files and the `common attribute` used is not set properly you can end up with isoforms having independent parent gene features. See below for more details.

Input:  

    chr12   HAVANA  exon    100 500 .   +   .   ID=exon1;Parent=transcript1
    chr12   HAVANA  CDS 100 500 .   +   0   ID=cds-1;Parent=transcript1
    chr12   HAVANA  exon    100 600 .   +   .   ID=exon2;Parent=transcript2
    chr12   HAVANA  CDS 100 600 .   +   0   ID=cds-2;Parent=transcript2
    chr12   HAVANA  exon    700 900 .   +   .   ID=exonb;Parent=transcriptb
    chr12   HAVANA  CDS 700 900 .   +   0   ID=cds-b;Parent=transcriptb

`agat_convert_sp_gxf2gxf.pl --gff input.gff`

    chr12   HAVANA  gene    100 500 .   +   .   ID=nbis-gene-1
    chr12   HAVANA  mRNA    100 500 .   +   .   ID=transcript1;Parent=nbis-gene-1
    chr12   HAVANA  exon    100 500 .   +   .   ID=exon1;Parent=transcript1
    chr12   HAVANA  CDS 100 500 .   +   0   ID=cds-1;Parent=transcript1
    chr12   HAVANA  gene    100 600 .   +   .   ID=nbis-gene-2
    chr12   HAVANA  mRNA    100 600 .   +   .   ID=transcript2;Parent=nbis-gene-2
    chr12   HAVANA  exon    100 600 .   +   .   ID=exon2;Parent=transcript2
    chr12   HAVANA  CDS 100 600 .   +   0   ID=cds-2;Parent=transcript2
    chr12   HAVANA  gene    700 900 .   +   .   ID=nbis-gene-3
    chr12   HAVANA  mRNA    700 900 .   +   .   ID=transcriptb;Parent=nbis-gene-3
    chr12   HAVANA  exon    700 900 .   +   .   ID=exonb;Parent=transcriptb
    chr12   HAVANA  CDS 700 900 .   +   0   ID=cds-b;Parent=transcriptb

/!\\ This is not correct if transcript1 and transcript2 are isoforms.
You need to use a `common attribute` (here `locus_id`) to group the feature properly:

`agat_convert_sp_gxf2gxf.pl --gff testB.gff --ct locus_id`

    chr12   HAVANA  gene    100 600 .   +   .   ID="gene1";locus_id="gene1"
    chr12   HAVANA  mRNA    100 500 .   +   .   ID=transcript1;Parent="gene1";locus_id="gene1"
    chr12   HAVANA  exon    100 500 .   +   .   ID=exon1;Parent=transcript1;locus_id="gene1"
    chr12   HAVANA  CDS 100 500 .   +   0   ID=cds-1;Parent=transcript1;locus_id="gene1"
    chr12   HAVANA  mRNA    100 600 .   +   .   ID=transcript2;Parent="gene1";locus_id="gene1"
    chr12   HAVANA  exon    100 600 .   +   .   ID=exon2;Parent=transcript2;locus_id="gene1"
    chr12   HAVANA  CDS 100 600 .   +   0   ID=cds-2;Parent=transcript2;locus_id="gene1"
    chr12   HAVANA  gene    700 900 .   +   .   ID="gene2";locus_id="gene2"
    chr12   HAVANA  mRNA    700 900 .   +   .   ID=transcriptb;Parent="gene2";locus_id="gene2"
    chr12   HAVANA  exon    700 900 .   +   .   ID=exonb;Parent=transcriptb;locus_id="gene2"
    chr12   HAVANA  CDS 700 900 .   +   0   ID=cds-b;Parent=transcriptb;locus_id="gene2"


* Case without Parent/ID transcript_id/gene_id relationships. Only `common attribute` approach to parse the file can be used (locus_tag and gene_id by defaut).

/!\\ Here the worse case that can append: only level3 features, no Parent/ID transcript_id/gene_id relationships, and the default `common attributes` are absent.

Input:  

    chr12   HAVANA  exon    100 500 .   +   .   ID=exon1;locus_id="gene1"
    chr12   HAVANA  CDS 100 500 .   +   0   ID=cds-1;locus_id="gene1"
    chr12   HAVANA  exon    510 600 .   +   .   ID=exon2;locus_id="gene1"
    chr12   HAVANA  CDS 510 600 .   +   0   ID=cds-2;locus_id="gene1"
    chr12   HAVANA  exon    700 900 .   +   .   ID=exonb;locus_id="gene2"
    chr12   HAVANA  CDS 700 900 .   +   0   ID=cds-b;locus_id="gene2"

Output:  

    chr12   HAVANA  gene    100 900 .   +   .   ID=nbis-gene-1;locus_id="gene1"
    chr12   HAVANA  mRNA    100 900 .   +   .   ID=nbisL2-exon-1;Parent=nbis-gene-1;locus_id="gene1"
    chr12   HAVANA  exon    100 600 .   +   .   ID=exon1;Parent=nbisL2-exon-1;locus_id="gene1"
    chr12   HAVANA  exon    700 900 .   +   .   ID=exonb;Parent=nbisL2-exon-1;plocus_id="gene2"
    chr12   HAVANA  CDS 100 500 .   +   0   ID=cds-1;Parent=nbisL2-exon-1;locus_id="gene1"
    chr12   HAVANA  CDS 100 600 .   +   0   ID=cds-2;Parent=nbisL2-exon-1;locus_id="gene1"
    chr12   HAVANA  CDS 700 900 .   +   0   ID=cds-b;Parent=nbisL2-exon-1;locus_id="gene2"

/!\\ All features are collected under a single gene and mRNA feature, which is WRONG.

As the default `common attribute` are absent (gene_id or locus_tag), you have to inform AGAT what attribute to use to group features together properly, here `locus_id`.  

`agat_convert_sp_gxf2gxf.pl --gff testC.gff --ct locus_id`
Output:

    chr12   HAVANA  gene    100 600 .   +   .   ID=nbis-gene-1;locus_id="gene1"
    chr12   HAVANA  mRNA    100 600 .   +   .   ID=nbisL2-exon-1;Parent=nbis-gene-1;locus_id="gene1"
    chr12   HAVANA  exon    100 600 .   +   .   ID=exon1;Parent=nbisL2-exon-1;locus_id="gene1"
    chr12   HAVANA  CDS 100 500 .   +   0   ID=cds-1;Parent=nbisL2-exon-1;locus_id="gene1"
    chr12   HAVANA  CDS 100 600 .   +   0   ID=cds-2;Parent=nbisL2-exon-1;locus_id="gene1"
    chr12   HAVANA  gene    700 900 .   +   .   ID=nbis-gene-2;locus_id="gene2"
    chr12   HAVANA  mRNA    700 900 .   +   .   ID=nbisL2-exon-2;Parent=nbis-gene-2;locus_id="gene2"
    chr12   HAVANA  exon    700 900 .   +   .   ID=exonb;Parent=nbisL2-exon-2;locus_id="gene2"
    chr12   HAVANA  CDS 700 900 .   +   0   ID=cds-b;Parent=nbisL2-exon-2;locus_id="gene2"

/!\\ In the case of presence of isoforms (for Eukaryote), it will result of isoforms merged in chimeric transcripts (It will be really unlucky to end up in such situation, because even a human cannot resolve such type of situation. There is no information about isoforms structure...).
In Eukaryote case (multi-exon CDS) and absence of isoform, it will work correctly.

* In the extreme case where you have only one type of feature, you may decide to use the ID as common attribute.

Input:

    chr12   HAVANA  CDS 100 200 .   +   0   ID=cds-1;
    chr12   HAVANA  CDS 510 600 .   +   0   ID=cds-2;
    chr12   HAVANA  CDS 700 900 .   +   0   ID=cds-b;


Output:

    chr12   HAVANA  gene    100 500 .   +   0   ID=nbis-gene-2
    chr12   HAVANA  mRNA    100 500 .   +   0   ID=nbisL2-cds-1;Parent=nbis-gene-2
    chr12   HAVANA  exon    100 500 .   +   .   ID=nbis-exon-1;Parent=nbisL2-cds-1
    chr12   HAVANA  CDS 100 500 .   +   0   ID=cds-1;Parent=nbisL2-cds-1
    chr12   HAVANA  gene    100 600 .   +   0   ID=nbis-gene-3
    chr12   HAVANA  mRNA    100 600 .   +   0   ID=nbisL2-cds-2;Parent=nbis-gene-3
    chr12   HAVANA  exon    100 600 .   +   .   ID=nbis-exon-2;Parent=nbisL2-cds-2
    chr12   HAVANA  CDS 100 600 .   +   0   ID=cds-2;Parent=nbisL2-cds-2
    chr12   HAVANA  gene    700 900 .   +   0   ID=nbis-gene-1
    chr12   HAVANA  mRNA    700 900 .   +   0   ID=nbisL2-cds-3;Parent=nbis-gene-1
    chr12   HAVANA  exon    700 900 .   +   .   ID=nbis-exon-3;Parent=nbisL2-cds-3
    chr12   HAVANA  CDS 700 900 .   +   0   ID=cds-b;Parent=nbisL2-cds-3

/!\\ In the case of Eukaryotes, multi-exon CDS can share or not the same ID (both are allowed by the GFF format). If the CDS chunks share the same ID, the CDS part will be collected properly, if the CDS chuncks do not share the same ID the AGAT will slice it and create a gene/mRNA feature by CDS chunk! If you are in this case (not other attribute to use as `common attribute` you are screwed. This is not suppose to append because such type of input is actually not a GFF or GTF compliant format...)
