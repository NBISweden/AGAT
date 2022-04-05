# How does AGAT work?

All tools taking GFF/GTF as input can be divided in two groups: \_sp\_ and \_sq\_.

* Tools with \_sp\_ prefix

\_sp\_ stands for SLURP. Those tools will charge the file in memory Omniscient data structure. It has a memory cost but makes life smoother. Indeed, it allows to perform complicated tasks in a more time efficient way ( Any features can be accessed at any time by AGAT).
Moreover, it allows to fix all potential errors in the limit of the possibilities given by the format itself.
See the Omniscient section for more information about it.  

* with \_sq\_ prefix

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

  1. Parsing approach 1: by Parent/child relationship  

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

  2. ELSE Parsing approach 2: by a common attribute/tag  

  a common attribute (or common tag) is an attribute value shared by feature that must be grouped together. AGAT uses default attributes (`gene_id` and `locus_tag`) displayed in the log but can be set by the user using the `--ct` parameter).  

Example of relationship made using a common tag (here locus_tag):

    chr12	HAVANA	gene	100	500	.	+	.	locus_tag="gene1"
    chr12	HAVANA	transcript	100	500	.	+	.	locus_tag="gene1";ID="transcript1"
    chr12	HAVANA	exon	100	500	.	+	.	locus_tag="gene1";ID=exon1;
    chr12	HAVANA	CDS	100	500	.	+	0	locus_tag="gene1";ID=cds-1;

  3. ELSE Parsing approach 3: sequentially.

  Reading from top to the botom of the file, level3 features (e.g. exon, CDS, UTR) are attached to the last level2 feature (e.g. mRNA) met, and level2 feature are attached to the last L1 feature (e.g. gene) met. To see the list of features of each level see the corresponding json file (In the share folder in the github repo or using `agat_convert_sp_gxf2gxf.pl --expose`).

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

#### A. Level1 feature type missing and no Parent/gene_id

If you have isoforms (for Eukaryote organism) in your files and the `common attribute` used is not set properly you can end up with isoforms having independent parent gene features. See below for more details.

Here an example of three transcripts from two different genes (isoforms exist - testA.gff):

    chr12   HAVANA  transcript  100 500 .   +   .   ID="bbb";common_tag="gene1";transcript_id="transcript1";gene_info="gene1"
    chr12   HAVANA  exon    100 500 .   +   .   ID="ccc";common_tag="gene1"
    chr12   HAVANA  CDS 100 500 .   +   0   ID="ddd";common_tag="gene1"
    chr12   HAVANA  transcript  100 600 .   +   .   ID="bbb2";common_tag="gene1";transcript_id="transcript2";gene_info="gene1"
    chr12   HAVANA  exon    100 600 .   +   .   ID="ccc2";common_tag="gene1"
    chr12   HAVANA  CDS 100 600 .   +   0   ID="ddd2";common_tag="gene1"
    chr12   HAVANA  transcript  1000    5000    .   +   .   ID="yyy";common_tag="gene2";transcript_id="transcript3";gene_info="gene2"
    chr12   HAVANA  exon    1000    5000    .   +   .   ID="zzz";common_tag="gene2"
    chr12   HAVANA  CDS 1000    5000    .   +   0   ID="www";common_tag="gene2"


* /!\ Be careful in Eukaryote annotation containing isoforms. Indeed AGAT will create a gene feature by transcript. As in the example this is wrong because transcript1 and transcript2 should be attached to the same gene:  

`agat_convert_sp_gxf2gxf.pl --gff testA.gff`

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

* ! A way to fix that is to use a common attribute. Here you could use `common_tag`, `transcript_id`, `gene_info`:

`agat_convert_sp_gxf2gxf.pl --gff testA.gff --ct common_tag`  

This will work well even if transcript isoforms exist. This will use the parsing approach 2 (only using common attribute).

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

`agat_convert_sp_gxf2gxf.pl --gff testA.gff --ct gene_info`

This will work well even if transcript isoforms exist. This will use the parsing approach 2 (common attribute gene_info) for transcript features and approach 3 (sequential) for subfeatures, which do not have the transcript_id attribute.

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

`agat_convert_sp_gxf2gxf.pl --gff testA.gff --ct transcript_id`  

/!\ In our case, using `transcript_id` is not a good choice. Indeed each transcript will have its own gene feature, so isoform will not be linked to the same gene feature as expected. This will use the parsing approach 2 (common attribute transcript_id) for transcript features and approach 3 (sequential) for subfeatures that do not have the transcript_id attribute.

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


#### B. Level1 and Level2 feature types missing (Only Level3 features!)

In such case the sequential approach cannot be used (Indeed no level1 (e.g. gene) and no lelve2 (e.g. mrna) feature is present in the file). So the presence of parent/ID transcript_id/gene_id relationships and/or a proper common attribute is crucial.

1. Case with Parent/ID transcript_id/gene_id relationships.

If you have isoforms (for Eukaryote organism) in your files and the `common attribute` used is not set properly you can end up with isoforms having independent parent gene features. See below for more details.

1.1

Input (testB.gff):  

		chr12	HAVANA	exon	100	500	.	+	.	ID=exon1;Parent=transcript1;locus_id="gene1"
		chr12	HAVANA	CDS	100	500	.	+	0	ID=cds-1;Parent=transcript1;locus_id="gene1"
		chr12	HAVANA	exon	100	600	.	+	.	ID=exon2;Parent=transcript2;locus_id="gene1"
		chr12	HAVANA	CDS	100	600	.	+	0	ID=cds-2;Parent=transcript2;locus_id="gene1"
		chr12	HAVANA	exon	700	900	.	+	.	ID=exonb;Parent=transcriptb;locus_id="gene2"
		chr12	HAVANA	CDS	700	900	.	+	0	ID=cds-b;Parent=transcriptb;locus_id="gene2"

* /!\ Be careful in Eukaryote annotation containing isoforms. Indeed there is no Leve2 feature (e.g. mRNA) to indicate to which parental gene to link isoforms to. By default (see below) the result will be wrong because transcript1 and transcript2 should be attached to the same gene:  

`agat_convert_sp_gxf2gxf.pl --gff testB.gff`

		chr12	HAVANA	gene	100	500	.	+	.	ID=nbis-gene-1;locus_id="gene1"
		chr12	HAVANA	mRNA	100	500	.	+	.	ID=transcript1;Parent=nbis-gene-1;locus_id="gene1"
		chr12	HAVANA	exon	100	500	.	+	.	ID=exon1;Parent=transcript1;locus_id="gene1"
		chr12	HAVANA	CDS	100	500	.	+	0	ID=cds-1;Parent=transcript1;locus_id="gene1"
		chr12	HAVANA	gene	100	600	.	+	.	ID=nbis-gene-2;locus_id="gene1"
		chr12	HAVANA	mRNA	100	600	.	+	.	ID=transcript2;Parent=nbis-gene-2;locus_id="gene1"
		chr12	HAVANA	exon	100	600	.	+	.	ID=exon2;Parent=transcript2;locus_id="gene1"
		chr12	HAVANA	CDS	100	600	.	+	0	ID=cds-2;Parent=transcript2;locus_id="gene1"
		chr12	HAVANA	gene	700	900	.	+	.	ID=nbis-gene-3;locus_id="gene2"
		chr12	HAVANA	mRNA	700	900	.	+	.	ID=transcriptb;Parent=nbis-gene-3;locus_id="gene2"
		chr12	HAVANA	exon	700	900	.	+	.	ID=exonb;Parent=transcriptb;locus_id="gene2"
		chr12	HAVANA	CDS	700	900	.	+	0	ID=cds-b;Parent=transcriptb;locus_id="gene2"

* ! A way to fix that is to use a `common attribute` to group the feature properly:. AGAT uses `locus_tag` and `gene_id` by default.
If you are lucky those attributes already exist. Here they are absent, you can use `locus_id` instead.

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


1.2.

Here we have only level3 features, Parent/ID transcript_id/gene_id relationships present, default `common attributes` ( `locus_tag` or `gene_id`) is set for some features.

Input testF.gff:

		chr12	HAVANA	exon	100	500	.	+	.	ID=exon1;Parent=transcript1;locus_tag="gene1"
		chr12	HAVANA	CDS	100	500	.	+	0	ID=cds-1;Parent=transcript1;locus_tag="gene1"
		chr12	HAVANA	exon	100	600	.	+	.	ID=exon2;Parent=transcript2;locus_tag="gene1"
		chr12	HAVANA	CDS	100	600	.	+	0	ID=cds-2;Parent=transcript2;locus_tag="gene1"
		chr12	HAVANA	exon	700	900	.	+	.	ID=exonb;Parent=transcriptb;locus_tag="gene2"
		chr12	HAVANA	CDS	700	900	.	+	0	ID=cds-b;Parent=transcriptb;locus_tag="gene2"
		chr12	HAVANA	exon	1000	1110	.	+	.	ID=exon4;Parent=transcript4
		chr12	HAVANA	CDS	1000	1110	.	+	0	ID=cds4;Parent=transcript4

`agat_convert_sp_gxf2gxf.pl --gff testF.gff`

		chr12	HAVANA	gene	100	600	.	+	.	ID="gene1";locus_tag="gene1"
		chr12	HAVANA	mRNA	100	500	.	+	.	ID=transcript1;Parent="gene1";locus_tag="gene1"
		chr12	HAVANA	exon	100	500	.	+	.	ID=exon1;Parent=transcript1;locus_tag="gene1"
		chr12	HAVANA	CDS	100	500	.	+	0	ID=cds-1;Parent=transcript1;locus_tag="gene1"
		chr12	HAVANA	mRNA	100	600	.	+	.	ID=transcript2;Parent="gene1";locus_tag="gene1"
		chr12	HAVANA	exon	100	600	.	+	.	ID=exon2;Parent=transcript2;locus_tag="gene1"
		chr12	HAVANA	CDS	100	600	.	+	0	ID=cds-2;Parent=transcript2;locus_tag="gene1"
		chr12	HAVANA	gene	700	900	.	+	.	ID="gene2";locus_tag="gene2"
		chr12	HAVANA	mRNA	700	900	.	+	.	ID=transcriptb;Parent="gene2";locus_tag="gene2"
		chr12	HAVANA	exon	700	900	.	+	.	ID=exonb;Parent=transcriptb;locus_tag="gene2"
		chr12	HAVANA	CDS	700	900	.	+	0	ID=cds-b;Parent=transcriptb;locus_tag="gene2"
		chr12	HAVANA	gene	1000	1110	.	+	.	ID=nbis-gene-1
		chr12	HAVANA	mRNA	1000	1110	.	+	.	ID=transcript4;Parent=nbis-gene-1
		chr12	HAVANA	exon	1000	1110	.	+	.	ID=exon4;Parent=transcript4
		chr12	HAVANA	CDS	1000	1110	.	+	0	ID=cds4;Parent=transcript4

The `common attributes` is used to attach isoforms to a common gene feature. As transcript4 has no common attribute, it will have its own parent features.

2. Case without Parent/ID transcript_id/gene_id relationships. Only `common attribute` approach to parse the file can be used.

2.1.  

Here we have only level3 features, no Parent/ID transcript_id/gene_id relationships, but a default `common attributes` ( `locus_tag` or `gene_id`) is present.

Input testE.gff:

		chr12	HAVANA	exon	100	300	.	+	.	ID=exon1;locus_tag="gene1"
		chr12	HAVANA	CDS	100	300	.	+	0	ID=cds-1;locus_tag="gene1"
		chr12	HAVANA	exon	500	600	.	+	.	ID=exon2;locus_tag="gene1"
		chr12	HAVANA	CDS	500	600	.	+	0	ID=cds-2;locus_tag="gene1"
		chr12	HAVANA	exon	700	900	.	+	.	ID=exonb;locus_tag="gene2"
		chr12	HAVANA	CDS	700	900	.	+	0	ID=cds-b;locus_tag="gene2"

`agat_convert_sp_gxf2gxf.pl --gff testE.gff`

		chr12	HAVANA	gene	100	600	.	+	.	ID=nbis-gene-1;locus_tag="gene1"
		chr12	HAVANA	mRNA	100	600	.	+	.	ID=nbisL2-exon-1;Parent=nbis-gene-1;locus_tag="gene1"
		chr12	HAVANA	exon	100	300	.	+	.	ID=exon1;Parent=nbisL2-exon-1;locus_tag="gene1"
		chr12	HAVANA	exon	500	600	.	+	.	ID=exon2;Parent=nbisL2-exon-1;locus_tag="gene1"
		chr12	HAVANA	CDS	100	300	.	+	0	ID=cds-1;Parent=nbisL2-exon-1;locus_tag="gene1"
		chr12	HAVANA	CDS	500	600	.	+	0	ID=cds-2;Parent=nbisL2-exon-1;locus_tag="gene1"
		chr12	HAVANA	gene	700	900	.	+	.	ID=nbis-gene-2;locus_tag="gene2"
		chr12	HAVANA	mRNA	700	900	.	+	.	ID=nbisL2-exon-2;Parent=nbis-gene-2;locus_tag="gene2"
		chr12	HAVANA	exon	700	900	.	+	.	ID=exonb;Parent=nbisL2-exon-2;locus_tag="gene2"
		chr12	HAVANA	CDS	700	900	.	+	0	ID=cds-b;Parent=nbisL2-exon-2;locus_tag="gene2"

/!\\ In Eukaryote annotation containing isoforms it will not work properly. Indeed, it will result of isoforms merged in chimeric transcripts (It will be really unlucky to end up in such situation, because even a human cannot resolve such type of situation. There is no information about isoforms structure...).
In Eukaryote cases (even for multi-exon CDS) with absence of isoforms, it will work correctly.

2.2

Here the worse case that can append: only level3 features, no Parent/ID transcript_id/gene_id relationships, and the default `common attributes` ( `locus_tag` and `gene_id`) are absent. Sequential approach will be used by AGAT but as there are only level3 features,
all will be linked to only one parent. See below for more details.

Input testC.gff:  

    chr12   HAVANA  exon    100 500 .   +   .   ID=exon1;locus_id="gene1"
    chr12   HAVANA  CDS 100 500 .   +   0   ID=cds-1;locus_id="gene1"
    chr12   HAVANA  exon    510 600 .   +   .   ID=exon2;locus_id="gene1"
    chr12   HAVANA  CDS 510 600 .   +   0   ID=cds-2;locus_id="gene1"
    chr12   HAVANA  exon    700 900 .   +   .   ID=exonb;locus_id="gene2"
    chr12   HAVANA  CDS 700 900 .   +   0   ID=cds-b;locus_id="gene2"

`agat_convert_sp_gxf2gxf.pl --gff testC.gff`  

    chr12   HAVANA  gene    100 900 .   +   .   ID=nbis-gene-1;locus_id="gene1"
    chr12   HAVANA  mRNA    100 900 .   +   .   ID=nbisL2-exon-1;Parent=nbis-gene-1;locus_id="gene1"
    chr12   HAVANA  exon    100 600 .   +   .   ID=exon1;Parent=nbisL2-exon-1;locus_id="gene1"
    chr12   HAVANA  exon    700 900 .   +   .   ID=exonb;Parent=nbisL2-exon-1;plocus_id="gene2"
    chr12   HAVANA  CDS 100 500 .   +   0   ID=cds-1;Parent=nbisL2-exon-1;locus_id="gene1"
    chr12   HAVANA  CDS 100 600 .   +   0   ID=cds-2;Parent=nbisL2-exon-1;locus_id="gene1"
    chr12   HAVANA  CDS 700 900 .   +   0   ID=cds-b;Parent=nbisL2-exon-1;locus_id="gene2"

/!\ All features are collected under a single gene and mRNA feature, which is wrong.

As the default `common attribute` are absent (gene_id or locus_tag), you have to inform AGAT what attribute to use to group features together properly, here `locus_id` is a good choice:  

`agat_convert_sp_gxf2gxf.pl --gff testC.gff --ct locus_id`  

    chr12   HAVANA  gene    100 600 .   +   .   ID=nbis-gene-1;locus_id="gene1"
    chr12   HAVANA  mRNA    100 600 .   +   .   ID=nbisL2-exon-1;Parent=nbis-gene-1;locus_id="gene1"
    chr12   HAVANA  exon    100 600 .   +   .   ID=exon1;Parent=nbisL2-exon-1;locus_id="gene1"
    chr12   HAVANA  CDS 100 500 .   +   0   ID=cds-1;Parent=nbisL2-exon-1;locus_id="gene1"
    chr12   HAVANA  CDS 100 600 .   +   0   ID=cds-2;Parent=nbisL2-exon-1;locus_id="gene1"
    chr12   HAVANA  gene    700 900 .   +   .   ID=nbis-gene-2;locus_id="gene2"
    chr12   HAVANA  mRNA    700 900 .   +   .   ID=nbisL2-exon-2;Parent=nbis-gene-2;locus_id="gene2"
    chr12   HAVANA  exon    700 900 .   +   .   ID=exonb;Parent=nbisL2-exon-2;locus_id="gene2"
    chr12   HAVANA  CDS 700 900 .   +   0   ID=cds-b;Parent=nbisL2-exon-2;locus_id="gene2"

/!\\ In Eukaryote annotation containing isoforms it will not work properly. Indeed, it will result of isoforms merged in chimeric transcripts (It will be really unlucky to end up in such situation, because even a human cannot resolve such type of situation. There is no information about isoforms structure...).
In Eukaryote cases (even for multi-exon CDS) with absence of isoforms, it will work correctly.

3. In the extreme case where you have only one type of feature, you may decide to use the ID as common attribute.

This is the same problem as seen previously. Here the worse case that can append: only level3 features, no Parent/ID transcript_id/gene_id relationships, and the default `common attributes` ( `locus_tag` and `gene_id`) are absent. Sequential approach will be used by AGAT but as there are only level3 features,
all will be linked to only one parent. See below for more details.

Input (testD.gff):

		chr10	Liftoff	CDS	100	300	.	+	0	ID=cds1
		chr10	Liftoff	CDS	600	900	.	+	0	ID=cds2
		chr10	Liftoff	CDS	400	490	.	-	0	ID=cds3

`agat_convert_sp_gxf2gxf.pl --gff testD.gff`  

		chr10	Liftoff	gene	100	900	.	+	.	ID=nbis-gene-1
		chr10	Liftoff	mRNA	100	900	.	+	.	ID=nbisL2-cds-1;Parent=nbis-gene-1
		chr10	Liftoff	exon	100	300	.	+	.	ID=nbis-exon-1;Parent=nbisL2-cds-1
		chr10	Liftoff	exon	400	490	.	+	.	ID=nbis-exon-2;Parent=nbisL2-cds-1
		chr10	Liftoff	exon	600	900	.	+	.	ID=nbis-exon-3;Parent=nbisL2-cds-1
		chr10	Liftoff	CDS	100	300	.	+	0	ID=cds1;Parent=nbisL2-cds-1
		chr10	Liftoff	CDS	400	490	.	-	0	ID=cds3;Parent=nbisL2-cds-1
		chr10	Liftoff	CDS	600	900	.	+	0	ID=cds2;Parent=nbisL2-cds-1

/!\ All features are collected under a single gene and mRNA feature, which is wrong.

`agat_convert_sp_gxf2gxf.pl --gff testD.gff --ct ID`

		chr10	Liftoff	gene	100	300	.	+	0	ID=nbis-gene-1
		chr10	Liftoff	mRNA	100	300	.	+	0	ID=nbisL2-cds-1;Parent=nbis-gene-1
		chr10	Liftoff	exon	100	300	.	+	.	ID=nbis-exon-1;Parent=nbisL2-cds-1
		chr10	Liftoff	CDS	100	300	.	+	0	ID=cds1;Parent=nbisL2-cds-1
		chr10	Liftoff	gene	400	490	.	-	0	ID=nbis-gene-3
		chr10	Liftoff	mRNA	400	490	.	-	0	ID=nbisL2-cds-3;Parent=nbis-gene-3
		chr10	Liftoff	exon	400	490	.	-	.	ID=nbis-exon-3;Parent=nbisL2-cds-3
		chr10	Liftoff	CDS	400	490	.	-	0	ID=cds3;Parent=nbisL2-cds-3
		chr10	Liftoff	gene	600	900	.	+	0	ID=nbis-gene-2
		chr10	Liftoff	mRNA	600	900	.	+	0	ID=nbisL2-cds-2;Parent=nbis-gene-2
		chr10	Liftoff	exon	600	900	.	+	.	ID=nbis-exon-2;Parent=nbisL2-cds-2
		chr10	Liftoff	CDS	600	900	.	+	0	ID=cds2;Parent=nbisL2-cds-2

This case is fine for Prokaryote annotation.  
/!\ For Eukaryote it might work is some conditions:  
A) The annotation should not contain isoforms (Indeed, there is no existing information to decipher to which isoform a CDS will be part of. If isoforms are present, each one will be linked to its own gene feature).  
B) If there are multi-exon CDS, CDS parts must share the same ID (Indeed multi-exon CDS can share or not the same ID. Both way are allowed by the GFF format. If the CDS parts share the same ID, the CDS parts will be collected properly. If the CDS parts do not share the same ID, AGAT will slice it and create a gene/mRNA feature by CDS part!).

4. Case where you have only one type of feature, and some feature have Parent attributes and some other have common attributes.

Input (testG.gff):

		chr12	HAVANA	exon	100	500	.	+	.	ID=exon1;Parent=transcript1
		chr12	HAVANA	CDS	100	500	.	+	0	ID=cds-1;Parent=transcript1
		chr12	HAVANA	exon	100	600	.	+	.	ID=exon2;Parent=transcript2
		chr12	HAVANA	CDS	100	600	.	+	0	ID=cds-2;Parent=transcript2
		chr12	HAVANA	exon	700	900	.	+	.	ID=exonb;locus_tag="gene1"
		chr12	HAVANA	CDS	700	900	.	+	0	ID=cds-b;locus_tag="gene1"
		chr12	HAVANA	exon	1000	1110	.	+	.	ID=exon4;locus_tag="gene2"
		chr12	HAVANA	CDS	1000	1110	.	+	0	ID=cds4;locus_tag="gene2"


`agat_convert_sp_gxf2gxf.pl --gff testG.gff`  

		chr12	HAVANA	gene	100	500	.	+	.	ID=nbis-gene-3
		chr12	HAVANA	mRNA	100	500	.	+	.	ID=transcript1;Parent=nbis-gene-3
		chr12	HAVANA	exon	100	500	.	+	.	ID=exon1;Parent=transcript1
		chr12	HAVANA	CDS	100	500	.	+	0	ID=cds-1;Parent=transcript1
		chr12	HAVANA	gene	100	600	.	+	.	ID=nbis-gene-4
		chr12	HAVANA	mRNA	100	600	.	+	.	ID=transcript2;Parent=nbis-gene-4
		chr12	HAVANA	exon	100	600	.	+	.	ID=exon2;Parent=transcript2
		chr12	HAVANA	CDS	100	600	.	+	0	ID=cds-2;Parent=transcript2
		chr12	HAVANA	gene	700	900	.	+	.	ID=nbis-gene-1;locus_tag="gene1"
		chr12	HAVANA	mRNA	700	900	.	+	.	ID=nbisL2-exon-1;Parent=nbis-gene-1;locus_tag="gene1"
		chr12	HAVANA	exon	700	900	.	+	.	ID=exonb;Parent=nbisL2-exon-1;locus_tag="gene1"
		chr12	HAVANA	CDS	700	900	.	+	0	ID=cds-b;Parent=nbisL2-exon-1;locus_tag="gene1"
		chr12	HAVANA	gene	1000	1110	.	+	.	ID=nbis-gene-2;locus_tag="gene2"
		chr12	HAVANA	mRNA	1000	1110	.	+	.	ID=nbisL2-exon-2;Parent=nbis-gene-2;locus_tag="gene2"
		chr12	HAVANA	exon	1000	1110	.	+	.	ID=exon4;Parent=nbisL2-exon-2;locus_tag="gene2"
		chr12	HAVANA	CDS	1000	1110	.	+	0	ID=cds4;Parent=nbisL2-exon-2;locus_tag="gene2"

/!\ For Eukaryote annotation with isoforms, features would need to have the Parent attribute along with a common attribute to help AGAT to properly reconstruct the parental features (a single gene feature for isoforms).
