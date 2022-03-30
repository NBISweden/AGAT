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
Below you will find more information about peculiarity of the Omniscient structure,
and the parsing approach used.

### Omniscient data structure

The method create a hash structure containing all the data in memory. We call it OMNISCIENT. The OMNISCIENT structure is a three levels structure:
```
$omniscient{level1}{tag_l1}{level1_id} = feature <= tag could be gene, match  
$omniscient{level2}{tag_l2}{idY} = @featureListL2 <= tag could be mRNA,rRNA,tRNA,etc. idY is a level1_id (know as Parent attribute within the level2 feature). The @featureList is a list to be able to manage isoform cases.  
$omniscient{level3}{tag_l3}{idZ} =  @featureListL3 <= tag could be exon,cds,utr3,utr5,etc. idZ is the ID of a level2 feature (know as Parent attribute within the level3 feature). The @featureList is a list to be able to put all the feature of a same tag together.  
```

### How does the Omniscient parser work

The Omniscient parser phylosophy:
  * 1) Parse by Parent/child relationship  

example of Parent/ID relationship used by the GFF format:

  chr12	HAVANA	gene	100	500	.	+	.	ID=gene1
  chr12	HAVANA	transcript	100	500	.	+	.	ID=transcript1;Parent=gene1
  chr12	HAVANA	exon	100	500	.	+	.	ID=exon1;Parent=transcript1  
  chr12	HAVANA	CDS	100	500	.	+	0	ID=cds-1;Parent=transcript1

example of gene_id/transcript_id relationship used by the GTF format:

  chr12	HAVANA	gene	100	500	.	+	.	gene_id "gene1";
  chr12	HAVANA	transcript	100	500	.	+	.	gene_id "gene1"; transcript_id "transcript1";
  chr12	HAVANA	exon	100	500	.	+	.	gene_id "gene1"; transcript_id "transcript1"; exon_id=exon1;
  chr12	HAVANA	CDS	100	500	.	+	0	gene_id "gene1"; transcript_id "transcript1"; cds_id=cds-1;

  * 2) ELSE Parse by a common tag  (an attribute value shared by feature that must be grouped together. AGAT uses default attributes displayed in the log but can be set by the user using the --ct parameter).  

example of relationship made using a commong tag (here locus_tag):

  chr12	HAVANA	gene	100	500	.	+	.	locus_tag="gene1"
  chr12	HAVANA	transcript	100	500	.	+	.	locus_tag="gene1";ID="transcript1"
  chr12	HAVANA	exon	100	500	.	+	.	locus_tag="gene1";ID=exon1;
  chr12	HAVANA	CDS	100	500	.	+	0	locus_tag="gene1";ID=cds-1;

  * 3) ELSE Parse sequentially (Level3 features as exon and CDS are grouped in a bucket, and each bucket change at each level2 feature, and bucket are join in a common tag at each new L1 feature).  

example of relationship made sequentially:

    chr12	HAVANA	gene	100	500	.	+	.	ID="aaa"
    chr12	HAVANA	transcript	100	500	.	+	.	ID="bbb"
    chr12	HAVANA	exon	100	500	.	+	.	ID="ccc"
    chr12	HAVANA	CDS	100	500	.	+	0	ID="ddd"
    chr12	HAVANA	gene	1000	5000	.	+	.	ID="xxx"
    chr12	HAVANA	transcript	1000	5000	.	+	.	ID="yyy"
    chr12	HAVANA	exon	1000	5000	.	+	.	ID="zzz"
    chr12	HAVANA	CDS	1000	5000	.	+	0	ID="www"

XXX update info

**/!\\** Case with only level3 features (i.e rast or some prokka files, sequential will not work as expected. Indeed all features will be the child of only one newly created Parent. To create a parent per feature or group of features, a common tag must be used to group them correctly. We use `gene_id` and `locus_tag` by default but you can set up the one of your choice)

To resume by priority of way to parse: **Parent/child relationship > locus_tag > sequential.**  
The parser may used only one or a mix of these approaches according of the peculiarity of the gtf/gff file you provide.

### What can the Omniscient parser do for you

* It creates missing parental features. (e.g if a level2 or level3 feature do not have parental feature(s) we create the missing level2 and/or level1 feature(s)).    
* It creates missing mandatory attributes (ID and/or Parent).  
* It fixes identifier to be uniq.  
* It removes duplicated features (same position, same ID, same Parent).  
* It expands level3 features sharing multiple parents (e.g  if one exon has list of multiple parent mRNA in its Parent attribute, one exon per parent with uniq ID will be created.  
* It fixes feature location errors (e.g an mRNA spanning over its gene location, we fix the gene location).  
* It adds UTR if possible (CDS and exon present).  
* It adds exon if possible (CDS has to be present).  
* It groups features together (if related features are spread at different places in the file).  



Below you will find more information about peculiarity of the Omniscient structure,
and the parsing approach used.

If you need to use the `--locus_tag` option you will have to process the file `agat_convert_sp_gxf2gxf.pl` first  before running any other tool.

The only problem that might occur is when you have a file containing only level3 feature (level1=gene, level2=mRNA,level3=CDS,exon,UTR,etc ) and no relationship described. That might be the case you described previously. Lets look into details:

 * Case1 only feature level3 but ID/Parent(GFF) or gene_id/transcript_id(GTF) relationship:  

Input:  

    chr12	HAVANA	exon	100	500	.	+	.	ID=exon1;Parent=transcript1  
    chr12	HAVANA	CDS	100	500	.	+	0	ID=cds-1;Parent=transcript1
    chr12	HAVANA	exon	700	900	.	+	.	ID=exonb;Parent=transcriptb  
    chr12	HAVANA	CDS	700	900	.	+	0	ID=cds-b;Parent=transcriptb

Output:

    chr12	HAVANA	gene	100	500	.	+	.	ID=gene1
    chr12	HAVANA	mRNA	100	500	.	+	.	ID=transcript1;Parent=gene1
    chr12	HAVANA	exon	100	500	.	+	.	ID=exon1;Parent=transcript1
    chr12	HAVANA	CDS	100	500	.	+	0	ID=cds-1;Parent=transcript1
    chr12	HAVANA	gene	700	900	.	+	.	ID=geneb
    chr12	HAVANA	mRNA	700	900	.	+	.	ID=transcriptb;Parent=geneb
    chr12	HAVANA	exon	700	900	.	+	.	ID=exon1;Parent=transcriptb
    chr12	HAVANA	CDS	700	900	.	+	0	ID=cds-1;Parent=transcriptb

=> Everything is fine

* Case2 only feature level3 but relationship made with a comon tag (e.g. locus_tag or gene_id by defaut):  

Input:  

    chr12	HAVANA	exon	100	500	.	+	.	ID=exon1;gene_id=gene1
    chr12	HAVANA	CDS	100	500	.	+	0	ID=cds-1;gene_id=gene1
    chr12	HAVANA	exon	700	900	.	+	.	ID=exonb;gene_id=geneb
    chr12	HAVANA	CDS	700	900	.	+	0	ID=cds-b;gene_id=geneb


Output:  

    chr12	HAVANA	gene	100	500	.	+	.	ID=gene1;gene_id=gene1
    chr12	HAVANA	mRNA	100	500	.	+	0	ID=transcript1;Parent=geneb;gene_id=gene1
    chr12	HAVANA	exon	100	500	.	+	.	ID=exon1;Parent=transcript1;gene_id=gene1
    chr12	HAVANA	CDS	100	500	.	+	0	ID=cds-1;Parent=transcript1;gene_id=gene1
    chr12	HAVANA	gene	700	900	.	+	.	ID=geneb;gene_id=geneb
    chr12	HAVANA	mRNA	700	900	.	+	0	ID=transcriptb;Parent=geneb;gene_id=geneb
    chr12	HAVANA	exon	700	900	.	+	.	ID=exonb;Parent=transcriptb;gene_id=geneb
    chr12	HAVANA	CDS	700	900	.	+	0	ID=cds-b;Parent=transcriptb;gene_id=geneb


=> Everything is fine.
If the comon tag is not gene_id or locus_tag, you have to inform AGAT which attribute to use to group features together

* Case3 only feature level3 but no ID/Parent relationship and no comon locus tag:  

Input:

    chr12	HAVANA	exon	100	500	.	+	.	ID=exon1
    chr12	HAVANA	CDS	100	500	.	+	0	ID=cds-1
    chr12	HAVANA	exon	700	900	.	+	.	ID=exonb
    chr12	HAVANA	CDS	700	900	.	+	0	ID=cds-b

Output:

    chr12	HAVANA	gene	100	900	.	+	.	ID=metaGENE
    chr12	HAVANA	mRNA	100	900	.	+	.	ID=metamRNA;Parent=metaGENE
    chr12	HAVANA	exon	100	500	.	+	.	ID=exon1;Parent=metamRNA
    chr12	HAVANA	CDS	100	500	.	+	0	ID=cds-1;Parent=metamRNA
    chr12	HAVANA	exon	700	900	.	+	.	ID=exonb;Parent=metamRNA
    chr12	HAVANA	CDS	700	900	.	+	0	ID=cds-b;Parent=metamRNA

=> This is WRONG. All features level 3 will be attached to a uniq gene.

So, I hope you didn't end up in this case3 which is very rare (It should not exists because this type of input is actually not a GFF or GTF but I have already seen it).
