# GFF to GTF conversion
## Review of the main conversion tools

It exists many GFF formats and many GTF formats 
(see [here](https://github.com/NBISweden/GAAS/blob/master/annotation/knowledge/gxf.md) for a complete review) and many tools
to perform the conversion. We will try to see in this review the main differences.

# Table of Contents

 * [Test summary](#test-summary)
 * [The GFF file to convert](#the-gff-file-to-convert)
 * [The converters](#the-converters)
   * [AGAT](#agat)
   * [gffread](#gffread) 
   * [GenomeTools](#genometools)
   * [ea-utils](#ea-utils)
   * [TransDecoder](#transdecoder)
   * [Kent utils](#kent-utils)
 * [Feature types in GTF versions](#feature-types-in-gtf-versions)

## Test summary

tool | respect GTF format | UTR conserved | attribute conserved | Stop codon removed from CDS | Comment
-- | -- | -- | -- | -- | -- |
[AGAT](https://github.com/NBISweden/AGAT) | Yes - All (default GTF3) | Yes it converts UTR terms to the appropriate ones according to the GTF version selected.| Yes - All | Yes | Can take any GTF GFF as input. The only one keeping comments at the beginning of the file.
[gffread](https://github.com/gpertea/gffread) | No - They say GTF2.2 but it is not: transcript should be removed; start_codon and stop_codon should stay. | No | No  |  No | 
[GenomeTools](https://github.com/genometools/genometools) | No - only CDS and exon kept | No | No |  No | gene_id and transcript_id get new identifiers.
[ea-utils](https://github.com/ExpressionAnalysis/ea-utils) |  No - only CDS and exon kept | No | No | No | 
[TransDecoder](https://github.com/TransDecoder/TransDecoder) |  No - start and stop codon removed | No | Name only | No | Needs the fasta file for the conversion. Location of the last CDS modified and incorrect
[Kent utils](http://hgdownload.cse.ucsc.edu/admin/exe/) | No - gene is missing or transcript is superfluous to be compliant to one of the GTF format | No | No |  Yes | Create a new attribute 'gene_name'. 

## The GFF file to convert

The test file is a GFF3 file:

```
##gff-version 3
# This is a test sample
scaffold625	maker	gene	337818	343277	.	+	.	ID=CLUHARG00000005458;Name=TUBB3_2
scaffold625	maker	mRNA	337818	343277	.	+	.	ID=CLUHART00000008717;Parent=CLUHARG00000005458
scaffold625	maker	tss	337916	337918	.	+	.	ID=CLUHART00000008717:tss;Parent=CLUHART00000008717
scaffold625	maker	start_codon	337916	337918	.	+	.	ID=CLUHART00000008717:start;Parent=CLUHART00000008717
scaffold625	maker	CDS	337915	337971	.	+	0	ID=CLUHART00000008717:cds;Parent=CLUHART00000008717
scaffold625	maker	CDS	340733	340841	.	+	0	ID=CLUHART00000008717:cds;Parent=CLUHART00000008717
scaffold625	maker	CDS	341518	341628	.	+	2	ID=CLUHART00000008717:cds;Parent=CLUHART00000008717
scaffold625	maker	CDS	341964	343033	.	+	2	ID=CLUHART00000008717:cds;Parent=CLUHART00000008717
scaffold625	maker	stop_codon	343031	343033	.	+	.	ID=CLUHART00000008717:stop;Parent=CLUHART00000008717
scaffold625	maker	exon	337818	337971	.	+	.	ID=CLUHART00000008717:exon1;Parent=CLUHART00000008717
scaffold625	maker	exon	340733	340841	.	+	.	ID=CLUHART00000008717:exon2;Parent=CLUHART00000008717
scaffold625	maker	exon	341518	341628	.	+	.	ID=CLUHART00000008717:exon3;Parent=CLUHART00000008717
scaffold625	maker	exon	341964	343277	.	+	.	ID=CLUHART00000008717:exon4;Parent=CLUHART00000008717
scaffold625	maker	five_prime_utr	337818	337914	.	+	.	ID=CLUHART00000008717:five_prime_utr;Parent=CLUHART00000008717
scaffold625	maker	three_prime_UTR	343034	343277	.	+	.	ID=CLUHART00000008717:three_prime_utr;Parent=CLUHART00000008717
```

## The converters

### AGAT

AGAT v0.5.1  

`agat_convert_sp_gff2gtf.pl --gff 1_test.gff -o 1_test_agat.gtf`

```
##gtf-version 3
##This is a test sample
scaffold625	maker	gene	337818	343277	.	+	.	gene_id "CLUHARG00000005458"; ID "CLUHARG00000005458"; Name "TUBB3_2";
scaffold625	maker	transcript	337818	343277	.	+	.	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717"; ID "CLUHART00000008717"; Parent "CLUHARG00000005458"; original_biotype "mrna";
scaffold625	maker	exon	337818	337971	.	+	.	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717"; ID "CLUHART00000008717:exon1"; Parent "CLUHART00000008717";
scaffold625	maker	exon	340733	340841	.	+	.	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717"; ID "CLUHART00000008717:exon2"; Parent "CLUHART00000008717";
scaffold625	maker	exon	341518	341628	.	+	.	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717"; ID "CLUHART00000008717:exon3"; Parent "CLUHART00000008717";
scaffold625	maker	exon	341964	343277	.	+	.	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717"; ID "CLUHART00000008717:exon4"; Parent "CLUHART00000008717";
scaffold625	maker	CDS	337915	337971	.	+	0	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717"; ID "CLUHART00000008717:cds"; Parent "CLUHART00000008717";
scaffold625	maker	CDS	340733	340841	.	+	0	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717"; ID "CLUHART00000008717:cds"; Parent "CLUHART00000008717";
scaffold625	maker	CDS	341518	341628	.	+	2	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717"; ID "CLUHART00000008717:cds"; Parent "CLUHART00000008717";
scaffold625	maker	CDS	341964	343030	.	+	2	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717"; ID "CLUHART00000008717:cds"; Parent "CLUHART00000008717";
scaffold625	maker	five_prime_utr	337818	337914	.	+	.	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717"; ID "CLUHART00000008717:five_prime_utr"; Parent "CLUHART00000008717";
scaffold625	maker	start_codon	337916	337918	.	+	.	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717"; ID "CLUHART00000008717:start"; Parent "CLUHART00000008717";
scaffold625	maker	stop_codon	343031	343033	.	+	.	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717"; ID "CLUHART00000008717:stop"; Parent "CLUHART00000008717";
scaffold625	maker	three_prime_utr	343034	343277	.	+	.	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717"; ID "CLUHART00000008717:three_prime_utr"; Parent "CLUHART00000008717"; original_biotype "three_prime_UTR";
```

### gffread

gffread 0.11.4

`gffread -E 1_test.gff -T -o  1_test_gffread.gtf` 

```
scaffold625	maker	transcript	337818	343277	.	+	.	transcript_id "CLUHART00000008717"; gene_id "CLUHARG00000005458";
scaffold625	maker	exon	337818	337971	.	+	.	transcript_id "CLUHART00000008717"; gene_id "CLUHARG00000005458";
scaffold625	maker	exon	340733	340841	.	+	.	transcript_id "CLUHART00000008717"; gene_id "CLUHARG00000005458";
scaffold625	maker	exon	341518	341628	.	+	.	transcript_id "CLUHART00000008717"; gene_id "CLUHARG00000005458";
scaffold625	maker	exon	341964	343277	.	+	.	transcript_id "CLUHART00000008717"; gene_id "CLUHARG00000005458";
scaffold625	maker	CDS	337915	337971	.	+	0	transcript_id "CLUHART00000008717"; gene_id "CLUHARG00000005458";
scaffold625	maker	CDS	340733	340841	.	+	0	transcript_id "CLUHART00000008717"; gene_id "CLUHARG00000005458";
scaffold625	maker	CDS	341518	341628	.	+	2	transcript_id "CLUHART00000008717"; gene_id "CLUHARG00000005458";
scaffold625	maker	CDS	341964	343033	.	+	2	transcript_id "CLUHART00000008717"; gene_id "CLUHARG00000005458";
```

### GenomeTools

GenomeTools 1.6.1  
The help says it convert into GTF2.2

`gt gff3_to_gtf 1_test.gff > 1_test_genometools.gtf`

```
scaffold625	maker	exon	337818	337971	.	+	.	gene_id "1"; transcript_id "1.1";
scaffold625	maker	exon	340733	340841	.	+	.	gene_id "1"; transcript_id "1.1";
scaffold625	maker	exon	341518	341628	.	+	.	gene_id "1"; transcript_id "1.1";
scaffold625	maker	exon	341964	343277	.	+	.	gene_id "1"; transcript_id "1.1";
scaffold625	maker	CDS	337915	337971	.	+	0	gene_id "1"; transcript_id "1.1";
scaffold625	maker	CDS	340733	340841	.	+	0	gene_id "1"; transcript_id "1.1";
scaffold625	maker	CDS	341518	341628	.	+	2	gene_id "1"; transcript_id "1.1";
scaffold625	maker	CDS	341964	343033	.	+	2	gene_id "1"; transcript_id "1.1";
```

### ea-utils

[ea-utils](https://github.com/ExpressionAnalysis/ea-utils) commit 2b3d8c5d148801c98a2b3f3d54009a72c5b99521

`./gff2gtf-eautils test_1.gff >  1_test_ea-utils.gtf`

```
scaffold625	maker	exon	337818	337971	0	+	.	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717:CLUHARG00000005458";
scaffold625	maker	CDS	337915	337971	0	+	0	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717:CLUHARG00000005458";
scaffold625	maker	CDS	340733	340841	0	+	0	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717:CLUHARG00000005458";
scaffold625	maker	exon	340733	340841	0	+	.	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717:CLUHARG00000005458";
scaffold625	maker	CDS	341518	341628	0	+	2	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717:CLUHARG00000005458";
scaffold625	maker	exon	341518	341628	0	+	.	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717:CLUHARG00000005458";
scaffold625	maker	CDS	341964	343033	0	+	2	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717:CLUHARG00000005458";
scaffold625	maker	exon	341964	343277	0	+	.	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717:CLUHARG00000005458";
```

### TransDecoder

Transdecoder  v5.5.0

`gff3_gene_to_gtf_format.pl test_1.gff test_1.fa > 1_test_transdecoder.gtf`

```
scaffold625	maker	gene	337818	343277	0	+	.	gene_id "CLUHARG00000005458"; Name "TUBB3_2";
scaffold625	maker	transcript	337818	343277	0	+	.	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717"; Name "TUBB3_2";
scaffold625	maker	exon	337818	337971	0	+	.	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717"; Name "TUBB3_2";
scaffold625	maker	CDS	337818	337971	0	+	.	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717"; Name "TUBB3_2";
scaffold625	maker	exon	340733	340841	0	+	.	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717"; Name "TUBB3_2";
scaffold625	maker	CDS	340733	340841	0	+	.	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717"; Name "TUBB3_2";
scaffold625	maker	exon	341518	341628	0	+	.	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717"; Name "TUBB3_2";
scaffold625	maker	CDS	341518	341628	0	+	.	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717"; Name "TUBB3_2";
scaffold625	maker	exon	341964	343277	0	+	.	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717"; Name "TUBB3_2";
scaffold625	maker	CDS	341964	343277	0	+	.	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717"; Name "TUBB3_2";
```

### Kent utils

version from 26-Feb-2020

`./gff3ToGenePred.dms 1_test.gff temp.genePred`
`./genePredToGtf.dms file temp.genePred 1_test_genePred.gtf`

```
scaffold625	temp.genePred	transcript	337818	343277	.	+	.	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717";  gene_name "CLUHARG00000005458";
scaffold625	temp.genePred	exon	337818	337971	.	+	.	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717"; exon_number "1"; exon_id "CLUHART00000008717.1"; gene_name "CLUHARG00000005458";
scaffold625	temp.genePred	CDS	337915	337971	.	+	0	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717"; exon_number "1"; exon_id "CLUHART00000008717.1"; gene_name "CLUHARG00000005458";
scaffold625	temp.genePred	exon	340733	340841	.	+	.	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717"; exon_number "2"; exon_id "CLUHART00000008717.2"; gene_name "CLUHARG00000005458";
scaffold625	temp.genePred	CDS	340733	340841	.	+	0	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717"; exon_number "2"; exon_id "CLUHART00000008717.2"; gene_name "CLUHARG00000005458";
scaffold625	temp.genePred	exon	341518	341628	.	+	.	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717"; exon_number "3"; exon_id "CLUHART00000008717.3"; gene_name "CLUHARG00000005458";
scaffold625	temp.genePred	CDS	341518	341628	.	+	2	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717"; exon_number "3"; exon_id "CLUHART00000008717.3"; gene_name "CLUHARG00000005458";
scaffold625	temp.genePred	exon	341964	343277	.	+	.	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717"; exon_number "4"; exon_id "CLUHART00000008717.4"; gene_name "CLUHARG00000005458";
scaffold625	temp.genePred	CDS	341964	343030	.	+	2	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717"; exon_number "4"; exon_id "CLUHART00000008717.4"; gene_name "CLUHARG00000005458";
scaffold625	temp.genePred	start_codon	337915	337917	.	+	0	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717"; exon_number "1"; exon_id "CLUHART00000008717.1"; gene_name "CLUHARG00000005458";
scaffold625	temp.genePred	stop_codon	343031	343033	.	+	0	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717"; exon_number "4"; exon_id "CLUHART00000008717.4"; gene_name "CLUHARG00000005458";
```

# Feature types in GTF versions

GTF version | feature type accepted |
-- | -- |
GTF3 | gene,  transcript,  exon,  CDS,  Selenocysteine,  start_codon,  stop_codon,  three_prime_utr,  five_prime_utr 
GTF2_5 | gene,  transcript,  exon,  CDS,  UTR,  start_codon,  stop_codon,  Selenocysteine 
GTF2_2 | CDS,  start_codon,  stop_codon,  5UTR,  3UTR,  inter,  inter_CNS,  intron_CNS,  exon 
GTF2_1 | CDS,  start_codon,  stop_codon,  exon,  5UTR,  3UTR 
GTF2 | CDS,  start_codon,  stop_codon,  exon 
GTF1 | CDS,  start_codon,  stop_codon,  exon,  intron 
