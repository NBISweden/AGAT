# Topological sorting of gff features

It might be critical to have a GFF/GTF file properly sorted:

 * Not properly sorted, a genome browser can bug or give wrong displays
 * Some tools require files sorted in a particular way (e.g.tabix tool from htslib need a GFF sorted by chromosomes and positions).
 * It makes it easy to ready for the human eye

Zhigang Lu has made a nice post about his experience trying to find a way to get a correct topological sorting. See [here](https://zhiganglu.com/post/sort-gff-topologically/).

### Table of Contents

 * [Tests summary](#test-summary)
 * [Example 1](#example-1)
   * [The GFF file to sort](#the-gff-file-to-sort)
   * [Results](results)
     * [AGAT](#agat)
     * [GenomeTools](#genometools)
     * [GFF3sort](#gff3sort)
     * [gffread](#gffread)
 * [Example 2](#example-2)
   * [The GFF file to sort](#the-gff-file-to-sort-2)
   * [Results](results-2)
     * [AGAT](#agat-2)
     * [GenomeTools](#genometools-2)
     * [GFF3sort](#gff3sort-2)
     * [gffread](#gffread-2)


#### Tests summary

tool | option in command line | Type of sorting | Comment
-- | -- | -- | -- |
[AGAT](https://github.com/NBISweden/AGAT) | / | by chromosomes, by gene position, by type (mRNAs then exon, then CDS then alphabetical feature types; then mRNA2 then exon2, then CDS2 then alphabetical feature2 types) | Fix GFF/GTF if needed
[GenomeTools](https://github.com/genometools/genometools) | -sortlines -tidy -retainids | by chromosomes and positions then random feature type | Lines with the same chromosomes and start positions would be placed randomly, so parent feature lines might sometimes be placed after their children lines.
[GenomeTools](https://github.com/genometools/genometools) | -retainids | by chromosomes, by gene position, by type (mRNA then children; then mRNA2 then children2), by position (children are sorted by positions) |
[GFF3sort](https://github.com/billzt/gff3sort) | --precise | by chromosomes and positions then attribute with Parent attribute first.  | move lines with "Parent=" attributes (case insensitive) behind lines without "Parent=" attributes. The goal of GFF3sort is not to obtain a topological sorting but rather getting something that could be indexed optimally by third part tools.
[gffread](https://github.com/gpertea/gffread) | | By default, chromosomes are kept in the order they were found. With --sort-alpha parameter the chromosomes (reference sequences) are sorted alphabetically | /!\ Some feature types are lost e.g. `gene`, `three_prime_UTR`, `five_prime_UTR`, etc...

### Example 1

This test is based on the file used by [Zhigang Lu](https://zhiganglu.com/post/sort-gff-topologically/)

#### The GFF file to sort

```
##gff-version 3
SM_V7_1	AUGUSTUS	exon	103403	103770	.	-	.	Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	exon	103403	103770	.	-	.	Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	gene	103403	151162	0.12	-	.	ID=Smp_315690
SM_V7_1	AUGUSTUS	mRNA	103403	151162	0.02	-	.	ID=Smp_315690.1;Parent=Smp_315690
SM_V7_1	AUGUSTUS	mRNA	103403	151162	0.1	-	.	ID=Smp_315690.2;Parent=Smp_315690
SM_V7_1	AUGUSTUS	three_prime_UTR	103403	103440	.	-	.	Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	three_prime_UTR	103403	103440	.	-	.	Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	CDS	103441	103770	0.93	-	0	ID=Smp_315690.1.cds;Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	CDS	103441	103770	0.96	-	0	ID=Smp_315690.2.cds;Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	CDS	105920	106144	1	-	0	ID=Smp_315690.2.cds;Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	exon	105920	106144	.	-	.	Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	CDS	106876	107159	0.93	-	2	ID=Smp_315690.2.cds;Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	exon	106876	107159	.	-	.	Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	CDS	140582	140849	0.85	-	0	ID=Smp_315690.2.cds;Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	exon	140582	140849	.	-	.	Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	CDS	142981	143205	1	-	0	ID=Smp_315690.1.cds;Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	CDS	142981	143205	1	-	0	ID=Smp_315690.2.cds;Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	exon	142981	143205	.	-	.	Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	exon	142981	143205	.	-	.	Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	CDS	145395	145678	1	-	2	ID=Smp_315690.1.cds;Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	CDS	145395	145678	1	-	2	ID=Smp_315690.2.cds;Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	exon	145395	145678	.	-	.	Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	exon	145395	145678	.	-	.	Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	CDS	151075	151132	1	-	0	ID=Smp_315690.1.cds;Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	CDS	151075	151132	1	-	0	ID=Smp_315690.2.cds;Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	exon	151075	151162	.	-	.	Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	exon	151075	151162	.	-	.	Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	five_prime_UTR	151133	151162	.	-	.	Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	five_prime_UTR	151133	151162	.	-	.	Parent=Smp_315690.2
```

#### Results

##### AGAT

AGAT v0.4.0  

`agat_convert_sp_gxf2gxf.pl --gff test.gff`

```
##gff-version 3
SM_V7_1	AUGUSTUS	gene	103403	151162	0.12	-	.	ID=Smp_315690
SM_V7_1	AUGUSTUS	mRNA	103403	151162	0.02	-	.	ID=Smp_315690.1;Parent=Smp_315690
SM_V7_1	AUGUSTUS	exon	103403	103770	.	-	.	ID=exon-1;Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	exon	142981	143205	.	-	.	ID=exon-6;Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	exon	145395	145678	.	-	.	ID=exon-8;Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	exon	151075	151162	.	-	.	ID=exon-10;Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	CDS	103441	103770	0.93	-	0	ID=Smp_315690.1.cds;Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	CDS	142981	143205	1	-	0	ID=Smp_315690.1.cds;Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	CDS	145395	145678	1	-	2	ID=Smp_315690.1.cds;Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	CDS	151075	151132	1	-	0	ID=Smp_315690.1.cds;Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	five_prime_UTR	151133	151162	.	-	.	ID=five_prime_utr-1;Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	three_prime_UTR	103403	103440	.	-	.	ID=three_prime_utr-1;Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	mRNA	103403	151162	0.1	-	.	ID=Smp_315690.2;Parent=Smp_315690
SM_V7_1	AUGUSTUS	exon	103403	103770	.	-	.	ID=exon-2;Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	exon	105920	106144	.	-	.	ID=exon-3;Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	exon	106876	107159	.	-	.	ID=exon-4;Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	exon	140582	140849	.	-	.	ID=exon-5;Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	exon	142981	143205	.	-	.	ID=exon-7;Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	exon	145395	145678	.	-	.	ID=exon-9;Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	exon	151075	151162	.	-	.	ID=exon-11;Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	CDS	103441	103770	0.96	-	0	ID=Smp_315690.2.cds;Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	CDS	105920	106144	1	-	0	ID=Smp_315690.2.cds;Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	CDS	106876	107159	0.93	-	2	ID=Smp_315690.2.cds;Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	CDS	140582	140849	0.85	-	0	ID=Smp_315690.2.cds;Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	CDS	142981	143205	1	-	0	ID=Smp_315690.2.cds;Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	CDS	145395	145678	1	-	2	ID=Smp_315690.2.cds;Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	CDS	151075	151132	1	-	0	ID=Smp_315690.2.cds;Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	five_prime_UTR	151133	151162	.	-	.	ID=five_prime_utr-2;Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	three_prime_UTR	103403	103440	.	-	.	ID=three_prime_utr-2;Parent=Smp_315690.2
```

##### GenomeTools

GenomeTools 1.6.1

`gt gff3 -sortlines -tidy -retainids test.gff `

```
##gff-version 3
##sequence-region   SM_V7_1 103403 151162
SM_V7_1	AUGUSTUS	three_prime_UTR	103403	103440	.	-	.	Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	three_prime_UTR	103403	103440	.	-	.	Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	exon	103403	103770	.	-	.	Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	exon	103403	103770	.	-	.	Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	gene	103403	151162	0.12	-	.	ID=Smp_315690
SM_V7_1	AUGUSTUS	mRNA	103403	151162	0.02	-	.	ID=Smp_315690.1;Parent=Smp_315690
SM_V7_1	AUGUSTUS	mRNA	103403	151162	0.1	-	.	ID=Smp_315690.2;Parent=Smp_315690
SM_V7_1	AUGUSTUS	CDS	103441	103770	0.96	-	0	ID=Smp_315690.2.cds;Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	CDS	103441	103770	0.93	-	0	ID=Smp_315690.1.cds;Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	CDS	105920	106144	1	-	0	ID=Smp_315690.2.cds;Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	exon	105920	106144	.	-	.	Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	CDS	106876	107159	0.93	-	2	ID=Smp_315690.2.cds;Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	exon	106876	107159	.	-	.	Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	CDS	140582	140849	0.85	-	0	ID=Smp_315690.2.cds;Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	exon	140582	140849	.	-	.	Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	CDS	142981	143205	1	-	0	ID=Smp_315690.1.cds;Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	exon	142981	143205	.	-	.	Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	CDS	142981	143205	1	-	0	ID=Smp_315690.2.cds;Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	exon	142981	143205	.	-	.	Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	CDS	145395	145678	1	-	2	ID=Smp_315690.1.cds;Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	CDS	145395	145678	1	-	2	ID=Smp_315690.2.cds;Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	exon	145395	145678	.	-	.	Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	exon	145395	145678	.	-	.	Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	CDS	151075	151132	1	-	0	ID=Smp_315690.2.cds;Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	CDS	151075	151132	1	-	0	ID=Smp_315690.1.cds;Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	exon	151075	151162	.	-	.	Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	exon	151075	151162	.	-	.	Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	five_prime_UTR	151133	151162	.	-	.	Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	five_prime_UTR	151133	151162	.	-	.	Parent=Smp_315690.2
```

`gt gff3 -retainids test.gff `

```
##gff-version 3
##sequence-region   SM_V7_1 103403 151162
SM_V7_1	AUGUSTUS	gene	103403	151162	0.12	-	.	ID=Smp_315690
SM_V7_1	AUGUSTUS	mRNA	103403	151162	0.02	-	.	ID=Smp_315690.1;Parent=Smp_315690
SM_V7_1	AUGUSTUS	three_prime_UTR	103403	103440	.	-	.	Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	exon	103403	103770	.	-	.	Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	CDS	103441	103770	0.93	-	0	ID=Smp_315690.1.cds;Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	CDS	142981	143205	1	-	0	ID=Smp_315690.1.cds;Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	exon	142981	143205	.	-	.	Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	CDS	145395	145678	1	-	2	ID=Smp_315690.1.cds;Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	exon	145395	145678	.	-	.	Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	CDS	151075	151132	1	-	0	ID=Smp_315690.1.cds;Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	exon	151075	151162	.	-	.	Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	five_prime_UTR	151133	151162	.	-	.	Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	mRNA	103403	151162	0.1	-	.	ID=Smp_315690.2;Parent=Smp_315690
SM_V7_1	AUGUSTUS	three_prime_UTR	103403	103440	.	-	.	Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	exon	103403	103770	.	-	.	Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	CDS	103441	103770	0.96	-	0	ID=Smp_315690.2.cds;Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	CDS	105920	106144	1	-	0	ID=Smp_315690.2.cds;Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	exon	105920	106144	.	-	.	Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	CDS	106876	107159	0.93	-	2	ID=Smp_315690.2.cds;Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	exon	106876	107159	.	-	.	Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	CDS	140582	140849	0.85	-	0	ID=Smp_315690.2.cds;Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	exon	140582	140849	.	-	.	Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	CDS	142981	143205	1	-	0	ID=Smp_315690.2.cds;Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	exon	142981	143205	.	-	.	Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	CDS	145395	145678	1	-	2	ID=Smp_315690.2.cds;Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	exon	145395	145678	.	-	.	Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	CDS	151075	151132	1	-	0	ID=Smp_315690.2.cds;Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	exon	151075	151162	.	-	.	Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	five_prime_UTR	151133	151162	.	-	.	Parent=Smp_315690.2
###
```

##### GFF3sort

GFF3sort 0.1.a1a2bc9

`gff3sort.pl --precise test.gff`

```
##gff-version 3
SM_V7_1	AUGUSTUS	gene	103403	151162	0.12	-	.	ID=Smp_315690
SM_V7_1	AUGUSTUS	mRNA	103403	151162	0.02	-	.	ID=Smp_315690.1;Parent=Smp_315690
SM_V7_1	AUGUSTUS	mRNA	103403	151162	0.1	-	.	ID=Smp_315690.2;Parent=Smp_315690
SM_V7_1	AUGUSTUS	three_prime_UTR	103403	103440	.	-	.	Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	exon	103403	103770	.	-	.	Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	three_prime_UTR	103403	103440	.	-	.	Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	exon	103403	103770	.	-	.	Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	CDS	103441	103770	0.93	-	0	ID=Smp_315690.1.cds;Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	CDS	103441	103770	0.96	-	0	ID=Smp_315690.2.cds;Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	CDS	105920	106144	1	-	0	ID=Smp_315690.2.cds;Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	exon	105920	106144	.	-	.	Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	CDS	106876	107159	0.93	-	2	ID=Smp_315690.2.cds;Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	exon	106876	107159	.	-	.	Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	CDS	140582	140849	0.85	-	0	ID=Smp_315690.2.cds;Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	exon	140582	140849	.	-	.	Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	exon	142981	143205	.	-	.	Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	exon	142981	143205	.	-	.	Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	CDS	142981	143205	1	-	0	ID=Smp_315690.2.cds;Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	CDS	142981	143205	1	-	0	ID=Smp_315690.1.cds;Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	CDS	145395	145678	1	-	2	ID=Smp_315690.1.cds;Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	exon	145395	145678	.	-	.	Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	CDS	145395	145678	1	-	2	ID=Smp_315690.2.cds;Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	exon	145395	145678	.	-	.	Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	CDS	151075	151132	1	-	0	ID=Smp_315690.1.cds;Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	exon	151075	151162	.	-	.	Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	exon	151075	151162	.	-	.	Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	CDS	151075	151132	1	-	0	ID=Smp_315690.2.cds;Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	five_prime_UTR	151133	151162	.	-	.	Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	five_prime_UTR	151133	151162	.	-	.	Parent=Smp_315690.2
```

##### gffread

gffread v0.11.4

`gffread test.gff`

```
# gffread test.gff
# gffread v0.11.4
##gff-version 3
SM_V7_1	AUGUSTUS	mRNA	103403	151162	.	-	.	ID=Smp_315690.1;geneID=Smp_315690
SM_V7_1	AUGUSTUS	exon	103403	103770	.	-	.	Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	exon	142981	143205	.	-	.	Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	exon	145395	145678	.	-	.	Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	exon	151075	151162	.	-	.	Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	CDS	103441	103770	.	-	0	Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	CDS	142981	143205	.	-	0	Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	CDS	145395	145678	.	-	2	Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	CDS	151075	151132	.	-	0	Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	mRNA	103403	151162	.	-	.	ID=Smp_315690.2;geneID=Smp_315690
SM_V7_1	AUGUSTUS	exon	103403	103770	.	-	.	Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	exon	105920	106144	.	-	.	Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	exon	106876	107159	.	-	.	Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	exon	140582	140849	.	-	.	Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	exon	142981	143205	.	-	.	Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	exon	145395	145678	.	-	.	Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	exon	151075	151162	.	-	.	Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	CDS	103441	103770	.	-	0	Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	CDS	105920	106144	.	-	0	Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	CDS	106876	107159	.	-	2	Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	CDS	140582	140849	.	-	0	Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	CDS	142981	143205	.	-	0	Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	CDS	145395	145678	.	-	2	Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	CDS	151075	151132	.	-	0	Parent=Smp_315690.2
```

### Example 2

This test is based on the file used by [GFF3sort](https://github.com/billzt/gff3sort)

#### The GFF file to sort

```
##gff-version 3
###
A01	Cufflinks	mRNA	473	6154	.	-	.	ID=XLOC_001154.41;description=Novel: Intergenic transcript
A01	Cufflinks	exon	473	814	.	-	.	Parent=XLOC_001154.41
A01	Cufflinks	exon	1626	2574	.	-	.	Parent=XLOC_001154.41
A01	Cufflinks	exon	2695	2721	.	-	.	Parent=XLOC_001154.41
A01	Cufflinks	exon	3637	3726	.	-	.	Parent=XLOC_001154.41
A01	Cufflinks	exon	5329	5408	.	-	.	Parent=XLOC_001154.41
A01	Cufflinks	exon	5994	6154	.	-	.	Parent=XLOC_001154.41
###
A01	Cufflinks	mRNA	473	6386	.	-	.	ID=XLOC_001154.42;description=Novel: Intergenic transcript
A01	Cufflinks	exon	473	2024	.	-	.	Parent=XLOC_001154.42
A01	Cufflinks	exon	2615	2721	.	-	.	Parent=XLOC_001154.42
A01	Cufflinks	exon	3637	3726	.	-	.	Parent=XLOC_001154.42
A01	Cufflinks	exon	5329	6386	.	-	.	Parent=XLOC_001154.42
```

#### Results

##### AGAT

AGAT v0.4.0  

`agat_convert_sp_gxf2gxf.pl --gff test2.gff --merge_loci`

```
##gff-version 3
###
A01	Cufflinks	gene	473	6386	.	-	.	ID=nbisL1-mrna-1;description=Novel: Intergenic transcript
A01	Cufflinks	mRNA	473	6154	.	-	.	ID=XLOC_001154.41;Parent=nbisL1-mrna-1;description=Novel: Intergenic transcript
A01	Cufflinks	exon	473	814	.	-	.	ID=exon-1;Parent=XLOC_001154.41
A01	Cufflinks	exon	1626	2574	.	-	.	ID=exon-2;Parent=XLOC_001154.41
A01	Cufflinks	exon	2695	2721	.	-	.	ID=exon-3;Parent=XLOC_001154.41
A01	Cufflinks	exon	3637	3726	.	-	.	ID=exon-4;Parent=XLOC_001154.41
A01	Cufflinks	exon	5329	5408	.	-	.	ID=exon-5;Parent=XLOC_001154.41
A01	Cufflinks	exon	5994	6154	.	-	.	ID=exon-6;Parent=XLOC_001154.41
A01	Cufflinks	mRNA	473	6386	.	-	.	ID=XLOC_001154.42;Parent=nbisL1-mrna-1;description=Novel: Intergenic transcript
A01	Cufflinks	exon	473	2024	.	-	.	ID=exon-7;Parent=XLOC_001154.42
A01	Cufflinks	exon	2615	2721	.	-	.	ID=exon-8;Parent=XLOC_001154.42
A01	Cufflinks	exon	3637	3726	.	-	.	ID=exon-9;Parent=XLOC_001154.42
A01	Cufflinks	exon	5329	6386	.	-	.	ID=exon-10;Parent=XLOC_001154.42
```

##### GenomeTools

GenomeTools 1.6.1

`gt gff3 -sortlines -tidy -retainids test2.gff `

```
##gff-version 3
##sequence-region   A01 473 6386
A01	Cufflinks	exon	473	814	.	-	.	Parent=XLOC_001154.41
A01	Cufflinks	exon	473	2024	.	-	.	Parent=XLOC_001154.42
A01	Cufflinks	mRNA	473	6154	.	-	.	ID=XLOC_001154.41;description=Novel: Intergenic transcript
A01	Cufflinks	mRNA	473	6386	.	-	.	ID=XLOC_001154.42;description=Novel: Intergenic transcript
A01	Cufflinks	exon	1626	2574	.	-	.	Parent=XLOC_001154.41
A01	Cufflinks	exon	2615	2721	.	-	.	Parent=XLOC_001154.42
A01	Cufflinks	exon	2695	2721	.	-	.	Parent=XLOC_001154.41
A01	Cufflinks	exon	3637	3726	.	-	.	Parent=XLOC_001154.42
A01	Cufflinks	exon	3637	3726	.	-	.	Parent=XLOC_001154.41
A01	Cufflinks	exon	5329	5408	.	-	.	Parent=XLOC_001154.41
A01	Cufflinks	exon	5329	6386	.	-	.	Parent=XLOC_001154.42
A01	Cufflinks	exon	5994	6154	.	-	.	Parent=XLOC_001154.41
###
```

`gt gff3 -retainids test2.gff `

```
##gff-version 3
##sequence-region   SM_V7_1 103403 151162
SM_V7_1	AUGUSTUS	gene	103403	151162	0.12	-	.	ID=Smp_315690
SM_V7_1	AUGUSTUS	mRNA	103403	151162	0.02	-	.	ID=Smp_315690.1;Parent=Smp_315690
SM_V7_1	AUGUSTUS	three_prime_UTR	103403	103440	.	-	.	Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	exon	103403	103770	.	-	.	Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	CDS	103441	103770	0.93	-	0	ID=Smp_315690.1.cds;Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	CDS	142981	143205	1	-	0	ID=Smp_315690.1.cds;Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	exon	142981	143205	.	-	.	Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	CDS	145395	145678	1	-	2	ID=Smp_315690.1.cds;Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	exon	145395	145678	.	-	.	Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	CDS	151075	151132	1	-	0	ID=Smp_315690.1.cds;Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	exon	151075	151162	.	-	.	Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	five_prime_UTR	151133	151162	.	-	.	Parent=Smp_315690.1
SM_V7_1	AUGUSTUS	mRNA	103403	151162	0.1	-	.	ID=Smp_315690.2;Parent=Smp_315690
SM_V7_1	AUGUSTUS	three_prime_UTR	103403	103440	.	-	.	Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	exon	103403	103770	.	-	.	Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	CDS	103441	103770	0.96	-	0	ID=Smp_315690.2.cds;Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	CDS	105920	106144	1	-	0	ID=Smp_315690.2.cds;Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	exon	105920	106144	.	-	.	Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	CDS	106876	107159	0.93	-	2	ID=Smp_315690.2.cds;Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	exon	106876	107159	.	-	.	Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	CDS	140582	140849	0.85	-	0	ID=Smp_315690.2.cds;Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	exon	140582	140849	.	-	.	Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	CDS	142981	143205	1	-	0	ID=Smp_315690.2.cds;Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	exon	142981	143205	.	-	.	Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	CDS	145395	145678	1	-	2	ID=Smp_315690.2.cds;Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	exon	145395	145678	.	-	.	Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	CDS	151075	151132	1	-	0	ID=Smp_315690.2.cds;Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	exon	151075	151162	.	-	.	Parent=Smp_315690.2
SM_V7_1	AUGUSTUS	five_prime_UTR	151133	151162	.	-	.	Parent=Smp_315690.2
###
```

##### GFF3sort

GFF3sort 0.1.a1a2bc9

`gff3sort.pl --precise test2.gff`

```
##gff-version 3
A01	Cufflinks	mRNA	473	6154	.	-	.	ID=XLOC_001154.41;description=Novel: Intergenic transcript
A01	Cufflinks	mRNA	473	6386	.	-	.	ID=XLOC_001154.42;description=Novel: Intergenic transcript
A01	Cufflinks	exon	473	2024	.	-	.	Parent=XLOC_001154.42
A01	Cufflinks	exon	473	814	.	-	.	Parent=XLOC_001154.41
A01	Cufflinks	exon	1626	2574	.	-	.	Parent=XLOC_001154.41
A01	Cufflinks	exon	2615	2721	.	-	.	Parent=XLOC_001154.42
A01	Cufflinks	exon	2695	2721	.	-	.	Parent=XLOC_001154.41
A01	Cufflinks	exon	3637	3726	.	-	.	Parent=XLOC_001154.41
A01	Cufflinks	exon	3637	3726	.	-	.	Parent=XLOC_001154.42
A01	Cufflinks	exon	5329	6386	.	-	.	Parent=XLOC_001154.42
A01	Cufflinks	exon	5329	5408	.	-	.	Parent=XLOC_001154.41
A01	Cufflinks	exon	5994	6154	.	-	.	Parent=XLOC_001154.41
```

##### gffread

gffread v0.11.4

`gffread test2.gff`

```
# gffread test2.gff
# gffread v0.11.4
##gff-version 3
A01	Cufflinks	mRNA	473	6154	.	-	.	ID=XLOC_001154.41
A01	Cufflinks	exon	473	814	.	-	.	Parent=XLOC_001154.41
A01	Cufflinks	exon	1626	2574	.	-	.	Parent=XLOC_001154.41
A01	Cufflinks	exon	2695	2721	.	-	.	Parent=XLOC_001154.41
A01	Cufflinks	exon	3637	3726	.	-	.	Parent=XLOC_001154.41
A01	Cufflinks	exon	5329	5408	.	-	.	Parent=XLOC_001154.41
A01	Cufflinks	exon	5994	6154	.	-	.	Parent=XLOC_001154.41
A01	Cufflinks	mRNA	473	6386	.	-	.	ID=XLOC_001154.42
A01	Cufflinks	exon	473	2024	.	-	.	Parent=XLOC_001154.42
A01	Cufflinks	exon	2615	2721	.	-	.	Parent=XLOC_001154.42
A01	Cufflinks	exon	3637	3726	.	-	.	Parent=XLOC_001154.42
A01	Cufflinks	exon	5329	6386	.	-	.	Parent=XLOC_001154.42
```
