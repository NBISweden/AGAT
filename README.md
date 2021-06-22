[![Build Status](https://travis-ci.org/NBISweden/AGAT.svg?branch=master)](https://travis-ci.org/NBISweden/AGAT)
[![Coverage Status](https://coveralls.io/repos/github/NBISweden/AGAT/badge.svg)](https://coveralls.io/github/NBISweden/AGAT)
[![Documentation Status](https://readthedocs.org/projects/agat/badge/?version=latest)](https://agat.readthedocs.io/en/latest/?badge=latest)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/agat/README.html)
[![DOI](https://zenodo.org/badge/222659741.svg)](https://zenodo.org/badge/latestdoi/222659741)
[<img alt="docker_agat" src="https://quay.io/repository/biocontainers/agat/status">](https://quay.io/repository/biocontainers/agat)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/agat/badges/license.svg)](https://anaconda.org/bioconda/agat)
[![Anaconda-Server Badge](https://img.shields.io/conda/dn/bioconda/agat.svg?style=flat)](https://anaconda.org/bioconda/agat)  

AGAT
=========================================
<h2><em>A</em>nother <em>G</em>tf/Gff <em>A</em>nalysis <i>T</i>oolkit</h2>  
Suite of tools to handle gene annotations in any GTF/GFF format.

[<img align="right" src="NBIS.png" width="200" height="100" />](https://nbis.se)

---------------------------

## Table of Contents

   * [What can AGAT do for you?](#what-can-agat-do-for-you)
   * [Installation](#installation)  
       * [Using bioconda](using-bioconda)
          * [Install AGAT](#install-agat)
          * [Update AGAT](#update-agat)
          * [Uninstall AGAT](#uninstall-agat)
       * [Old school - Manually](#old-school---manually)
          * [Install prerequisites](#install-prerequisites)
          * [Install AGAT](#install-agat-1)
          * [Update AGAT](#update-agat-1)
          * [Change to a specific version](#change-to-a-specific-version)
          * [Uninstall AGAT](#uninstall-agat-1)
   * [Usage](#usage)
   * [List of tools](#list-of-tools)
   * [More about the tools](#more-about-the-tools)
   * [Omniscient - Standardisation for a full GFF3 compliant to any tool](#omniscient---standardisation-for-a-full-gff3-compliant-to-any-tool)
      * [Omniscient data structure](#omniscient-data-structure)
      * [How does the Omniscient parser work](#how-does-the-omniscient-parser-work)
      * [What can the Omniscient parser do for you](#what-can-the-omniscient-parser-do-for-you)
      * [examples](#examples)
   * [How to cite?](#how-to-cite)
   * [Publication using AGAT](#publication-using-agat)
   * [Troubleshooting](#troubleshooting)
---------------------------

## What can AGAT do for you?  

It has the power to check, fix, pad missing information (features/attributes) of any kind of GTF and GFF to create complete, sorted and standardised gff3 format.  
The GTF/GFF formats are 9-column text formats used to describe and represent genomic features.
The formats have quite evolved since 1997, and despite well-defined specifications existing nowadays they have a great flexibility allowing holding wide variety of information.
This flexibility has a drawback aspect, there is an incredible amount of flavour of the formats, that can result in problems when using downstream programs.  
For a complete overview of the GTF/GFF formats have a look [here](https://agat.readthedocs.io/en/latest/gxf.html).

Some examples **what AGAT can do**:  
  * standardise any GTF/GFF file into a comprehensive GFF3 format (script with `agat_sp` prefix):  
    * add missing parent features (e.g. gene and mRNA if only CDS/exon exist).  
    * add missing features (e.g. exon and UTR).  
    * add missing mandatory attributes (i.e. ID, Parent).  
    * fix identifier to be uniq.  
    * fix feature location.  
    * remove duplicated features.  
    * group related features (if spread in different places in the file).  
    * sort features.  
    * merge overlapping loci into one single locus (only if option activated).  

  * perform different tasks (using different AGAT's tools):

| task | tool |
| --- | --- |
| **check, fix, pad** missing information into sorted and standardised gff3 | `agat_convert_sp_gxf2gxf.pl`  |
| make feature **statistics** | `agat_sp_statistics.pl`  |
| make **function statistics** | `agat_sp_functional_statistics.pl`  |
| **extract** any type of sequence | `agat_sp_extract_sequences.pl`  |
| **extract** attributes | `agat_sp_extract_attributes.pl`  |
| **complement** annotations (non-overlapping loci) | `agat_sp_complement_annotations.pl`  |
| **merge** annotations | `agat_sp_merge_annotations.pl`  |
| **filter** gene models by ORF size | `agat_sp_filter_by_ORF_size.pl`  |
| **filter** to keep only longest isoforms | `agat_sp_keep_longest_isoform.pl`  |
| **create** introns features | `agat_sp_add_introns.pl`  |
| **fix** cds phases | `agat_sp_fix_cds_phases.pl`  |
| **manage** IDs | `agat_sp_manage_IDs.pl`  |
| **manage** UTRs | `agat_sp_manage_UTRs.pl`  |
| **manage** introns | `agat_sp_manage_introns.pl`  |
| **manage** functional annotation | `agat_sp_manage_functional_annotation.pl`  |
| **specificity sensitivity** | `agat_sp_sensitivity_specificity.pl`  |
| **fusion / split** analysis between two annotations | `agat_sp_compare_two_annotations.pl`  |
| analyze differences between **BUSCO** results | `agat_sp_compare_two_BUSCOs.pl`   |
| convert any **GTF/GFF** into **tabulated format** | `agat_sp_to_tabulated.pl`  |
| convert any **GTF/GFF** into **BED** format | `agat_convert_sp_gff2bed.pl`  |
| convert any **GTF/GFF** into **GTF** format | `agat_convert_sp_gff2gtf.pl`  |
| convert any **GTF/GFF** into any **GTF/GFF** (bioperl) format | `agat_convert_sp_gxf2gxf.pl`  |
| convert **BED** format into **GFF3** format | `agat_convert_bed2gff.pl`  |
| convert **EMBL** format into **GFF3** format | `agat_convert_embl2gff.pl`  |
| convert **genscan** format into **GFF3** format | `agat_convert_genscan2gff.pl`  |
| convert **mfannot** format into **GFF3** format | `agat_convert_mfannot2gff.pl`  |
| ... and much more ...| ... see [here](https://agat.readthedocs.io/en/latest/) ...|


## Installation

### Using conda

#### Install AGAT

  ```
  conda install -c bioconda agat
  ```

#### Update AGAT

  ```
  conda update agat
  ```

#### Uninstall AGAT
  ```
  conda uninstall agat  
  ```

### Old school - Manually

You will have to install all prerequisites and AGAT manually.

#### Install prerequisites
  * R  
    You can install it by conda (`conda install r-base`), through [CRAN](https://cran.r-project.org) ([See here for a nice tutorial](https://www.datacamp.com/community/tutorials/installing-R-windows-mac-ubuntu)) or using your package management tool (e.g apt for Debian, Ubuntu, and related Linux distributions).
  * Perl >= 5.8  
    It should already be available on your computer. If you are unlucky [perl.org](https://www.perl.org/get.html) is the place to go. 

  * Perl modules  
    They can be installed in different ways:
    
    * using cpan or cpanm
  
    ```
    cpanm install bioperl Clone Graph::Directed LWP::UserAgent Statistics::R JSON Carp Sort::Naturally File::Share File::ShareDir::Install Moose
    ```
    
    * using conda
    
      * using the provided yaml file
    
      ```
      conda env create -f conda_environment_AGAT.yml
      conda activate agat
      ``` 
    
      * manually  
    
      ```
      conda install perl-bioperl perl-clone perl-graph perl-lwp-simple perl-statistics-r perl-json perl-carp perl-sort-naturally perl-file-share perl-file-sharedir-install perl-moose
      ```
      
    * using your package management tool (e.g apt for Debian, Ubuntu, and related Linux distributions)
      
    ```
    apt install libbio-perl-perl libclone-perl libgraph-perl liblwp-useragent-determined-perl libstatistics-r-perl libjson-perl libcarp-clan-perl libsort-naturally-perl libfile-share-perl libfile-sharedir libfile-sharedir-install-perl
    ```

#### Install AGAT

  ```
  git clone https://github.com/NBISweden/AGAT.git # Clone AGAT
  cd AGAT                                         # move into AGAT folder
  perl Makefile.PL                                # Check all the dependencies*
  make                                            # Compile
  make test                                       # Test
  make install                                    # Install
  ```

<sup>*</sup>If dependencies are missing you will be warn. Please refer to the [Install prerequisites](#install-prerequisites) section.

**Remark**: On MS Windows, instead of make you'd probably have to use dmake or nmake depending the toolchain you have.

#### Update AGAT
From the folder where the repository is located.

  ```
  git pull                                        # Update to last AGAT
  perl Makefile.PL                                # Check all the dependencies*
  make                                            # Compile
  make test                                       # Test
  make install                                    # Install
  ```
<sup>*</sup>If dependencies are missing you will be warn. Please refer to the [Install prerequisites](#install-prerequisites) section.

#### Change to a specific version
From the folder where the repository is located.  

  ```
  git pull                                        # Update the code
  git checkout v0.1                               # use version v0.1 (See releases tab for a list of available versions)
  perl Makefile.PL                                # Check all the dependencies*
  make                                            # Compile
  make test                                       # Test
  make install                                    # Install
  ```
<sup>*</sup>If dependencies are missing you will be warn. Please refer to the [Install prerequisites](#install-prerequisites) section.

#### Uninstall AGAT

  ```
  perl uninstall_AGAT
  ```

## Usage

  ```
  script_name.pl -h
  ```

## List of tools
See [here](https://github.com/NBISweden/AGAT/wiki#list-of-agat-tools-v021) for a list of tools.  
As AGAT is a toolkit, it contains a lot of tools. The main one is `agat_convert_sp_gxf2gxf.pl` that allows to check, fix, pad missing information (features/attributes) of any kind of gtf and gff to create complete, sorted and standardised gff3 format.  
All the installed scripts have the `agat_` prefix.  
Typing `agat_` in your terminal followed by the <TAB> key to activate the autocompletion will display the complete list of available tool installed.

### More about the tools

#### with \_sp\_ prefix => Means SLURP

The gff file will be charged in memory Omniscient data structure that is way to facilitate access to desired features at any time.
It has a memory cost but make life smoother. Indeed, it allows to perform complicated tasks in a more time efficient way.
Moreover, it allows to fix all potential errors in the limit of the possibilities given by the format itself.
See the Omniscient section for more information about it.  

#### with \_sq\_ prefix => Means SEQUENTIAL

The gff file is read and processed from its top to the end line by line without sanity check. This is memory efficient.

## Omniscient - Standardisation for a full GFF3 compliant to any tool  

All tools with `agat_sp_` prefix will parse and slurps the entire data into a data structure called Omniscient.
Below you will find more information about peculiarity of the Omniscient structure,
and the parsing approach used.

#### Omniscient data structure

The method create a hash structure containing all the data in memory. We call it OMNISCIENT. The OMNISCIENT structure is a three levels structure:
```
$omniscient{level1}{tag_l1}{level1_id} = feature <= tag could be gene, match  
$omniscient{level2}{tag_l2}{idY} = @featureListL2 <= tag could be mRNA,rRNA,tRNA,etc. idY is a level1_id (know as Parent attribute within the level2 feature). The @featureList is a list to be able to manage isoform cases.  
$omniscient{level3}{tag_l3}{idZ} =  @featureListL3 <= tag could be exon,cds,utr3,utr5,etc. idZ is the ID of a level2 feature (know as Parent attribute within the level3 feature). The @featureList is a list to be able to put all the feature of a same tag together.  
```

#### How does the Omniscient parser work

The Omniscient parser phylosophy:
  * 1) Parse by Parent/child relationship  
  * 2) ELSE Parse by a common tag  (an attribute value shared by feature that must be grouped together. By default we are using locus_tag but can be set by parameter).  
  * 3) ELSE Parse sequentially (mean group features in a bucket, and the bucket change at each level2 feature, and bucket are join in a common tag at each new L1 feature).  

**/!\\** Case with only level3 features (i.e rast or some prokka files, sequential will not work as expected. Indeed all features will be the child of only one newly created Parent. To create a parent per feature or group of features, a common tag must be used to group them correctly. We use `gene_id` and `locus_tag` by default but you can set up the one of your choice)

To resume by priority of way to parse: **Parent/child relationship > locus_tag > sequential.**  
The parser may used only one or a mix of these approaches according of the peculiarity of the gtf/gff file you provide.

#### What can the Omniscient parser do for you

* It creates missing parental features. (e.g if a level2 or level3 feature do not have parental feature(s) we create the missing level2 and/or level1 feature(s)).    
* It creates missing mandatory attributes (ID and/or Parent).  
* It fixes identifier to be uniq.  
* It removes duplicated features (same position, same ID, same Parent).  
* It expands level3 features sharing multiple parents (e.g  if one exon has list of multiple parent mRNA in its Parent attribute, one exon per parent with uniq ID will be created.  
* It fixes feature location errors (e.g an mRNA spanning over its gene location, we fix the gene location).  
* It adds UTR if possible (CDS and exon present).  
* It adds exon if possible (CDS has to be present).  
* It groups features together (if related features are spread at different places in the file).  



#### examples
AGAT has been tested on 36 different peculiar GTF/GFF formats being different flavours or/and containing errors.
Below few are listed but you can find the full list of them into the `t/gff_syntax` directory.

example 8 - only CDS defined:  
```
##gff-version 3
Tob1_contig1	Prodigal:2.60	CDS	476	670	.	-	0	ID=Tob1_00001;locus_tag=Tob1_00001;product=hypothetical protein
Tob1_contig1	Prodigal:2.60	CDS	34266	35222	.	+	0	ID=Tob1_00024;locus_tag=Tob1_00024;product=hypothetical protein
Tob1_contig1	SignalP:4.1	sig_peptide	34266	34298	.	+	0	inference=ab initio prediction:SignalP:4.1;note=predicted cleavage at residue 33;product=putative signal peptide
Tob1_contig1	Prodigal:2.60	CDS	35267	37444	.	-	0	ID=Tob1_00025;locus_tag=Tob1_00025;
Tob1_contig1	SignalP:4.1	sig_peptide	37420	37444	.	-	0	inference=ab initio prediction:SignalP:4.1;note=predicted cleavage at residue 25;product=putative signal peptide
Tob1_contig1	Prodigal:2.60	CDS	38304	39338	.	-	0	ID=Tob1_00026;locus_tag=Tob1_00026;
```

`agat_convert_sp_gxf2gxf.pl --gff 8_test.gff`:  

```
##gff-version 3
Tob1_contig1	Prodigal:2.60	gene	476	670	.	-	0	ID=nbis_NEW-gene-1;locus_tag=Tob1_00001;product=hypothetical protein
Tob1_contig1	Prodigal:2.60	mRNA	476	670	.	-	0	ID=nbis_nol2id-cds-1;Parent=nbis_NEW-gene-1;locus_tag=Tob1_00001;product=hypothetical protein
Tob1_contig1	Prodigal:2.60	exon	476	670	.	-	.	ID=nbis_NEW-exon-1;Parent=nbis_nol2id-cds-1;locus_tag=Tob1_00001;product=hypothetical protein
Tob1_contig1	Prodigal:2.60	CDS	476	670	.	-	0	ID=Tob1_00001;Parent=nbis_nol2id-cds-1;locus_tag=Tob1_00001;product=hypothetical protein
Tob1_contig1	Prodigal:2.60	gene	34266	35222	.	+	0	ID=nbis_NEW-gene-2;locus_tag=Tob1_00024;product=hypothetical protein
Tob1_contig1	Prodigal:2.60	mRNA	34266	35222	.	+	0	ID=nbis_nol2id-cds-2;Parent=nbis_NEW-gene-2;locus_tag=Tob1_00024;product=hypothetical protein
Tob1_contig1	Prodigal:2.60	exon	34266	35222	.	+	.	ID=nbis_NEW-exon-2;Parent=nbis_nol2id-cds-2;locus_tag=Tob1_00024;product=hypothetical protein
Tob1_contig1	Prodigal:2.60	CDS	34266	35222	.	+	0	ID=Tob1_00024;Parent=nbis_nol2id-cds-2;locus_tag=Tob1_00024;product=hypothetical protein
Tob1_contig1	SignalP:4.1	sig_peptide	34266	34298	.	+	0	ID=sig_peptide-1;Parent=nbis_nol2id-cds-2;inference=ab initio prediction:SignalP:4.1;note=predicted cleavage at residue 33;product=putative signal peptide
Tob1_contig1	Prodigal:2.60	gene	35267	37444	.	-	0	ID=nbis_NEW-gene-3;locus_tag=Tob1_00025
Tob1_contig1	Prodigal:2.60	mRNA	35267	37444	.	-	0	ID=nbis_nol2id-cds-3;Parent=nbis_NEW-gene-3;locus_tag=Tob1_00025
Tob1_contig1	Prodigal:2.60	exon	35267	37444	.	-	.	ID=nbis_NEW-exon-3;Parent=nbis_nol2id-cds-3;locus_tag=Tob1_00025
Tob1_contig1	Prodigal:2.60	CDS	35267	37444	.	-	0	ID=Tob1_00025;Parent=nbis_nol2id-cds-3;locus_tag=Tob1_00025
Tob1_contig1	SignalP:4.1	sig_peptide	37420	37444	.	-	0	ID=sig_peptide-2;Parent=nbis_nol2id-cds-3;inference=ab initio prediction:SignalP:4.1;note=predicted cleavage at residue 25;product=putative signal peptide
Tob1_contig1	Prodigal:2.60	gene	38304	39338	.	-	0	ID=nbis_NEW-gene-4;locus_tag=Tob1_00026
Tob1_contig1	Prodigal:2.60	mRNA	38304	39338	.	-	0	ID=nbis_nol2id-cds-4;Parent=nbis_NEW-gene-4;locus_tag=Tob1_00026
Tob1_contig1	Prodigal:2.60	exon	38304	39338	.	-	.	ID=nbis_NEW-exon-4;Parent=nbis_nol2id-cds-4;locus_tag=Tob1_00026
Tob1_contig1	Prodigal:2.60	CDS	38304	39338	.	-	0	ID=Tob1_00026;Parent=nbis_nol2id-cds-4;locus_tag=Tob1_00026

```

example 9 - level2 feature missing (mRNA) and level3 features missing (UTRs):  
```
##gff-version 3
#!gff-spec-version 1.14
#!source-version NCBI C++ formatter 0.2
##Type DNA NC_003070.9
NC_003070.9	RefSeq	source	1	30427671	.	+	.	organism=Arabidopsis thaliana;mol_type=genomic DNA;db_xref=taxon:3702;chromosome=1;ecotype=Columbia
NC_003070.9	RefSeq	gene	3631	5899	.	+	.	ID=NC_003070.9:NAC001;locus_tag=AT1G01010;
NC_003070.9	RefSeq	exon	3631	3913	.	+	.	ID=NM_099983.2;Parent=NC_003070.9:NAC001;gbkey=mRNA;locus_tag=AT1G01010;
NC_003070.9	RefSeq	exon	3996	4276	.	+	.	ID=NM_099983.2;Parent=NC_003070.9:NAC001;gbkey=mRNA;locus_tag=AT1G01010;
NC_003070.9	RefSeq	exon	4486	4605	.	+	.	ID=NM_099983.2;Parent=NC_003070.9:NAC001;gbkey=mRNA;locus_tag=AT1G01010;
NC_003070.9	RefSeq	exon	4706	5095	.	+	.	ID=NM_099983.2;Parent=NC_003070.9:NAC001;gbkey=mRNA;locus_tag=AT1G01010;
NC_003070.9	RefSeq	exon	5174	5326	.	+	.	ID=NM_099983.2;Parent=NC_003070.9:NAC001;gbkey=mRNA;locus_tag=AT1G01010;
NC_003070.9	RefSeq	exon	5439	5899	.	+	.	ID=NM_099983.2;Parent=NC_003070.9:NAC001;gbkey=mRNA;locus_tag=AT1G01010;
NC_003070.9	RefSeq	CDS	3760	3913	.	+	0	ID=NM_099983.2;Parent=NC_003070.9:NAC001;locus_tag=AT1G01010;
NC_003070.9	RefSeq	CDS	3996	4276	.	+	2	ID=NM_099983.2;Parent=NC_003070.9:NAC001;locus_tag=AT1G01010;
NC_003070.9	RefSeq	CDS	4486	4605	.	+	0	ID=NM_099983.2;Parent=NC_003070.9:NAC001;locus_tag=AT1G01010;
NC_003070.9	RefSeq	CDS	4706	5095	.	+	0	ID=NM_099983.2;Parent=NC_003070.9:NAC001;locus_tag=AT1G01010;
NC_003070.9	RefSeq	CDS	5174	5326	.	+	0	ID=NM_099983.2;Parent=NC_003070.9:NAC001;locus_tag=AT1G01010;
NC_003070.9	RefSeq	CDS	5439	5627	.	+	0	ID=NM_099983.2;Parent=NC_003070.9:NAC001;locus_tag=AT1G01010;
NC_003070.9	RefSeq	start_codon	3760	3762	.	+	0	ID=NM_099983.2;Parent=NC_003070.9:NAC001;locus_tag=AT1G01010;
NC_003070.9	RefSeq	stop_codon	5628	5630	.	+	0	ID=NM_099983.2;Parent=NC_003070.9:NAC001;locus_tag=AT1G01010;
```

`agat_convert_sp_gxf2gxf.pl --gff 8_test.gff`:  

```
##gff-version 3
#!gff-spec-version 1.14
#!source-version NCBI C++ formatter 0.2
##Type DNA NC_003070.9
NC_003070.9	RefSeq	source	1	30427671	.	+	.	ID=source-1;chromosome=1;db_xref=taxon:3702;ecotype=Columbia;mol_type=genomic DNA;organism=Arabidopsis thaliana
NC_003070.9	RefSeq	gene	3631	5899	.	+	.	ID=nbis_NEW-gene-1;locus_tag=AT1G01010
NC_003070.9	RefSeq	mRNA	3631	5899	.	+	.	ID=NC_003070.9:NAC001;Parent=nbis_NEW-gene-1;locus_tag=AT1G01010
NC_003070.9	RefSeq	exon	3631	3913	.	+	.	ID=NM_099983.2;Parent=NC_003070.9:NAC001;gbkey=mRNA;locus_tag=AT1G01010
NC_003070.9	RefSeq	exon	3996	4276	.	+	.	ID=nbis_NEW-exon-1;Parent=NC_003070.9:NAC001;gbkey=mRNA;locus_tag=AT1G01010
NC_003070.9	RefSeq	exon	4486	4605	.	+	.	ID=nbis_NEW-exon-2;Parent=NC_003070.9:NAC001;gbkey=mRNA;locus_tag=AT1G01010
NC_003070.9	RefSeq	exon	4706	5095	.	+	.	ID=nbis_NEW-exon-3;Parent=NC_003070.9:NAC001;gbkey=mRNA;locus_tag=AT1G01010
NC_003070.9	RefSeq	exon	5174	5326	.	+	.	ID=nbis_NEW-exon-4;Parent=NC_003070.9:NAC001;gbkey=mRNA;locus_tag=AT1G01010
NC_003070.9	RefSeq	exon	5439	5899	.	+	.	ID=nbis_NEW-exon-5;Parent=NC_003070.9:NAC001;gbkey=mRNA;locus_tag=AT1G01010
NC_003070.9	RefSeq	CDS	3760	3913	.	+	0	ID=nbis_NEW-cds-1;Parent=NC_003070.9:NAC001;locus_tag=AT1G01010
NC_003070.9	RefSeq	CDS	3996	4276	.	+	2	ID=nbis_NEW-cds-1;Parent=NC_003070.9:NAC001;locus_tag=AT1G01010
NC_003070.9	RefSeq	CDS	4486	4605	.	+	0	ID=nbis_NEW-cds-1;Parent=NC_003070.9:NAC001;locus_tag=AT1G01010
NC_003070.9	RefSeq	CDS	4706	5095	.	+	0	ID=nbis_NEW-cds-1;Parent=NC_003070.9:NAC001;locus_tag=AT1G01010
NC_003070.9	RefSeq	CDS	5174	5326	.	+	0	ID=nbis_NEW-cds-1;Parent=NC_003070.9:NAC001;locus_tag=AT1G01010
NC_003070.9	RefSeq	CDS	5439	5627	.	+	0	ID=nbis_NEW-cds-1;Parent=NC_003070.9:NAC001;locus_tag=AT1G01010
NC_003070.9	RefSeq	five_prime_UTR	3631	3759	.	+	.	ID=nbis_NEW-five_prime_utr-1;Parent=NC_003070.9:NAC001;gbkey=mRNA;locus_tag=AT1G01010
NC_003070.9	RefSeq	start_codon	3760	3762	.	+	0	ID=nbis_NEW-start_codon-1;Parent=NC_003070.9:NAC001;locus_tag=AT1G01010
NC_003070.9	RefSeq	stop_codon	5628	5630	.	+	0	ID=nbis_NEW-stop_codon-1;Parent=NC_003070.9:NAC001;locus_tag=AT1G01010
NC_003070.9	RefSeq	three_prime_UTR	5628	5899	.	+	.	ID=nbis_NEW-three_prime_utr-1;Parent=NC_003070.9:NAC001;gbkey=mRNA;locus_tag=AT1G01010
```

example 18 - related features spread within the file:  
```
##gff-version 3
scaffold625	maker	gene	337818	343277	.	+	.	ID=CLUHARG00000005458;Name=TUBB3_2
scaffold625	maker	mRNA	337818	343277	.	+	.	ID=CLUHART00000008717;Parent=CLUHARG00000005458
scaffold625	maker	exon	337818	337971	.	+	.	ID=CLUHART00000008717:exon:1404;Parent=CLUHART00000008717
scaffold625	maker	exon	340733	340841	.	+	.	ID=CLUHART00000008717:exon:1405;Parent=CLUHART00000008717
scaffold789	maker	three_prime_UTR	564589	564780	.	+	.	ID=CLUHART00000006146:three_prime_utr;Parent=CLUHART00000006146
scaffold789	maker	mRNA	558184	564780	.	+	.	ID=CLUHART00000006147;Parent=CLUHARG00000003852
scaffold625	maker	CDS	337915	337971	.	+	0	ID=CLUHART00000008717:cds;Parent=CLUHART00000008717
scaffold625	maker	CDS	340733	340841	.	+	0	ID=CLUHART00000008717:cds;Parent=CLUHART00000008717
scaffold625	maker	CDS	341518	341628	.	+	2	ID=CLUHART00000008717:cds;Parent=CLUHART00000008717
scaffold625	maker	CDS	341964	343033	.	+	2	ID=CLUHART00000008717:cds;Parent=CLUHART00000008717
scaffold625	maker	five_prime_UTR	337818	337914	.	+	.	ID=CLUHART00000008717:five_prime_utr;Parent=CLUHART00000008717
scaffold625	maker	three_prime_UTR	343034	343277	.	+	.	ID=CLUHART00000008717:three_prime_utr;Parent=CLUHART00000008717
scaffold789	maker	gene	558184	564780	.	+	.	ID=CLUHARG00000003852;Name=PF11_0240
scaffold789	maker	mRNA	558184	564780	.	+	.	ID=CLUHART00000006146;Parent=CLUHARG00000003852
scaffold789	maker	exon	558184	560123	.	+	.	ID=CLUHART00000006146:exon:995;Parent=CLUHART00000006146
scaffold789	maker	exon	561401	561519	.	+	.	ID=CLUHART00000006146:exon:996;Parent=CLUHART00000006146
scaffold789	maker	exon	564171	564235	.	+	.	ID=CLUHART00000006146:exon:997;Parent=CLUHART00000006146
scaffold789	maker	exon	564372	564780	.	+	.	ID=CLUHART00000006146:exon:998;Parent=CLUHART00000006146
scaffold789	maker	CDS	558191	560123	.	+	0	ID=CLUHART00000006146:cds;Parent=CLUHART00000006146
scaffold789	maker	CDS	561401	561519	.	+	2	ID=CLUHART00000006146:cds;Parent=CLUHART00000006146
scaffold625	maker	exon	341518	341628	.	+	.	ID=CLUHART00000008717:exon:1406;Parent=CLUHART00000008717
scaffold625	maker	exon	341964	343277	.	+	.	ID=CLUHART00000008717:exon:1407;Parent=CLUHART00000008717
scaffold789	maker	CDS	564171	564235	.	+	0	ID=CLUHART00000006146:cds;Parent=CLUHART00000006146
scaffold789	maker	CDS	564372	564588	.	+	1	ID=CLUHART00000006146:cds;Parent=CLUHART00000006146
scaffold789	maker	five_prime_UTR	558184	558190	.	+	.	ID=CLUHART00000006146:five_prime_utr;Parent=CLUHART00000006146
scaffold789	maker	exon	558184	560123	.	+	.	ID=CLUHART00000006147:exon:997;Parent=CLUHART00000006147
scaffold789	maker	exon	561401	561519	.	+	.	ID=CLUHART00000006147:exon:998;Parent=CLUHART00000006147
scaffold789	maker	exon	562057	562121	.	+	.	ID=CLUHART00000006147:exon:999;Parent=CLUHART00000006147
scaffold789	maker	exon	564372	564780	.	+	.	ID=CLUHART00000006147:exon:1000;Parent=CLUHART00000006147
scaffold789	maker	CDS	558191	560123	.	+	0	ID=CLUHART00000006147:cds;Parent=CLUHART00000006147
scaffold789	maker	CDS	561401	561519	.	+	2	ID=CLUHART00000006147:cds;Parent=CLUHART00000006147
scaffold789	maker	CDS	562057	562121	.	+	0	ID=CLUHART00000006147:cds;Parent=CLUHART00000006147
scaffold789	maker	CDS	564372	564588	.	+	1	ID=CLUHART00000006147:cds;Parent=CLUHART00000006147
scaffold789	maker	five_prime_UTR	558184	558190	.	+	.	ID=CLUHART00000006147:five_prime_utr;Parent=CLUHART00000006147
scaffold789	maker	three_prime_UTR	564589	564780	.	+	.	ID=CLUHART00000006147:three_prime_utr;Parent=CLUHART00000006147
```

`agat_convert_sp_gxf2gxf.pl --gff 18_test.gff`:  
```
##gff-version 3
scaffold625	maker	gene	337818	343277	.	+	.	ID=CLUHARG00000005458;Name=TUBB3_2
scaffold625	maker	mRNA	337818	343277	.	+	.	ID=CLUHART00000008717;Parent=CLUHARG00000005458
scaffold625	maker	exon	337818	337971	.	+	.	ID=CLUHART00000008717:exon:1404;Parent=CLUHART00000008717
scaffold625	maker	exon	340733	340841	.	+	.	ID=CLUHART00000008717:exon:1405;Parent=CLUHART00000008717
scaffold625	maker	exon	341518	341628	.	+	.	ID=CLUHART00000008717:exon:1406;Parent=CLUHART00000008717
scaffold625	maker	exon	341964	343277	.	+	.	ID=CLUHART00000008717:exon:1407;Parent=CLUHART00000008717
scaffold625	maker	CDS	337915	337971	.	+	0	ID=CLUHART00000008717:cds;Parent=CLUHART00000008717
scaffold625	maker	CDS	340733	340841	.	+	0	ID=CLUHART00000008717:cds;Parent=CLUHART00000008717
scaffold625	maker	CDS	341518	341628	.	+	2	ID=CLUHART00000008717:cds;Parent=CLUHART00000008717
scaffold625	maker	CDS	341964	343033	.	+	2	ID=CLUHART00000008717:cds;Parent=CLUHART00000008717
scaffold625	maker	five_prime_UTR	337818	337914	.	+	.	ID=CLUHART00000008717:five_prime_utr;Parent=CLUHART00000008717
scaffold625	maker	three_prime_UTR	343034	343277	.	+	.	ID=CLUHART00000008717:three_prime_utr;Parent=CLUHART00000008717
scaffold789	maker	gene	558184	564780	.	+	.	ID=CLUHARG00000003852;Name=PF11_0240
scaffold789	maker	mRNA	558184	564780	.	+	.	ID=CLUHART00000006146;Parent=CLUHARG00000003852
scaffold789	maker	exon	558184	560123	.	+	.	ID=CLUHART00000006146:exon:995;Parent=CLUHART00000006146
scaffold789	maker	exon	561401	561519	.	+	.	ID=CLUHART00000006146:exon:996;Parent=CLUHART00000006146
scaffold789	maker	exon	564171	564235	.	+	.	ID=CLUHART00000006146:exon:997;Parent=CLUHART00000006146
scaffold789	maker	exon	564372	564780	.	+	.	ID=CLUHART00000006146:exon:998;Parent=CLUHART00000006146
scaffold789	maker	CDS	558191	560123	.	+	0	ID=CLUHART00000006146:cds;Parent=CLUHART00000006146
scaffold789	maker	CDS	561401	561519	.	+	2	ID=CLUHART00000006146:cds;Parent=CLUHART00000006146
scaffold789	maker	CDS	564171	564235	.	+	0	ID=CLUHART00000006146:cds;Parent=CLUHART00000006146
scaffold789	maker	CDS	564372	564588	.	+	1	ID=CLUHART00000006146:cds;Parent=CLUHART00000006146
scaffold789	maker	five_prime_UTR	558184	558190	.	+	.	ID=CLUHART00000006146:five_prime_utr;Parent=CLUHART00000006146
scaffold789	maker	three_prime_UTR	564589	564780	.	+	.	ID=CLUHART00000006146:three_prime_utr;Parent=CLUHART00000006146
scaffold789	maker	mRNA	558184	564780	.	+	.	ID=CLUHART00000006147;Parent=CLUHARG00000003852
scaffold789	maker	exon	558184	560123	.	+	.	ID=CLUHART00000006147:exon:997;Parent=CLUHART00000006147
scaffold789	maker	exon	561401	561519	.	+	.	ID=CLUHART00000006147:exon:998;Parent=CLUHART00000006147
scaffold789	maker	exon	562057	562121	.	+	.	ID=CLUHART00000006147:exon:999;Parent=CLUHART00000006147
scaffold789	maker	exon	564372	564780	.	+	.	ID=CLUHART00000006147:exon:1000;Parent=CLUHART00000006147
scaffold789	maker	CDS	558191	560123	.	+	0	ID=CLUHART00000006147:cds;Parent=CLUHART00000006147
scaffold789	maker	CDS	561401	561519	.	+	2	ID=CLUHART00000006147:cds;Parent=CLUHART00000006147
scaffold789	maker	CDS	562057	562121	.	+	0	ID=CLUHART00000006147:cds;Parent=CLUHART00000006147
scaffold789	maker	CDS	564372	564588	.	+	1	ID=CLUHART00000006147:cds;Parent=CLUHART00000006147
scaffold789	maker	five_prime_UTR	558184	558190	.	+	.	ID=CLUHART00000006147:five_prime_utr;Parent=CLUHART00000006147
scaffold789	maker	three_prime_UTR	564589	564780	.	+	.	ID=CLUHART00000006147:three_prime_utr;Parent=CLUHART00000006147
```

## How to cite?

This work has not been published (I will think about it). But if you wish to cite AGAT you could probably do it as follow (Adapt the version for the one you have used): 

```
Dainat J. AGAT: Another Gff Analysis Toolkit to handle annotations in any GTF/GFF format.  
(Version v0.4.0). Zenodo. https://www.doi.org/10.5281/zenodo.3552717
```

## Publication using AGAT
Some examples of publications that have used AGAT
  * [A long read optimized de novo transcriptome pipeline reveals novel ocular developmentally regulated gene isoforms and disease targets](https://www.biorxiv.org/content/10.1101/2020.08.21.261644v2.full.pdf)
  * [A telomere to telomere assembly of Oscheius tipulae and the evolution of rhabditid nematode
chromosomes](https://www.biorxiv.org/content/10.1101/2020.09.04.283127v1.full.pdf)
  * [In vitro resynthesis of lichenization reveals the genetic background of symbiosis-specific fungal-algal interaction in Usnea hakonensis](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-020-07086-9)
  * [Application of an optimized annotation pipeline to the Cryptococcus deuterogattii genome reveals dynamic primary metabolic gene clusters and genomic impact of RNAi loss](https://www.biorxiv.org/content/10.1101/2020.09.01.278374v1.full)
  * [Genomics of an avian neo-sex chromosome reveals the evolutionary dynamics of recombination suppression and sex-linked genes](https://www.biorxiv.org/content/10.1101/2020.09.25.314088v1.full)
  * [Four novel Picornaviruses detected in Magellanic Penguins (Spheniscus magellanicus) in Chile](https://www.biorxiv.org/content/10.1101/2020.10.26.356485v1.full.pdf)
  * [The Crown Pearl: a draft genome assembly of the European freshwater pearl mussel Margaritifera margaritifera (Linnaeus, 1758)](https://www.biorxiv.org/content/10.1101/2020.12.06.413450v1.full)
  * [Investigating the impact of reference assembly choice on genomic analyses in a cattle breed](https://www.biorxiv.org/content/10.1101/2021.01.15.426838v1.full.pdf)
  * [Two novel loci underlie natural differences in Caenorhabditis elegans abamectin responses](https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1009297)
  * [Butterfly eyespots evolved via co-option of the antennal gene-regulatory network](https://www.biorxiv.org/content/10.1101/2021.03.01.429915v2.full)
  * [Transcript- and annotation-guided genome assembly of the European starling](https://www.biorxiv.org/content/10.1101/2021.04.07.438753v1)
  * [...]

## Troubleshooting

### AGAT throws features out, because the feature type is not yet taken into account
Feature types (primary_tag) handled by AGAT are defined within json files. Most common features are already defined in those files. If you encounter files with feature types not accepted, AGAT will inform you and throw the features out. To keep those feature you must inform properly AGAT how to handle them.
First access the json files by running:
```
			agat_convert_sp_gxf2gxf.pl --expose
```

Then open the file corresponding to the feature type you want to add:
* Feature level1 (e.g. gene, match, region):
  My feature has no parent
  => features_level1.json
* Feature level2 (e.g. mrna, match_part, trna):
  My feature has one parent and the parent is a level 1 feature.
  => features_level2.json.
* Feature level3 (e.g. exon, intron, cds):
  My feature has one parent (the parent has also a parent) and no children
  => features_level3.json.
* Feature level3 discontinuous (e.g. cds, utr):
  A single feature that exists over multiple genomic locations
  => features_spread.json.

Then add the feature type information by adding a paired-value like this:
```
	"bigRNA":"gene",
```
Where `bigRNA`is the feature type and `gene`the parent feature type expected.
/!\\ For level1 feature type the second value can be:
 * topfeature: feature does not expect children, and will be written first in the sequence
 * standalone: feature does not expect children
 * other values do not have any meaning but a value is required, write whatever you want. 

### AGAT throws features out, because child features are not provided
Features level1 (e.g. gene, match, chromosome) may require to have child features or not depending of the information stored into the `features_level1.json` file. If a child is required, and the GFF file does not contain it, the level1 feature will be thrown away. You must modify the json file to add the the term `standalone` to inform AGAT that this feature level1 do not require any child. (This work only on feature level1, not level2 or level3). To access the json files run the following command:
```
# export the json files
agat_convert_sp_gxf2gxf.pl --expose
```
Then open the `features_level1.json` and put the value `standalone` as value to the required feature.
Finally run your scripts in the same folder as the modified json files are standing.

### Use a version of AGAT from a specific branch
```
# install AGAT dependencies
conda install -c bioconda agat
# clone the repo
git clone https://github.com/NBISweden/AGAT.git
# if the branch you want is not the master (replace BRANCHE_NAME by the one you wish to use)
git checkout BRANCHE_NAME
# move into AGAT folder
cd AGAT 
# Check all the dependencies*
perl Makefile.PL
# Compile
make
# Test
make test
# Install
make install                                    
```
