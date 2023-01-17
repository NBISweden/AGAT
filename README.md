[![Build Status](https://travis-ci.org/NBISweden/AGAT.svg?branch=master)](https://travis-ci.org/NBISweden/AGAT)
[![Coverage Status](https://coveralls.io/repos/github/NBISweden/AGAT/badge.svg)](https://coveralls.io/github/NBISweden/AGAT)
[![Documentation Status](https://readthedocs.org/projects/agat/badge/?version=latest)](https://agat.readthedocs.io/en/latest/?badge=latest)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/agat/README.html)
[<img alt="docker_agat" src="https://img.shields.io/badge/container-Docker-blue">](https://quay.io/repository/biocontainers/agat)
[<img alt="singularity_agat" src="https://img.shields.io/badge/container-Singularity-orange">](https://quay.io/repository/biocontainers/agat)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/agat/badges/license.svg)](https://anaconda.org/bioconda/agat)
[![Anaconda-Server Badge](https://img.shields.io/conda/dn/bioconda/agat.svg?style=flat)](https://anaconda.org/bioconda/agat)
[![DOI](https://zenodo.org/badge/222659741.svg)](https://zenodo.org/badge/latestdoi/222659741)  

AGAT
=========================================
[<img align="right" src="docs/img/NBIS.png" width="200" height="100" />](https://nbis.se)
<h2><em>A</em>nother <em>G</em>tf/Gff <em>A</em>nalysis <i>T</i>oolkit</h2>  

Suite of tools to handle gene annotations in any GTF/GFF format.  
\>>[docs](https://agat.readthedocs.io/en/latest/index.html)<<

---------------------------

## Table of Contents

   * [What can AGAT do for you?](#what-can-agat-do-for-you)
   * [Installation](#installation)
       * [Using Docker](#using-docker)
       * [Using Singularity](#using-singularity)
       * [Using Bioconda](#using-bioconda)
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
   * [The AGAT parser - Standardisation to create GXF files compliant to any tool](#the-agat-parser---standardisation-to-create-gxf-files-compliant-to-any-tool)
      * [The data structure](#the-data-structure)
      * [How does the AGAT parser work](#how-does-the-agat-parser-work)
      * [What can the AGAT parser do for you](#what-can-the-agat-parser-do-for-you)
      * [examples](#examples)
   * [How to cite?](#how-to-cite)
   * [Publication using AGAT](#publication-using-agat)
   * [Troubleshooting](#troubleshooting)
---------------------------

## What can AGAT do for you?  

AGAT has the power to check, fix, pad missing information (features/attributes) of any kind of GTF and GFF to create complete, sorted and standardised gff3 format. Over the years it has been enriched by many many tools to perform just about any tasks that is possible related to GTF/GFF format files (sanitizing, conversions, merging, modifying, filtering, FASTA sequence extraction, adding information, etc). Comparing to other methods AGAT is robust to even the most despicable GTF/GFF files.

  * Standardize/sanitize any GTF/GFF file into a comprehensive GFF3 format (script with `_sp_` prefix)

    <details>
      <summary>See standardization/sanitization tool</summary>

      | task | tool |
      | --- | --- |
      | **check, fix, pad** missing information into sorted and standardised gff3 | `agat_convert_sp_gxf2gxf.pl`  |

      * add missing parent features (e.g. gene and mRNA if only CDS/exon exists).  
      * add missing features (e.g. exon and UTR).  
      * add missing mandatory attributes (i.e. ID, Parent).  
      * fix identifiers to be uniq.  
      * fix feature locations.  
      * remove duplicated features.  
      * group related features (if spread in different places in the file).  
      * sort features (tabix optional).  
      * merge overlapping loci into one single locus (only if option activated).  
    </details>

  * Convert many formats

    <details>
      <summary>See conversion tools</summary>

      | task | tool |
      | --- | --- |
      | convert any **GTF/GFF** into **BED** format | `agat_convert_sp_gff2bed.pl`  |
      | convert any **GTF/GFF** into **GTF** format | `agat_convert_sp_gff2gtf.pl`  |
      | convert any **GTF/GFF** into **tabulated format** | `agat_sp_gff2tsv.pl`  |
      | convert any **BAM** from minimap2 into **GFF** format | `agat_convert_sp_minimap2_bam2gff.pl`  |
      | convert any **GTF/GFF** into **ZFF** format | `agat_sp_gff2zff.pl`  |
      | convert any **GTF/GFF** into any **GTF/GFF** (bioperl) format | `agat_convert_sp_gxf2gxf.pl`  |
      | convert **BED** format into **GFF3** format | `agat_convert_bed2gff.pl`  |
      | convert **EMBL** format into **GFF3** format | `agat_convert_embl2gff.pl`  |
      | convert **genscan** format into **GFF3** format | `agat_convert_genscan2gff.pl`  |
      | convert **mfannot** format into **GFF3** format | `agat_convert_mfannot2gff.pl`  |
    </details>


  * Perform numerous tasks (Just about anything that is possible)

    <details>
      <summary>See tools</summary>

      | task | tool |
      | --- | --- |
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
      | ... and much more ...| ... see [here](https://agat.readthedocs.io/en/latest/) ...|
    </details>

**About the GTF/GFF fromat**  
The GTF/GFF formats are 9-column text formats used to describe and represent genomic features.
The formats have quite evolved since 1997, and despite well-defined specifications existing nowadays they have a great flexibility allowing holding wide variety of information.
This flexibility has a drawback aspect, there is an incredible amount of flavour of the formats, that can result in problems when using downstream programs.  
For a complete overview of the GTF/GFF formats have a look [here](https://agat.readthedocs.io/en/latest/gxf.html).

## Installation

### Using Docker
   <details>
      <summary>See details</summary>
First you must have [Docker](https://docs.docker.com/get-docker/) installed and running.  
Secondly have look at the availabe AGAT biocontainers at [quay.io](https://quay.io/repository/biocontainers/agat?tab=tags).  
Then:
  ```
# get the chosen AGAT container version
docker pull quay.io/biocontainers/agat:0.8.0--pl5262hdfd78af_0
# use an AGAT's tool e.g. agat_convert_sp_gxf2gxf.pl
docker run quay.io/biocontainers/agat:0.8.0--pl5262hdfd78af_0 agat_convert_sp_gxf2gxf.pl --help
  ```
   </details>
 
### Using Singularity
   <details>
      <summary>See details</summary>
First you must have [Singularity](https://sylabs.io/guides/3.5/user-guide/quick_start.html) installed and running.
Secondly have look at the availabe AGAT biocontainers at [quay.io](https://quay.io/repository/biocontainers/agat?tab=tags).  
Then:
```
# get the chosen AGAT container version
singularity pull docker://quay.io/biocontainers/agat:0.8.0--pl5262hdfd78af_0
# run the container
singularity run agat_0.8.0--pl5262hdfd78af_0.sif
```

You are now in the container. You can use an AGAT's tool e.g. agat_convert_sp_gxf2gxf.pl doing
```
agat_convert_sp_gxf2gxf.pl --help
```
   </details>

### Using Bioconda
   <details>
      <summary>See details</summary>
      
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
   </details>
   
### Old school - Manually
   <details>
      <summary>See details</summary>
You will have to install all prerequisites and AGAT manually.

#### Install prerequisites
  * R (optional)  
    You can install it by conda (`conda install r-base`), through [CRAN](https://cran.r-project.org) ([See here for a nice tutorial](https://www.datacamp.com/community/tutorials/installing-R-windows-mac-ubuntu)) or using your package management tool (e.g apt for Debian, Ubuntu, and related Linux distributions). R is optional and can be used to perform some plots. You will need to install the perl depency Statistics::R

  * Perl >= 5.8  
    It should already be available on your computer. If you are unlucky [perl.org](https://www.perl.org/get.html) is the place to go.

  * Perl modules  
    They can be installed in different ways:

    * using cpan or cpanm

    ```
    cpanm install bioperl Clone Graph::Directed LWP::UserAgent Carp Sort::Naturally File::Share File::ShareDir::Install Moose YAML LWP::Protocol::https
    ```

    * using conda

      * using the provided yaml file

      ```
      conda env create -f conda_environment_AGAT.yml
      conda activate agat
      ```

      * manually  

      ```
      conda install perl-bioperl perl-clone perl-graph perl-lwp-simple perl-carp perl-sort-naturally perl-file-share perl-file-sharedir-install perl-moose perl-yaml perl-lwp-protocol-https
      ```

    * using your package management tool (e.g apt for Debian, Ubuntu, and related Linux distributions)

    ```
    apt install libbio-perl-perl libclone-perl libgraph-perl liblwp-useragent-determined-perl libstatistics-r-perl libcarp-clan-perl libsort-naturally-perl libfile-share-perl libfile-sharedir libfile-sharedir-install-perl libyaml-perl liblwp-protocol-https-perl
    ```

  * Optional
    Some scripts offer the possibility to perform plots. You will need R and Statistics::R which are not included by default.

    * R   
      You can install it by conda (`conda install r-base`), through [CRAN](https://cran.r-project.org) ([See here for a nice tutorial](https://www.datacamp.com/community/tutorials/installing-R-windows-mac-ubuntu)) or using your package management tool (e.g apt for Debian, Ubuntu, and related Linux distributions).

    * Statistics::R
        You can install it through conda  (`conda install perl-statistics-r`), using cpan/cpanm (`cpanm install Statistics::R`), or your package management tool  (`apt install libstatistics-r-perl`)



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
   </details>

## Usage

  ```
  script_name.pl -h
  ```

## List of tools

As AGAT is a toolkit, it contains a lot of tools. The main one is `agat_convert_sp_gxf2gxf.pl` that allows to check, fix, pad missing information (features/attributes) of any kind of gtf and gff to create complete, sorted and standardised gff3 format.
All the installed scripts have the `agat_` prefix.  

To have a look to the available tools you have several approaches:  
  * `agat --tools`
  * Typing `agat_` in your terminal followed by the <TAB> key to activate the autocompletion will display the complete list of available tools.
  * [The documentation](https://agat.readthedocs.io/en/latest/?badge=latest).  


### More about the tools

#### with \_sp\_ prefix => Means SLURP

The gff file will be charged in memory in a specific data structure facilitating the access to desired features at any time.
It has a memory cost but make life smoother. Indeed, it allows to perform complicated tasks in a more time efficient way.
Moreover, it allows to fix all potential errors in the limit of the possibilities given by the format itself.
See the AGAT parser section for more information about it.  

#### with \_sq\_ prefix => Means SEQUENTIAL

The gff file is read and processed from its top to the end line by line without sanity check. This is memory efficient.

## The AGAT parser - Standardisation to create GXF files compliant to any tool  

All tools with `agat_sp_` prefix will parse and slurps the entire data into a specific data structure called.
Below you will find more information about peculiarity of the data structure,
and the parsing approach used.

#### the data structure

The method create a hash structure containing all the data in memory. We can call it OMNISCIENT. The OMNISCIENT structure is a three levels structure:
```
$omniscient{level1}{tag_l1}{level1_id} = feature <= tag could be gene, match  
$omniscient{level2}{tag_l2}{idY} = @featureListL2 <= tag could be mRNA,rRNA,tRNA,etc. idY is a level1_id (know as Parent attribute within the level2 feature). The @featureListL2 is a list to be able to manage isoform cases.  
$omniscient{level3}{tag_l3}{idZ} =  @featureListL3 <= tag could be exon,cds,utr3,utr5,etc. idZ is the ID of a level2 feature (know as Parent attribute within the level3 feature). The @featureListL3 is a list to be able to put all the feature of a same tag together.  
```

#### How does the AGAT parser work

The AGAT parser phylosophy:
  * 1) Parse by Parent/child relationship or gene_id/transcript_id relationship.
  * 2) ELSE Parse by a common tag  (an attribute value shared by feature that must be grouped together. By default we are using locus_tag but can be set by parameter).  
  * 3) ELSE Parse sequentially (mean group features in a bucket, and the bucket change at each level2 feature, and bucket are join in a common tag at each new L1 feature).  

**/!\\** Case with only level3 features (i.e rast or some prokka files, sequential will not work as expected. Indeed all features will be the child of only one newly created Parent. To create a parent per feature or group of features, a common tag must be used to group them correctly. We use `gene_id` and `locus_tag` by default but you can set up the one of your choice)

To resume by priority of way to parse: **Parent/child relationship > locus_tag > sequential.**  
The parser may used only one or a mix of these approaches according of the peculiarity of the gtf/gff file you provide.

#### What can the AGAT parser do for you

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
AGAT was tested on 42 different types of GTF/GFF of different flavours or/and containing errors.
Below few are listed but you can find the full list of them into the `t/gff_syntax` directory.

##### example 8 - only CDS defined

<details>
  <summary>See example</summary>

```
##gff-version 3
Tob1_contig1	Prodigal:2.60	CDS	476	670	.	-	0	ID=Tob1_00001;locus_tag=Tob1_00001;product=hypothetical protein
Tob1_contig1	Prodigal:2.60	CDS	34266	35222	.	+	0	ID=Tob1_00024;locus_tag=Tob1_00024;product=hypothetical protein
Tob1_contig1	SignalP:4.1	sig_peptide	34266	34298	.	+	0	inference=ab initio prediction:SignalP:4.1;note=predicted cleavage at residue 33;product=putative signal peptide
Tob1_contig1	Prodigal:2.60	CDS	35267	37444	.	-	0	ID=Tob1_00025;locus_tag=Tob1_00025;
Tob1_contig1	SignalP:4.1	sig_peptide	37420	37444	.	-	0	inference=ab initio prediction:SignalP:4.1;note=predicted cleavage at residue 25;product=putative signal peptide
Tob1_contig1	Prodigal:2.60	CDS	38304	39338	.	-	0	ID=Tob1_00026;locus_tag=Tob1_00026;
```
</details>

`agat_convert_sp_gxf2gxf.pl --gff 8_test.gff`  

<details>
  <summary>See result</summary>

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
</details>

##### example 9 - level2 feature missing (mRNA) and level3 features missing (UTRs)

<details>
  <summary>See example</summary>

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
</details>

`agat_convert_sp_gxf2gxf.pl --gff 8_test.gff`  

<details>
  <summary>See result</summary>

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
</details>

##### example 18 - related features spread within the file  

<details>
  <summary>See example</summary>

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
</details>

`agat_convert_sp_gxf2gxf.pl --gff 18_test.gff`

<details>
  <summary>See result</summary>

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
</details>

## How to cite?

This work has not been published (I will think about it). But if you wish to cite AGAT you could probably do it as follow (Adapt the version for the one you have used):

```
Dainat J. AGAT: Another Gff Analysis Toolkit to handle annotations in any GTF/GFF format.  
(Version v0.8.0). Zenodo. https://www.doi.org/10.5281/zenodo.3552717
```

## Publication using AGAT

Some examples of publications that have used AGAT  

<details>
<summary>See publications</summary>

| Journal | Title |
| --- | --- |
| Genome Biology and Evolution | [Ancestral Physical Stress and Later Immune Gene Family Expansions Shaped Bivalve Mollusc Evolution](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8382680/)
| Preprint | [A long read optimized de novo transcriptome pipeline reveals novel ocular developmentally regulated gene isoforms and disease targets](https://www.biorxiv.org/content/10.1101/2020.08.21.261644v2.full.pdf)
| G3 Genes Genomes Genetics | [A telomere to telomere assembly of Oscheius tipulae and the evolution of rhabditid nematode chromosomes](https://academic.oup.com/g3journal/article/11/1/jkaa020/6026964)
| BMC genomics | [In vitro resynthesis of lichenization reveals the genetic background of symbiosis-specific fungal-algal interaction in Usnea hakonensis](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-020-07086-9)
| G3 Genes Genomes Genetics | [Application of an optimized annotation pipeline to the Cryptococcus deuterogattii genome reveals dynamic primary metabolic gene clusters and genomic impact of RNAi loss](https://www.biorxiv.org/content/10.1101/2020.09.01.278374v1.full)
| Mol. Biol. Evol. | [Genomics of an avian neo-sex chromosome reveals the evolutionary dynamics of recombination suppression and sex-linked genes](https://academic.oup.com/mbe/advance-article/doi/10.1093/molbev/msab277/6372697)
| Virology | [Four novel Picornaviruses detected in Magellanic Penguins (Spheniscus magellanicus) in Chile](https://www.sciencedirect.com/science/article/pii/S0042682221001148)
| DNA Research | [The Crown Pearl: a draft genome assembly of the European freshwater pearl mussel Margaritifera margaritifera (Linnaeus, 1758)](https://academic.oup.com/dnaresearch/article/28/2/dsab002/6182681)
| BMC genomics | [Investigating the impact of reference assembly choice on genomic analyses in a cattle breed](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-021-07554-w)
| Plos pathogens | [Two novel loci underlie natural differences in Caenorhabditis elegans abamectin responses](https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1009297)
| Preprint | [Butterfly eyespots evolved via co-option of the antennal gene-regulatory network](https://www.biorxiv.org/content/10.1101/2021.03.01.429915v2.full)
| Preprint | [Transcript- and annotation-guided genome assembly of the European starling](https://www.biorxiv.org/content/10.1101/2021.04.07.438753v1)
| Microbiol Resour Announc. | [LGAAP: Leishmaniinae Genome Assembly and Annotation Pipeline](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8297458/)
| Genome Biology and Evolution | [A Chromosome-level Genome Assembly of the Reed Warbler (Acrocephalus scirpaceus) ](https://academic.oup.com/gbe/article/13/9/evab212/6367782)
| Preprint | [Barcoded RH-seq illuminates the complex genetic basis of yeast thermotolerance](https://www.biorxiv.org/content/10.1101/2021.07.26.453780v1.full)
| Gygabyte | [A high-quality draft genome for Melaleuca alternifolia (tea tree): a new platform for evolutionary genomics of myrtaceous terpene-rich species](https://gigabytejournal.com/articles/28)
| Nature | [Chromosome-scale genome sequencing, assembly and annotation of six genomes from subfamily Leishmaniinae](https://www.nature.com/articles/s41597-021-01017-3#citeas)
| Preprint |[High quality, phased genomes of Phytophthora ramorum clonal lineages NA1 and EU1](https://www.biorxiv.org/content/10.1101/2021.06.23.449625v1.full)
| Elife | [Analysis of meiosis in Pristionchus pacificus reveals plasticity in homolog pairing and synapsis in the nematode lineage](https://elifesciences.org/articles/70990)
| MDPI | [Transcriptome Comparison of Secondary Metabolite Biosynthesis Genes Expressed in Cultured and Lichenized Conditions of Cladonia rangiferina](https://www.mdpi.com/1424-2818/13/11/529/html)
| MDPI | [FA-nf: A Functional Annotation Pipeline for Proteins from Non-Model Organisms Implemented in Nextflow](https://www.mdpi.com/2073-4425/12/10/1645/htm)
| Preprint | [De Novo Whole Genome Assembly of the Roborovski Dwarf Hamster (Phodopus roborovskii) Genome, an Animal Model for Severe/Critical COVID-19](https://www.biorxiv.org/content/10.1101/2021.10.02.462569v2.full)
| Preprint | [Using historical museum samples to examine divergent and parallel evolution in the invasive starling](https://www.biorxiv.org/content/10.1101/2021.08.22.457241v1.full)|
| GBE | [A Chromosome-Level Genome Assembly of the Reed Warbler (Acrocephalus scirpaceus)](https://helda.helsinki.fi/bitstream/handle/10138/336322/evab212.pdf?sequence=1&isAllowed=y)|
| Preprint | [A genome assembly of the Atlantic chub mackerel (Scomber colias): a valuable teleost fishing resource](https://www.biorxiv.org/content/10.1101/2021.11.19.468211v1.full.pdf)|
| Current Protocols | [BUSCO: Assessing Genomic Data Quality and Beyond](https://currentprotocols.onlinelibrary.wiley.com/doi/full/10.1002/cpz1.323)
| [...] | [...]
</details>

## Troubleshooting

See Troubleshooting section form the doc [here](https://agat.readthedocs.io/en/latest/troubleshooting.html).
