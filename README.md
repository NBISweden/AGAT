

[<img align="center" src="NBIS.png" width="200" height="100" />](https://nbis.se)
<h2><em>A</em>nother <em>G</em>ff <em>A</em>nalysis <i>T</i>oolkit (AGAT)</h2>  
Suite of tools to handle gene annotations in any GTF/GFF format

---------------------------


# What AGAT can do for you?  

## The Handler (gxf_to_gff3.pl)
  => Standardized any GTF/GFF file into a proper GFF3 format
  =>
  =>

## Installation

    * Clone AGAT

    ```
    git clone https://github.com/NBISweden/AGAT.git
    cd AGAT
    ```

    *  Check all the dependencies

    ```
    perl Makefile.PL     
    ```
    => If dependencies are missing you can install them using cpan/cpanm or use the conda environment agat_environment.yml

    * Compile

    ```
    make
    ```

    * Test

    ```
    make test
    ```

    * Install

    ```
    make install
    ```

* On MS Windows, instead of make you'd probably have to use dmake or nmake depending the toolchain you have.

## UnInstall

```
perl uninstall_AGAT
```


## Scripts


### with) \_sp\_ prefix => Means SLURP

The gff file will be charged in memory in a way to facilitate access to desired features at any time. It has a memory cost but make life smoother. Indeed, it allows to perform complicated tasks and more time efficient. Moreover, it allows to fix all potential errors in the limit of the possibilities given by the format itself. See SLURP section for more information about it.  

### with) \_sq\_ prefix => Means SEQUENTIAL

The gff file is read and processed from its top to the end line by line without sanity check. This is memory efficient.

### SLURP - GFF3 Standardization for a full GFF3 compliant to any tool  

#### Internal data structure:  

The method create a hash structure containing all the data in memory. We call it OMNISCIENT. The OMNISCIENT structure is a three levels structure :

$omniscient{level1}{tag_l1}{level1_id} = feature <= tag could be gene, match  
$omniscient{level2}{tag_l2}{idY} = @featureListL2 <= tag could be mRNA,rRNA,tRNA,etc. idY is a level1_id (know as Parent attribute within the level2 feature). The @featureList is a list to be able to manage isoform cases.  
$omniscient{level3}{tag_l3}{idZ} =  @featureListL3 <= tag could be exon,cds,utr3,utr5,etc. idZ is the ID of a level2 feature (know as Parent attribute within the level3 feature). The @featureList is a list to be able to put all the feature of a same tag together.  

#### How the slurp method works:  

The slurp parser phylosophy:
  * 1) Parse by Parent/child relationship
  * 2) ELSE Parse by a common tag  (an attribute value shared by feature that must be grouped together. By default we are using locus_tag but can be set by parameter)
  * 3) ELSE Parse sequentially (mean group features in a bucket, and the bucket change at each level2 feature, and bucket are join in a common tag at each new L1 feature)

/!\ Case with only level3 features (i.e rast or some prokka files, sequential will not work as expected. Indeed all features will be the child of only one newly created Parent.
    To create a parent per feature or group of feature, a common tag must be used to group them correctly. We use )

To resume by priority of way to parse: Parent/child relationship > locus_tag > sequential

#### What the slurp method do for you:

=> It creates missing parental features. (e.g if a level2 or level3 feature do not have parental feature(s) we create the missing level2 and/or level1 feature(s))  
=> It creates missing mandatory attributes (ID and/or Parent)
=> It fixes identifier to be uniq  
=> It removes duplicated features (same position, same ID, same Parent).  
=> It expands level3 features sharing multiple parents (e.g  if one exon has list of multiple parent mRNA in its Parent attribute, one exon per parent with uniq ID will be created.  
=> It fixes feature location errors (e.g an mRNA spanning over its gene location, we fix the gene location).
=> It adds UTR if possible (CDS and exon present)
=> It add exon if possible (CDS has to be present)
=> It group features together (if related features are spread at different place in the file)



#### examples:
