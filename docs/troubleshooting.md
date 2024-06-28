# Troubleshooting

## AGAT throws features out, because the feature type is not yet taken into account
Feature types (primary_tag) handled by AGAT are defined within `feature_levels.yaml` file. Most common features are already defined in this file. If you encounter GFF/GTF files with feature types not accepted, AGAT will inform you and throw the features out. To keep those feature you must inform AGAT how to handle them:
First access the feature_levels.yaml files by running:
```
# expose the yaml files
agat levels --expose
```

Then open the file with your favorite text editor.

Now choose which section you want to modify:
* `level1` (e.g. gene, match, region):
  For features that do not have parent
* `level2` (e.g. mrna, match_part, trna):
  For features that have one parent and the parent is a level 1 feature.
* `level3` (e.g. exon, intron, cds):
  For features that have one parent (the parent has also a parent) and no children

For features that are discontinuous (i.e. when a single feature exists over multiple genomic locations like cds, utr) you must also fil the `spread` section.

Then add the feature type information by adding a paired-value like this:
```
	"bigRNA":"gene",
```
Where `bigRNA`is the feature type and `gene` the parent feature type expected.
/!\\ For level1 feature type the second value can be:
 * topfeature: feature does not expect children, and will be written first in the sequence
 * standalone: feature does not expect children
 * other values do not have any meaning but a value is required, write whatever you want.

## AGAT throws features out, because child features are not provided
Features level1 (e.g. gene, match, chromosome) may require to have child features or not depending of the information stored into the `feature_levels.yaml` file. If a child is required, and the GFF file does not contain it, the level1 feature will be thrown away. You must modify the `feature_levels.yaml` file to add the the term `standalone` to inform AGAT that this feature level1 do not require any child. (This work only on feature level1, not level2 or level3). To access the `feature_levels.yaml` file run the following command:
```
# expose the yaml files
agat levels --expose
```
Then open the `feature_levels.yaml` and put the value `standalone` as value to the required feature.
Finally run your scripts in the same folder as the modified `feature_levels.yaml` file is standing.

## How to use a version of AGAT from a specific branch
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

## How to fix Bio::Root::Exception met in AGAT

They are bioperl error messages. Encountered error can be:  
```
MSG: Failed validation of sequence '[unidentified sequence]'. Invalid characters 
```
or  
```
MSG: Each line of the file must be less than 65,536 characters.
```

Bio::DB::Fasta from Bioperl cannot handle line with more than 65,536 characters. So you must fold you fasta sequence before to run AGAT's scripts:
```
# Fold to 80 characters by line. 
# Be careful if you have long headers that can be folded over several lines. You must first shorten them, or fold with higher value.
fold input.fa > output.fa
```


## How to use codon table 0 (codon table 1 is used instead)?

Several scripts need to use a codon table:

```
agat_sp_add_start_and_stop.pl  
agat_sp_extract_sequences.pl  
agat_sp_filter_incomplete_gene_coding_models.pl  
agat_sp_fix_fusion.pl
agat_sp_fix_longest_ORF.pl
agat_sp_fix_small_exon_from_extremities.pl
agat_sp_flag_premature_stop_codons.pl
agat_sp_prokka_fix_fragmented_gene_annotations.pl
```

By default AGAT uses codon table 1 wich is the standard table. 

* What is the difference between table 1 and table 0?  
  The codon table 0 is strict and uses ATG-only start codon, while codon table 1 uses TTG, CTG and ATG possible start codon.

* What are the possible codon table?   
  In top of the table 0 which is specific to Bioperl many other tables are available. Their description can be found [here](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG1).

## Why when asking for table 0 AGAT keep using the table 1 ?

There are two possible reasons for that problem.

  * AGAT: Originally the problem comes from a bug in Bioperl. AGAT was trying to pass by the problem but the fix was not efficient until version 1.4.1. Please be sure to use a version >= 1.4.1 to avoid any problem from the AGAT side.

  * Bioperl: The problem has been present for a while and has been definitly fixed in the [commit fa9366f from the 24th of April 2024](https://github.com/bioperl/bioperl-live/tree/fa9366f3a2f48fd051343d488cfce26655f842b3).  
  So to fix the problem you need to use a bioperl version equal or later to that point. If not possible (e.g. not yet available for installation via conda) you can follow this procedure:
      * run AGAT once
      * catch the location of bioperl used from the prompt e.g.:
        ```
         ------------------------------------------------------------------------------
        |   Another GFF Analysis Toolkit (AGAT) - Version: v1.4.0                      |
        |   https://github.com/NBISweden/AGAT                                          |
        |   National Bioinformatics Infrastructure Sweden (NBIS) - www.nbis.se         |
        ------------------------------------------------------------------------------
        
        ...

        => Machine information:
              This script is being run by perl v5.32.1
              Bioperl location being used: /usr/local/lib/perl5/site_perl/Bio/
              Operating system being used: linux
        ```
        Here the bioperl path is here: `/usr/local/lib/perl5/site_perl/Bio/`
      * Move into the directory found in the previous step minus `/Bio`:
        `cd /usr/local/lib/perl5/site_perl/`
      * Copy paste locally the file and the folder from the bioperl-live repository (here)[https://github.com/bioperl/bioperl-live/tree/master/lib]:  
        ```
        git clone https://github.com/bioperl/bioperl-live
        cp -r bioperl-live/lib/* .
        ```
        
      Now you should be able to use the codon table 0. If not check your AGAT version (see above).



