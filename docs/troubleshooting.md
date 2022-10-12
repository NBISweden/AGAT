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