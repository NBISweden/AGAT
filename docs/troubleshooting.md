# Troubleshooting

## AGAT throws features out, because the feature type is not yet taken into account
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

## AGAT throws features out, because child features are not provided
Features level1 (e.g. gene, match, chromosome) may require to have child features or not depending of the information stored into the `features_level1.json` file. If a child is required, and the GFF file does not contain it, the level1 feature will be thrown away. You must modify the json file to add the the term `standalone` to inform AGAT that this feature level1 do not require any child. (This work only on feature level1, not level2 or level3). To access the json files run the following command:
```
# export the json files
agat_convert_sp_gxf2gxf.pl --expose
```
Then open the `features_level1.json` and put the value `standalone` as value to the required feature.
Finally run your scripts in the same folder as the modified json files are standing.

### How to use a version of AGAT from a specific branch
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

### How to fix Bio::Root::Exception met in AGAT

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