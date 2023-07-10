# agat\_sp\_complement\_annotations.pl

## DESCRIPTION

The script allows to complement a reference annotation with other annotations.
A l1 feature from the addfile.gff that does not overlap a l1 feature from the reference annotation will be added.
A l1 feature from the addfile.gff without a CDS that overlaps a l1 feature with a CDS from the reference annotation will be added.
A l1 feature from the addfile.gff with a CDS that overlaps a l1 feature without a CDS from the reference annotation will be added.
A l1 feature from the addfile.gff with a CDS that overlaps a l1 feature with a CDS from the reference annotation will be added only if the CDSs don't overlap.
A l1 feature from the addfile.gff without a CDS that overlaps a l1 feature without a CDS from the reference annotation will be added only if none of the l3 features overlap.
/!\\ It is sufficiant that only one isoform is overlapping to prevent the whole gene (l1 feature) from the addfile.gff to be added in the output.

## SYNOPSIS

```
agat_sp_complement_annotations.pl --ref annotation_ref.gff --add addfile1.gff --add addfile2.gff --out outFile
agat_sp_complement_annotations.pl --help
```

## OPTIONS

- **--ref**,  **-r** or **-i**

    Input GTF/GFF file used as reference.

- **--add** or **-a**

    Annotation(s) file you would like to use to complement the reference annotation. You can specify as much file you want like so: -a addfile1 -a addfile2 -a addfile3
    /!\\ The order you provide these files matter. Once the reference file has been complemented by file1, this new annotation becomes the new reference that will be complemented by file2 etc.
    /!\\ The result with -a addfile1 -a addfile2 will differ to the result from -a addfile2 -a addfile1. So, be aware of what you want if you use several addfiles.

- **--size\_min** or **-s**

    Option to keep the non-overlping gene only if the CDS size (in nucleotide) is over the minimum size defined. Default = 0 that means all of them are kept.

- **--out**, **--output**, **--outfile** or **-o**

    Output gff3 containing the reference annotation with all the non-overlapping newly added genes from addfiles.gff.

- **-c** or **--config**

    String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any,
    otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
    The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

- **--help** or **-h**

    Display this helpful text.

