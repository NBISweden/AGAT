# agat\_sp\_fix\_cds\_phases.pl

## DESCRIPTION

This script aims to fix the CDS phases.
The script is compatible with incomplete gene models (Missing start, CDS
multiple of 3 or not, i.e. with offset of 1 or 2) and + and - strand.

How this script works?  
AGAT uses the fasta sequence to verify the CDS frame.
In case the CDS start by a start codon the phase of the first CDS piece is set to 0.
In the case there is no start codon and: 
  - If there is only one stop codon in the sequence and it is located at the last position, the phase of the first CDS piece is set to 0.
  - If there is no stop codon, the phase of the first CDS piece is set to 0 (because sequence can be translated without premature stop codon).
  - If there is/are stop codon(s) in the middle of the sequence we re-execute the check with an offset of +2 nucleotides:
    - If there is only one stop codon in the sequence and it is located at the last position, the phase of the first CDS piece is set to 0.
    - If there is no stop codon, the phase of the first CDS piece is set to 0 (because sequence can be translated without premature stop codon).
    - If there is/are stop codon(s) in the middle of the sequence we re-execute the check with an offset of +1 nucleotide:
        - If there is only one stop codon in the sequence and it is located at the last position, the phase of the first CDS piece is set to 0.
        - If there is no stop codon, the phase of the first CDS piece is set to 0 (because sequence can be translated without premature stop codon).
        - If there is/are still stop codon(s) we keep original phase and throw a warning. In this last case it means we never succeded to make a translation without premature stop codon in all the 3 possible phases.
Then in case of CDS made of multiple CDS pieces (i.e. discontinuous feature), the rest of the CDS pieces will be checked accordingly to the first CDS piece.

What is a phase?  
For features of type "CDS", the phase indicates where the next codon begins
relative to the 5' end (where the 5' end of the CDS is relative to the strand
of the CDS feature) of the current CDS feature. For clarification the 5' end
for CDS features on the plus strand is the feature's start and and the 5' end
for CDS features on the minus strand is the feature's end. The phase is one of
the integers 0, 1, or 2, indicating the number of bases forward from the start
of the current CDS feature the next codon begins. A phase of "0" indicates that
a codon begins on the first nucleotide of the CDS feature (i.e. 0 bases forward),
a phase of "1" indicates that the codon begins at the second nucleotide of this
CDS feature and a phase of "2" indicates that the codon begins at the third
nucleotide of this region. Note that ‘Phase’ in the context of a GFF3 CDS
feature should not be confused with the similar concept of frame that is also a
common concept in bioinformatics. Frame is generally calculated as a value for
a given base relative to the start of the complete open reading frame (ORF) or
the codon (e.g. modulo 3) while CDS phase describes the start of the next codon
relative to a given CDS feature.  
The phase is REQUIRED for all CDS features.

## SYNOPSIS

```
agat_sp_fix_cds_phases.pl --gff infile.gff -f fasta [ -o outfile ]
agat_sp_fix_cds_phases.pl --help
```

## OPTIONS

- **-g**, **--gff** or **-ref**

    Input GTF/GFF file.

- **-fa** or **--fasta**

    Input fasta file.

- **-v** or **--verbose**

    Add verbosity.

- **-o** or **--output**

    Output GFF file. If no output file is specified, the output will be
    written to STDOUT.

- **-c** or **--config**

    String - Input agat config file. By default AGAT takes as input agat_config.yaml file from the working directory if any,
    otherwise it takes the orignal agat_config.yaml shipped with AGAT. To get the agat_config.yaml locally type: "agat config --expose".
    The --config option gives you the possibility to use your own AGAT config file (located elsewhere or named differently).

- **-h** or **--help**

    Display this helpful text.

