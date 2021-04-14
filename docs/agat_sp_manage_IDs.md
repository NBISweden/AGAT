# NAME

agat\_sp\_manage\_IDs.pl

# DESCRIPTION

The script takes a gff3 file as input and will go through all feature to overwrite
the value of the ID attribute.
By default the ID is built as follow: primary\_tag(i.e. 3rd column)-Number.
If you provide a specific prefix the ID is built as follow: $opt\_prefix.$letterCode.Number.
By default the numbering start at 1, but you can decide to change this value using the --nb option.
The $letterCode is the first letter of the feature type (3rd colum). It is uniq for each feature type,
i.e. when two feature types start with the same letter, the second one met will have the two first letter as $letterCode (and so one).

# SYNOPSIS

```
agat_sp_manage_IDs.pl --gff file.gff -p level2 -p cds -p exon [ -o outfile ]
agat_sp_manage_IDs.pl --help
```

# OPTIONS

- **--gff** or **-f**

    Input GTF/GFF file.

- **--gap**

    Integer - Increment the next gene (level1 feature) suffix with this value. Defauft 0.

- **--ensembl**

    Boolean - For an ID Ensembl like (e.g PREFIXG00000000022). The ID is built as follow:
    $opt\_prefix.$letterCode.0\*.Number where the number of 0 is adapted in order to have 11 digits.

- **--prefix**

    String - Add a specific prefix to the ID. By defaut if will be the feature type (3rd column).

- **--type\_dependent**

    Boolean - Activate type\_dependent numbering. The number is depedendent of the feature type.
    i.e instead of:
    NbV1Ch01        AUGUSTUS        gene    97932   99714   0.06    -       .       ID=gene1
    NbV1Ch01        AUGUSTUS        mRNA    97932   99714   0.06    -       .       ID=mRNA2
    NbV1Ch01        AUGUSTUS        exon    97932   98571   .       -       .       ID=exon3
    NbV1Ch01        AUGUSTUS        exon    98679   98844   .       -       .       ID=exon4
    You will get:
    NbV1Ch01        AUGUSTUS        gene    97932   99714   0.06    -       .       ID=gene1
    NbV1Ch01        AUGUSTUS        mRNA    97932   99714   0.06    -       .       ID=mRNA1
    NbV1Ch01        AUGUSTUS        exon    97932   98571   .       -       .       ID=exon1
    NbV1Ch01        AUGUSTUS        exon    98679   98844   .       -       .       ID=exon2

- **--collective**

    Boolean - In the case of discontinuous features (i.e. a single feature that
    exists over multiple genomic locations like CDS, UTR) we set a uniq ID by default.
    If you wish to set the a collective ID for those feature, please activate this option.

- **--tair**

    Boolean - Tair like Output:

    NbV1Ch01    TAIR10  gene    5928    8737    .       -       .       ID=AT1G01020
    NbV1Ch01    TAIR10  mRNA    5928    8737    .       -       .       ID=AT1G01020.1
    NbV1Ch01    TAIR10  exon    5928    8737   .       -       .        ID=AT1G01020.1-exon1

- **--nb**

    Integer - Start numbering to this value. Default 1.

- **-p**,  **-t** or  **-l**

    primary tag option, case insensitive, list. Allow to specied the feature types that will be handled.
    You can specified a specific feature by given its primary tag name (column 3) as: cds, Gene, MrNa
    You can specify directly all the feature of a particular level:
          level2=mRNA,ncRNA,tRNA,etc
          level3=CDS,exon,UTR,etc
    By default all feature are taken into account. fill the option by the value "all" will have the same behaviour.

- **-o** , **--output** , **--out** or **--outfile**

    String - Output GFF file. If no output file is specified, the output will be
    written to STDOUT.

- **-h** or **--help**

    Boolean - Display this helpful text.

