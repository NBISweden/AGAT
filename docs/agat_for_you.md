# What can AGAT do for you?

*AGAT a GFF/GTF toolkit allowing you to perform almost everything you might want to achieve ^^*

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
