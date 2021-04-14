# AGAT - **A**nother **G**tf/gff **A**nalysis **T**oolkit
## Suite of tools to handle gene annotations in any GTF/GFF format.
---------------------------------------------

# Table of Contents

* [Foreword](#foreword)
* [List of AGAT tools (v0.6.0)](#list-of-agat-tools-v060)
* [Topological sorting of gff features](https://github.com/NBISweden/AGAT/wiki/Topological-sorting-of-gff-features)

## Foreword
Providing support in genome annotation within [NBIS](https://nbis.se) the GTF/GFF format is the main format I handle. I receive from customers file in GTF/GFF format coming from a broad range of sources. Even sometimes files from mixed sources (concatenated in the same file), or manually edited.  
The problem is that often those files do not follow the official specifications or even if they do, they are not even be sure to be compatible we the inputs expected by the tools.  

* The main idea was **first** to be able to **parse all possible cases** that can be met (I listed more than 30 cases). To my knowledge AGAT is the only one able to handle all of them.

* The **second** idea was to be able to **create a full standardised GFF3** file that could actually fit in any tool.
Once again AGAT is the only one recreating fully the missing information:
   * missing features (gene, mRNA, tRNA, exon, UTRs, etc...)
   * missing attributes (ID, Parent).

   and fixing wrong information:
   * identifier to be uniq.
   * feature location (e.g mRNA will be stretched if shorter than its exons).
   * remove duplicated features.
   * merge overlapping loci (if option activate because for prokaryote is not something we would like)

* The **third** idea was to have a **correct topological sorting output**. To my knowledge AGAT is the only one dealing properly with this task. More information about it [here](https://github.com/NBISweden/AGAT/wiki/Topological-sorting-of-gff-features).

* **Finally**, based on the abilities described previously I have developed a **toolkit to perform different tasks**. Some are originals, some are similar than what other tools could offer, but within AGAT they will always have the strength of the 3 first points.


**A final word**  
AGAT can solve lot of complicated cases and save headaches.  
Enjoy!!

## List of AGAT tools (v0.6.1)
[agat_convert_bed2gff.pl](https://github.com/NBISweden/AGAT/wiki/agat_convert_bed2gff)  
[agat_convert_embl2gff.pl](https://github.com/NBISweden/AGAT/wiki/agat_convert_embl2gff)  
[agat_convert_genscan2gff.pl](https://github.com/NBISweden/AGAT/wiki/agat_convert_genscan2gff)  
[agat_convert_mfannot2gff.pl](https://github.com/NBISweden/AGAT/wiki/agat_convert_mfannot2gff)  
[agat_convert_minimap2_bam2gff.pl](https://github.com/NBISweden/AGAT/wiki/agat_convert_minimap2_bam2gff)  
[agat_convert_sp_gff2bed.pl](https://github.com/NBISweden/AGAT/wiki/agat_convert_sp_gff2bed)  
[agat_convert_sp_gff2gtf.pl](https://github.com/NBISweden/AGAT/wiki/agat_convert_sp_gff2gtf)  
[agat_convert_sp_gff2tsv.pl](https://github.com/NBISweden/AGAT/wiki/agat_convert_sp_gff2tsv)  
[agat_convert_sp_gff2zff.pl](https://github.com/NBISweden/AGAT/wiki/agat_convert_sp_gff2zff)  
[agat_convert_sp_gxf2gxf.pl](https://github.com/NBISweden/AGAT/wiki/agat_convert_sp_gxf2gxf)  
[agat_sp_Prokka_inferNameFromAttributes.pl](https://github.com/NBISweden/AGAT/wiki/agat_sp_Prokka_inferNameFromAttributes)  
[agat_sp_add_introns.pl](https://github.com/NBISweden/AGAT/wiki/agat_sp_add_introns)  
[agat_sp_add_start_and_stop.pl](https://github.com/NBISweden/AGAT/wiki/agat_sp_add_start_and_stop)  
[agat_sp_alignment_output_style.pl](https://github.com/NBISweden/AGAT/wiki/agat_sp_alignment_output_style)  
[agat_sp_clipN_seqExtremities_and_fixCoordinates.pl](https://github.com/NBISweden/AGAT/wiki/agat_sp_clipN_seqExtremities_and_fixCoordinates)  
[agat_sp_compare_two_BUSCOs.pl](https://github.com/NBISweden/AGAT/wiki/agat_sp_compare_two_BUSCOs)  
[agat_sp_compare_two_annotations.pl](https://github.com/NBISweden/AGAT/wiki/agat_sp_compare_two_annotations)  
[agat_sp_complement_annotations.pl](https://github.com/NBISweden/AGAT/wiki/agat_sp_complement_annotations)  
[agat_sp_ensembl_output_style.pl](https://github.com/NBISweden/AGAT/wiki/agat_sp_ensembl_output_style)  
[agat_sp_extract_attributes.pl](https://github.com/NBISweden/AGAT/wiki/agat_sp_extract_attributes)  
[agat_sp_extract_sequences.pl](https://github.com/NBISweden/AGAT/wiki/agat_sp_extract_sequences)  
[agat_sp_filter_by_ORF_size.pl](https://github.com/NBISweden/AGAT/wiki/agat_sp_filter_by_ORF_size)  
[agat_sp_filter_by_locus_distance.pl](https://github.com/NBISweden/AGAT/wiki/agat_sp_filter_by_locus_distance)  
[agat_sp_filter_by_mrnaBlastValue.pl](https://github.com/NBISweden/AGAT/wiki/agat_sp_filter_by_mrnaBlastValue)  
[agat_sp_filter_feature_by_attribute_presence.pl](https://github.com/NBISweden/AGAT/wiki/agat_sp_filter_feature_by_attribute_presence)  
[agat_sp_filter_feature_by_attribute_value.pl](https://github.com/NBISweden/AGAT/wiki/agat_sp_filter_feature_by_attribute_value)  
[agat_sp_filter_feature_from_keep_list.pl](https://github.com/NBISweden/AGAT/wiki/agat_sp_filter_feature_from_keep_list)  
[agat_sp_filter_feature_from_kill_list.pl](https://github.com/NBISweden/AGAT/wiki/agat_sp_filter_feature_from_kill_list)  
[agat_sp_filter_gene_by_intron_numbers.pl](https://github.com/NBISweden/AGAT/wiki/agat_sp_filter_gene_by_intron_numbers)  
[agat_sp_filter_gene_by_length.pl](https://github.com/NBISweden/AGAT/wiki/agat_sp_filter_gene_by_length)  
[agat_sp_filter_incomplete_gene_coding_models.pl](https://github.com/NBISweden/AGAT/wiki/agat_sp_filter_incomplete_gene_coding_models)  
[agat_sp_filter_record_by_coordinates.pl](https://github.com/NBISweden/AGAT/wiki/agat_sp_filter_record_by_coordinates)  
[agat_sp_fix_cds_phases.pl](https://github.com/NBISweden/AGAT/wiki/agat_sp_fix_cds_phases)  
[agat_sp_fix_features_locations_duplicated.pl](https://github.com/NBISweden/AGAT/wiki/agat_sp_fix_features_locations_duplicated)  
[agat_sp_fix_fusion.pl](https://github.com/NBISweden/AGAT/wiki/agat_sp_fix_fusion)  
[agat_sp_fix_longest_ORF.pl](https://github.com/NBISweden/AGAT/wiki/agat_sp_fix_longest_ORF)  
[agat_sp_fix_overlaping_genes.pl](https://github.com/NBISweden/AGAT/wiki/agat_sp_fix_overlaping_genes)  
[agat_sp_fix_small_exon_from_extremities.pl](https://github.com/NBISweden/AGAT/wiki/agat_sp_fix_small_exon_from_extremities)  
[agat_sp_flag_premature_stop_codons.pl](https://github.com/NBISweden/AGAT/wiki/agat_sp_flag_premature_stop_codons)  
[agat_sp_flag_short_introns.pl](https://github.com/NBISweden/AGAT/wiki/agat_sp_flag_short_introns)  
[agat_sp_functional_statistics.pl](https://github.com/NBISweden/AGAT/wiki/agat_sp_functional_statistics)  
[agat_sp_keep_longest_isoform.pl](https://github.com/NBISweden/AGAT/wiki/agat_sp_keep_longest_isoform)  
[agat_sp_kraken_assess_liftover.pl](https://github.com/NBISweden/AGAT/wiki/agat_sp_kraken_assess_liftover)  
[agat_sp_list_short_introns.pl](https://github.com/NBISweden/AGAT/wiki/agat_sp_list_short_introns)  
[agat_sp_load_function_from_protein_align.pl](https://github.com/NBISweden/AGAT/wiki/agat_sp_load_function_from_protein_align)  
[agat_sp_manage_IDs.pl](https://github.com/NBISweden/AGAT/wiki/agat_sp_manage_IDs)  
[agat_sp_manage_UTRs.pl](https://github.com/NBISweden/AGAT/wiki/agat_sp_manage_UTRs)  
[agat_sp_manage_attributes.pl](https://github.com/NBISweden/AGAT/wiki/agat_sp_manage_attributes)  
[agat_sp_manage_functional_annotation.pl](https://github.com/NBISweden/AGAT/wiki/agat_sp_manage_functional_annotation)  
[agat_sp_manage_introns.pl](https://github.com/NBISweden/AGAT/wiki/agat_sp_manage_introns)  
[agat_sp_merge_annotations.pl](https://github.com/NBISweden/AGAT/wiki/agat_sp_merge_annotations)  
[agat_sp_prokka_fix_fragmented_gene_annotations.pl](https://github.com/NBISweden/AGAT/wiki/agat_sp_prokka_fix_fragmented_gene_annotations)  
[agat_sp_sensitivity_specificity.pl](https://github.com/NBISweden/AGAT/wiki/agat_sp_sensitivity_specificity)  
[agat_sp_separate_by_record_type.pl](https://github.com/NBISweden/AGAT/wiki/agat_sp_separate_by_record_type)  
[agat_sp_statistics.pl](https://github.com/NBISweden/AGAT/wiki/agat_sp_statistics)  
[agat_sp_webApollo_compliant.pl](https://github.com/NBISweden/AGAT/wiki/agat_sp_webApollo_compliant)  
[agat_sq_add_attributes_from_tsv.pl](https://github.com/NBISweden/AGAT/wiki/agat_sq_add_attributes_from_tsv)  
[agat_sq_add_hash_tag.pl](https://github.com/NBISweden/AGAT/wiki/agat_sq_add_hash_tag)  
[agat_sq_add_locus_tag.pl](https://github.com/NBISweden/AGAT/wiki/agat_sq_add_locus_tag)  
[agat_sq_keep_annotation_from_fastaSeq.pl](https://github.com/NBISweden/AGAT/wiki/agat_sq_keep_annotation_from_fastaSeq)  
[agat_sq_list_attributes.pl](https://github.com/NBISweden/AGAT/wiki/agat_sq_list_attributes)  
[agat_sq_manage_IDs.pl](https://github.com/NBISweden/AGAT/wiki/agat_sq_manage_IDs)  
[agat_sq_manage_attributes.pl](https://github.com/NBISweden/AGAT/wiki/agat_sq_manage_attributes)  
[agat_sq_mask.pl](https://github.com/NBISweden/AGAT/wiki/agat_sq_mask)  
[agat_sq_remove_redundant_entries.pl](https://github.com/NBISweden/AGAT/wiki/agat_sq_remove_redundant_entries)  
[agat_sq_repeats_analyzer.pl](https://github.com/NBISweden/AGAT/wiki/agat_sq_repeats_analyzer)  
[agat_sq_rfam_analyzer.pl](https://github.com/NBISweden/AGAT/wiki/agat_sq_rfam_analyzer)  
[agat_sq_split.pl](https://github.com/NBISweden/AGAT/wiki/agat_sq_split)  
[agat_sq_stat_basic.pl](https://github.com/NBISweden/AGAT/wiki/agat_sq_stat_basic)  
