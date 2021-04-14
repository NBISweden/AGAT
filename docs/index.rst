Welcome to AGAT's documentation!
================================
<h2><em>A</em>nother <em>G</em>tf/Gff <em>A</em>nalysis <i>T</i>oolkit</h2>
Suite of tools to handle gene annotations in any GTF/GFF format.

[<img src="workcloud.png" width="600" height="300" />]

The GTF/GFF format
==================

The GTF/GFF formats are 9-column text formats used to describe and represent genomic features.
The formats have quite evolved since 1997, and despite well-defined specifications existing nowadays they have a great flexibility allowing holding wide variety of information.
This flexibility has a drawback aspect, there is an incredible amount of flavour of the formats, that can result in problems when using downstream programs.
For a complete overview of the GTF/GFF formats have a look [here](https://github.com/NBISweden/GAAS/blob/master/annotation/knowledge/gxf.md).

## What can AGAT do for you?

AGAT has the power to check, fix, pad missing information (features/attributes) of any kind of GTF and GFF to create complete, sorted and standardised gff3 format.
The toolkit comes with an exhaustive list of tools allowing to perform almost everything you might want to achieve ^^

Why this tool?
=============

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

Contents
========

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   topological-sorting-of-gff-features.md
   troubleshooting.md

   agat_convert_bed2gff.md
   agat_convert_embl2gff.md
   agat_convert_genscan2gff.md
   agat_convert_mfannot2gff.md
   agat_convert_minimap2_bam2gff.md
   agat_convert_sp_gff2bed.md
   agat_convert_sp_gff2gtf.md
   agat_convert_sp_gff2tsv.md
   agat_convert_sp_gff2zff.md
   agat_convert_sp_gxf2gxf.md
   agat_sp_Prokka_inferNameFromAttributes.md
   agat_sp_add_introns.md
   agat_sp_add_start_and_stop.md
   agat_sp_alignment_output_style.md
   agat_sp_clipN_seqExtremities_and_fixCoordinates.md
   agat_sp_compare_two_BUSCOs.md
   agat_sp_compare_two_annotations.md
   agat_sp_complement_annotations.md
   agat_sp_ensembl_output_style.md
   agat_sp_extract_attributes.md
   agat_sp_extract_sequences.md
   agat_sp_filter_by_ORF_size.md
   agat_sp_filter_by_locus_distance.md
   agat_sp_filter_by_mrnaBlastValue.md
   agat_sp_filter_feature_by_attribute_presence.md
   agat_sp_filter_feature_by_attribute_value.md
   agat_sp_filter_feature_from_keep_list.md
   agat_sp_filter_feature_from_kill_list.md
   agat_sp_filter_gene_by_intron_numbers.md
   agat_sp_filter_gene_by_length.md
   agat_sp_filter_incomplete_gene_coding_models.md
   agat_sp_filter_record_by_coordinates.md
   agat_sp_fix_cds_phases.md
   agat_sp_fix_features_locations_duplicated.md
   agat_sp_fix_fusion.md
   agat_sp_fix_longest_ORF.md
   agat_sp_fix_overlaping_genes.md
   agat_sp_fix_small_exon_from_extremities.md
   agat_sp_flag_premature_stop_codons.md
   agat_sp_flag_short_introns.md
   agat_sp_functional_statistics.md
   agat_sp_gxf_to_gff3.md
   agat_sp_keep_longest_isoform.md
   agat_sp_kraken_assess_liftover.md
   agat_sp_list_short_introns.md
   agat_sp_load_function_from_protein_align.md
   agat_sp_manage_IDs.md
   agat_sp_manage_UTRs.md
   agat_sp_manage_attributes.md
   agat_sp_manage_functional_annotation.md
   agat_sp_manage_introns.md
   agat_sp_merge_annotations.md
   agat_sp_prokka_fix_fragmented_gene_annotations.md
   agat_sp_sensitivity_specificity.md
   agat_sp_separate_by_record_type.md
   agat_sp_split_by_level2_feature.md
   agat_sp_statistics.md
   agat_sp_to_tabulated.md
   agat_sp_webApollo_compliant.md
   agat_sq_add_attributes_from_tsv.md
   agat_sq_add_hash_tag.md
   agat_sq_add_locus_tag.md
   agat_sq_keep_annotation_from_fastaSeq.md
   agat_sq_list_attributes.md
   agat_sq_manage_ID.md
   agat_sq_manage_IDs.md
   agat_sq_manage_attributes.md
   agat_sq_mask.md
   agat_sq_remove_redundant_entries.md
   agat_sq_repeats_analyzer.md
   agat_sq_rfam_analyzer.md
   agat_sq_split.md
   agat_sq_stat_basic.md


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
