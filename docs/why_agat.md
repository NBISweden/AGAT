Why AGAT?
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
