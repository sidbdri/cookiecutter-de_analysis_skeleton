# Basic RNA-seq analysis results structure

At the top-level, the main results folder should contain the following files and sub-folders:

- `*multiqc_report.html` – One or more web pages containing sequence-level and analysis pipeline quality control data.
- `differential_expression` – Contains gene-level differential expression results, graphs, and downstream analyses.
- `differential_expression_tx` – Contains transcript-level differential expression results and graphs.
- (plus `sessionInfo.txt` – reproducibility information, for our records).

### Differential gene expression

We use the "DESeq2" tool for differential expression (_Love et al._, "[Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8)", Genome Biology (2014)). Differential gene expression results sheets are in the sub-folder `differential_expression/de_gene`:

- `*deseq2_results_fpkm*.csv` - One or more CSV files (these can be opened directly in Excel) containing results of differential gene expression analysis.
- `*de_summary*.csv` - One or more CSV files containing summaries of the numbers of genes detected as differentially expressed (at a particular false discovery rate - unless otherwise requested, this will normally be 5%).

Each of the one or more main differential gene expression results files (there may be more than one if there are results for multiple species, or if a single file was too unwieldy) contains the following columns:

* **gene**: Ensembl gene ID.
* **gene_name**: Ensemble gene name.
* **chromosome**: Chromosome on which the gene lies.
* **description**: Ensembl gene description.
* **entrez_id**: NCBI Gene ID (this is mainly for our internal use).
* **gene_type**: e.g. protein coding, pseudogene, LINC RNA etc. etc.
* **gene_length**: number of bases contained in the union of all transcripts for the gene.
* **max\_transcript\_length**: maximum length of any one isoform of the gene - this is used to calculate FPKMs.
* **\*_fpkm**: Per-sample gene abundances measured in **f**ragments **p**er **k**ilobase per **m**illion mapped reads (n.b. if Sargasso has been used, there may be separate pre- and post-Sargasso FPKM columns for single-species samples).
* **\*\_fpkm\_avg**: Per-condition FPKMs averaged over the samples in that condition.

and for each differential expression comparison:

* **\<comparison\_name\>_l2fc**: DESeq2-moderated log2 fold change in gene expression between conditions.
* **\<comparison\_name\>_pval**: raw p-value for statistical significance of differential expression.
* **\<comparison\_name\>_padj**: adjusted p-value (i.e. false discovery rate) after correcting for multiple testing.
* **\<comparison\_name\>\_raw\_l2fc**: raw log2 fold change in gene expression between conditions, calculated directly from normalised read counts (this may be different to the fold change calculated by DESeq2, which does some further processing beyond the raw value).

_n.b._:

* In the first instance, the **l2fc** and **padj** columns are the ones to look at.
* In a comparison named "B\_vs\_A", a positive log2 fold change means that gene expression is higher in the "comparison" condition B than in the "base" condition A, and vice-versa for a negative fold change. If there is any confusion about which are the "comparison" and "base" conditions (e.g. if we have named the comparison incorrectly) then consulting the appropriate entry in the `*de_summary*.csv` file will indicate exactly which conditions were used as comparison and base.
* p-values and fold changes may be missing for some genes. If the raw p-value is missing, this means that the expression of the gene across all samples is negligible, and the gene was filtered before differential expression was performed (in this case, the fold change and adjusted p-value columns will be empty too). If the raw p-value is present, but the adjusted p-value is missing, this means that the gene was filtered by DESeq2 prior to multiple testing correction (DESeq2 tries to reduce the burden of this correction by excluding genes from its final calculation that it believes, due to their expression being below a comparison-specific expression threshold, are unlikely to be discovered as significantly differentially expressed).
* The raw log2 fold change column may contain the values `Inf` or `-Inf` – these indicate that expression in one of the experimental conditions being compared was zero in all samples, and hence a finite raw fold change could not be calculated.
* When testing the interaction between two factors (if such tests exist, you'll have either asked for them, or we'll explain why we've added them!), only the raw log2 fold change column is included, since DESeq2 does not further process fold changes in this case.

### Graphs

The directory `differential_expression/graphs` contains a range of plots describing data quality control and differential expression results. While some of these are mainly for internal quality control, a number may be of general interest:

- For each main differential expression comparison made there are two plots describing the samples included in the comparison. (i) The PCA plot (`pca_<comparison_name>_<species>.pdf`) graphs samples along the two axes of gene expression by which samples differ by the greatest amount (thus if there are large differences in expression between conditions, you would expect samples from each condition to cluster together - however, a lack of such clustering does not necessarily mean that there will be no differentially expressed genes); accompanying each PCA plot is a CSV file containing sample co-ordinates, in case you wish to re-plot this data yourself. (ii) The heatmaps (`heatmap_<comparison_name>_<species>.pdf`) show inter-sample "distance" - where each sample is considered as a point in an n-dimensional space (n=number of genes, sample location in the space determine by read count for each gene), and the distance is then just the Euclidean distance between points in that space; again, when there are large differences in expression between conditions, you would expect clustering of samples from each condition. 
- There are also PCA (`pca_all_samples_<species>.pdf`) and heatmap (`heatmap_all_samples_<species>.pdf`) plots showing all samples. In addition an all-sample PCA plot is produced (`pca_all_features_all_samples_<mouse>.pdf`) which is faceted by the different meta-data that is known about the samples (e.g. this might be by genotype and treatment in the case where those are the factors of interest in the experiment). 
- For a number of different brain cell types, graphs are included which plot per-sample expression of a couple of cell type marker genes (`cell_type_check_<type>_<species>.pdf`), which can be useful to spot potential contamination if data is expected to derive from a certain cell type.
- For each differential expression comparison, a scatter plot is produced showing average gene expression in the two conditions (`scatter_fpkm_<comparison_name>_<species>.pdf`); significantly up-regulated genes are marked in blue, and significantly down-regulated in red.
- For each differential expression comparison, plots are produced showing the expression of the four most significantly differentially-expressed genes across all samples (`top_de_genes_fpkm_<comparison_name>_<species>.pdf`); these can be useful to spot unexpected expression patterns of (potentially unexpected) genes.
- Finally, for each differential expression comparison, a volcano plot (`volcano_plot_<comparison_name>_<species>.pdf`) is produced showing p-values and log2 fold changes for each gene. The five most significantly up-regulated and five most significantly down-regulated genes are labelled.

_n.b._ these plots are mainly used to spot gross QC problems and sample mix-ups, but frequently nothing conclusive can be drawn from them.

### Gene Ontology enrichment analysis

In the sub-folder `differential_expression/enrichment_tests`, there are Gene Ontology (GO) enrichment analyses for each differential expression comparison (these may be in further species-specific sub-folders). 

For each differential expression comparison we take (i) all genes called differentially expressed at false discovery rate < 0.05, (ii) just the up-regulated genes and (iii) just the down-regulated genes. Then for each of the GO categories of "biological process", "cellular compartment" and "molecular function", we use an algorithm ("topGO" – _Alexa et al._, ["Improved scoring of functional groups from gene expression data by decorrelating GO graph structure"](https://academic.oup.com/bioinformatics/article/22/13/1600/193669)), Bioinformatics (2006)) to test whether the differentially expressed genes are enriched in genes annotated with particular GO terms, as compared to the background set of all genes expressed in this data. These results can start to give some clue as to the particular biological processes that are being affected, without having to just scan through huge lists of gene expression data.

_n.b._ some of these GO analysis files might not be present, if there were no, or only a few, differentially expressed genes for a particular differential expression comparison.

Each GO enrichment results CSV file contains the following columns:

* **GO.ID**: Gene Ontology term ID.
* **Term**: Gene Ontology term description.
* **annotated\_in\_background**: Number of genes in the background set annotated with this term.
* **annotated\_in\_gene\_set**: Number of genes in the differentially expressed set annotated with this term.
* **expected\_annotated\_in\_gene\_set**: The number of genes you would expect to be annotated with this term in a random set of genes of the same sizes as the differentially expressed set.
* **p.value**: A p-value for enrichment of this Gene Ontology term in the differentially expressed set.
* **Genes**: The differentially expressed genes annotated with this Gene Ontology term. 

Note that the p-values here are _not_ corrected for multiple testing (it's not straightforward to do this, since the statistical tests for different Gene Ontology terms are not independent). A common convention, however, is to consider GO terms with p-value < 0.01 as "potentially interesting".

### Gene set enrichment analysis

While GO analyses are useful, they can suffer from our having imposed an arbitrary significance cutoff for the sets of genes considered; whereas sub-threshold, yet coherent, shifts in the expression of groups of genes between conditions can also be biologically meaningful. These coherent shifts in expression can be investigated using "gene set enrichment analysis"; results are found in the sub-folder `differential_expression/gene_set_tests` (and may be in further species-specific sub-folders).

The analysis method used here is called "Camera" (_Wu & Smyth_, ["Camera: a competitive gene set test accounting for inter-gene correlation"](https://academic.oup.com/nar/article/40/17/e133/2411151), Nucleic Acids Research (2012)). This implements a "competitive gene set test" - that is, for a particular set of genes of interest, it takes the expression values of all genes and tests whether the fold change in expression between experimental conditions for the genes in the set is different (as a whole) to the genes not in the set. (It also takes into account that fold changes of different genes are not necessarily independent - e.g. in cases where a bunch of genes in a pathway are all coherently differentially regulated, so the p-value obtained should be more robust than some other competitive gene set test methods that are out there.)

The main gene sets we use are divided into three categories, and consist of gene sets that have been compiled by the Broad Institute:

* **CURATED**: gene sets curated from various sources such as online pathway databases, the biomedical literature, and knowledge of domain experts. 
* **MOTIF**: gene sets representing potential targets of regulation by transcription factors or microRNAs.
* **GO**: gene sets that contain genes annotated by the same GO term.
* **MSIGDB\_CELL\_TYPE**: gene sets that contain curated cluster markers for cell types identified in single-cell sequencing studies of human tissue.

For historical reasons (i.e. the **MSIGDB\_CELL\_TYPE** category has only recently been added by the Broad Institute), we also use a **CELL_TYPE** gene set category, which contains gene sets we derived using the Barres "[Brain RNA-seq](https://www.brainrnaseq.org)" data set, representing sets of marker genes for specific cell types in human and mouse. Sets with names containing `5_times` or `10_times` represent genes whose expression is at least 5 or 10 times greater in that cell type than any other cell type. Sets with `top100` in the name represent the top 100 genes ranked by ratio of expression in that cell type compared to all other cell types.

The results for each differential expression comparison are contained in a sub-folder under the `gene_set_tests` folder (these may also be in further species-specific sub-folders). Each sub-folder contains up to five CSV files, i.e. `*CELL_TYPE_sets.csv`, `*CURATED_sets.csv`, `*MOTIF_sets.csv`, `*GO_sets.csv` and `*-MSIGDB_CELL_TYPE_sets.csv`. These detail, for the particular comparison of experimental conditions, those genes sets in each category which were found by Camera to be significantly differentially expressed, as a whole, when compared to all other genes (with a False Discovery Rate cut-off of **10%**). If there were no significant differentially expressed gene sets for a category, then the file for that category will be missing.

The columns in these CSV files are:

* **GeneSet**: Name of the set of interesting genes.
* **NGenes**: Number of genes in the set.
* **Direction**: Whether expression of the genes in this set is shifted up or down in this particular comparison of experimental conditions.
* **PValue**: Raw p-value for the significance of this shift.
* **FDR**: p-value corrected for multiple testing (i.e. a false discovery rate). Only gene sets with FDR < 0.1 are included (hence all entries in these spreadsheets can be considered significant).

In addition to the CSV files – *if requested* – each differential expression comparison subfolder has further subfolders corresponding to each of the gene set categories, which contain heatmaps of genes in each significant gene set. This allows visualisation of gene expression across genes in each significant set (relative to mean expression for each gene), and how they change as a group across samples included in the comparison. In addition (again, if requested), files can be provided which allow the behaviour of individual genes in each significant gene set to be analysed in depth – however, because these files are very large, they are not provided by default.

_n.b._ for more information about the provenance of a particular gene set, see the [MSigDb](https://software.broadinstitute.org/gsea/msigdb/) database provided by the Broad Institute (free to access, but registration required I think).

### Reactome

In the sub-folder `differential_expression/reactome`, there are enrichment analyses for pathways defined in the [REACTOME pathway database](https://reactome.org).

These are very similar analyses to those performed for GO annotations. Again we take (i) all genes called differentially expressed at false discovery rate < 0.05, (ii) just the up-regulated genes and (iii) just the down-regulated genes. But now we test whether the differentially expressed genes are enriched in genes annotated as belonging to particular pathways defined in REACTOME.

_n.b._ some REACTOME analysis files might not be present, if there were no, or only a few, differentially expressed genes for a particular differential expression comparison.

### Differential transcript expression

Differential transcript expression results are in the folder `differential_expression_tx`. Analogously to the gene-level differential expression results contained in `differential_expression/de_gene`, the sub-folder `differential_expression_tx/de_tx` contains the main transcript-level differential expression results files:

- `*deseq2_results_tpm*.csv` - One or more CSV files containing results of differential transcript expression analysis.
- `*de_summary_mouse*.csv` - A CSV file containing a summary of the numbers of transcripts detected as differentially expressed (at a false discovery rate of 5%).

Each of the one or more main differential transcript expression results files (there may be more than one if there are results for multiple species, or a single file was too unwieldy) contains the following columns:

* **transcript**: Ensembl transcript ID.
* **transcript_length**: Length of transcript in bases.
* **gene**: Ensembl gene ID.
* **number\_of\_transcript**: Number of isoforms defined for this gene in Ensembl.
* **gene_name**: Ensemble gene name.
* **chromosome**: Chromosome on which the gene lies.
* **description**: Ensembl gene description.
* **entrez_id**: NCBI Gene ID (this is mainly for our internal use).
* **gene_type**: e.g. protein coding, pseudogene, LINC RNA etc. etc.
* **\*_tpm**: Per-sample transcript abundances measured in **t**ranscripts **p**er **m**illion (n.b. if Sargasso has been used, there may be separate pre- and post-Sargasso TPM columns for single-species samples).
* **\*\_avg\_tpm**: Per-condition TPMs averaged over the samples in that condition.

and for each differential expression comparison:

* **\<comparison\_name\>_l2fc**: log2 fold change in transcript expression between conditions.
* **\<comparison\_name\>_pval**: raw p-value for statistical significance of differential expression.
* **\<comparison\_name\>_padj**: adjusted p-value (i.e. false discovery rate) after correcting for multiple testing.
* **\<comparison\_name\>\_raw\_l2fc**: raw log2 fold change in transcript expression between conditions, calculated directly from normalised read counts (this may be different to the fold change calculated by DESeq2, which does some further processing beyond the raw value).

The directory `differential_expression_tx/graphs` contains a range of plots describing the transcript-level QC and differential expression results, analogous to the gene-level differential expression plots in `differential_expression/graphs`.

### rMATS

_Differential splicing analysis is optional - let us know if you would like it to be performed._

[rmats](http://rnaseq-mats.sourceforge.net/) (_Shen et al._, "[rMATS: Robust and flexible detection of differential alternative splicing from replicate RNA-Seq data](http://www.pnas.org/content/111/51/E5593)", PNAS (2014)) is a tool to detect differential alternative splicing events from RNA-Seq data. 
 
The statistical model of MATS calculates the p-value and false discovery rate that the difference in the "isoform ratio" of a gene between two experimental conditions exceeds a given user-defined threshold. We run rMATS for each comparison between conditions for which we produce differential gene expression results (e.g. mutant vs control). rMATS takes the samples from the two conditions in the comparison and works out the statistically significant (FDR < 0.05) alternative splicing events.

In the rMATS results folder (`differential_expression/de_rmats`), you will find one or more `AS_summary_{species}.csv` files, which contain brief summaries of all the comparisons made. Each comparison has 5 rows of summary data, corresponding to the 5 different types of splicing events supported by rMATS:

  * skipped exons (SE) 
  * alternative 5′ splice sites (A5SS)
  * alternative 3′ splice sites (A5SS)
  * mutually exclusive exons (MXE)
  * retained introns (RI) 
 
Most of the columns are fairly self-explanatory. The **Up_regulated**, **Down_regulated** and **D.E.total** colums show the number of significant up- and down- regulated AS events, plus the total (at the false discovery rate specified by the column **p.adj.cutoff**). The **Up\_regulated_gene**, **Down\_regulated_gene** and **D.E.total_gene** columns show the corresponding number of genes in which the above events are found (i.e. there may be more than one alternative splicing event per-gene). Note that each AS "event" comprises constitutive exon sequence (i.e. isoform sequence which is present whether the splicing event has taken place or not) and alternatively spliced exon sequence (i.e. isoform sequence which is present when the splicing event has taken place):

<img align="center" src="alternative_splicing_events.png">

In the rMATS result folder you will also find a sub-folder for each comparison (these may be in further species-specific sub-folders if data from multiple species are being sequenced), which in turn holds CSV files containing the details of the results for each of the five supported types of AS events. The columns in these files are:

* **gene**, **gene_name**, **chromosome**, **description**, **entrez_id**, **gene_type**, **gene_length**, **max\_transcript_length**: similar columns to those that appear in the differential gene expression results files
* **avg_fpkm**: Average gene FPKM measured across all samples involved in this comparison (may be useful for filtering rMATS results)
* **ID**: rMATS defines a database of "known" alternative splicing events for which it calculates p-values - this ID identifies the splicing event in the rMATS database.
* **strand**: the strand of the splicing event (+ or -)
* **exonStart_0base**, **exonEnd** etc., or equivalent: there follow a number of columns specific to the particular type of alternative splicing event which define its genomic location; for example, for skipped exons, these define the locations of the upstream, downstream, and (potentially) skipped exon.
* **IJC\_SAMPLE_1**: Numbers correspond to the samples listed for the "base" condition in the rMATS summary file; each is the number of reads in that sample that support the alternatively-spliced isoform sequence.
* **SJC\_SAMPLE_1**: Numbers correspond to the samples listed for the "base" condition in the rMATS summary file; each is the number of reads in that sample that support the non-alternatively-spliced isoform sequence.
* **IJC\_SAMPLE_2**: Numbers correspond to the samples listed for the "comparison" condition in the rMATS summary file; each is the number of reads in that sample that support the alternatively-spliced isoform sequence.
* **SJC\_SAMPLE_2**: Numbers correspond to the samples listed for the "comparison" condition in the rMATS summary file; each is the number of reads in that sample that support the non-alternatively-spliced isoform sequence.
* **IncLevel1**: Numbers correspond to the samples listed for the "base" condition in the rMATS summary file; each is an "inclusion level" (this nomenclature only really makes sense for the skipped exon type of AS event, but is used for all types) - that is a measure of the frequency with which the alternative splicing event takes place in isoforms which contain the constitutive isoform sequence. These values are calculated from the read counts in fields **IJC\_SAMPLE_1** and **SJC\_SAMPLE_1**, but include a correction for the length of the respective constitutive and AS sequences. A value of 0 means the alternatively-spliced sequence is never used, while 1 means the AS sequence is always used. For a good explanation, see [here](https://www.biostars.org/p/256949/).
* **IncLevel2**: Numbers correspond to the samples listed for the "comparison" condition in the rMATS summary file; these are inclusion levels for the comparison condition samples. 
* **IncLevelDifference**: Equals average(**IncLevel1**) - average(**IncLevel2**). That is, a positive value means more use of the alternative-spliced sequence in the base condition, while a negative value means more use of the alternatively-spliced sequence in the comparison condition.
* **PValue**: A likelihood-ratio test is used to calculate a p-value that **IncLevelDifference** is greater than some user-specified value; we have selected the value 0.1, with the belief that inclusion level differences smaller than this are unlikely to be of great biological significance.
* **FDR**: A false discovery rate calculated from the p-value.

By default, rMATS outputs _a lot_ of data, and we have observed that it has a tendency to assign very low p-values to AS events supported by very few reads. To attempt to remove some of this noise, we filter out all events supported by an average read count (across the columns **IJC\_SAMPLE_1**, **SJC\_SAMPLE_1**, **IJC\_SAMPLE_2** and **SJC\_SAMPLE_2**) of less than 5.
 
