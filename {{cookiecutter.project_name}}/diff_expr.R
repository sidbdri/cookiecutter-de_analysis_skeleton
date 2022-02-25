SPECIES <- "unknown_species"
META_DATA=stringr::str_c('meta_data_',SPECIES,'.R')
source(META_DATA)

# Note that when comparisons are run in parallel in RStudio, the output is silent and the R session 
# will be hung until all sub-processes finish or are terminated. When running on the command line, 
# we will see the output but in a random order from each core. Thus, we might want to turn off 
# parallel when debugging in RStudio.
start_parallel(NUM_CORES)
#stop_parallel()

{% if cookiecutter.qSVA !="no" %}
qSVA <- TRUE
{% else %}
qSVA <- FALSE
{% endif %}
PLOT_TO_FILE <- TRUE

MISASSIGNMENT_PERCENTAGE <- MISASSIGNMENT_SAMPLE_REFERENCE_TABLE %>% nrow() > 0

OUTPUT_DIR <- 'results/differential_expression/'
dir.create(file.path(OUTPUT_DIR, "de_gene"), recursive = TRUE)

GRAPHS_DIR <- 'results/differential_expression/graphs/'
dir.create(GRAPHS_DIR, recursive = TRUE)

#### Quality control ####

## we check read distribution for sample region
## for a large region or large number of samples, increase the bin_width to 10, 100 or 1000
## For multiple samples, the plot will order by the input samples order.
# g <- check_sample_bam(samples=SAMPLE_DATA  %>% arrange(CultureType) %>% pull(sample_name),
#                       species=SPECIES,chr='2',start= 75505857,end=75534985,bin_width = 1000)
# plot(g)
# # g + scale_y_continuous(trans='log10')


total_dds_data <- get_total_dds(SAMPLE_DATA, SPECIES, qSVA = qSVA, design_formula = ~1)
total_vst <- total_dds_data %>% varianceStabilizingTransformation(blind = TRUE)

start_plot("pca_all_samples")
total_vst %>% plot_pca(intgroup = FEATURES_FOR_ALL_SAMPLES_PCA,output_data_table_path=file.path(GRAPHS_DIR,str_c('pca_all_samples_',SPECIES,'.csv'))) %>% print()
end_plot()

# scale the pdf base on number of features to be plotted
num_features <- SAMPLE_DATA %>% dplyr::select(-contains('species'), -contains('sample_name')) %>% colnames() %>% length()
start_plot("pca_all_features_all_samples",num_plots=num_features)

if (global_exists('patchworkplot')) {
  rm_global('patchworkplot')
}

# This is to plot individually every feature defined in the SAMPLE_DATA table
SAMPLE_DATA %>% dplyr::select(-species,-sample_name) %>% colnames() %>%
  walk(function(feature) {
    total_vst %>%
      plot_pca(intgroup = c(feature), FALSE) %>%
      add_to_patchwork(plot_var_name = 'patchworkplot')
})

patchworkplot %>% print()
end_plot()

start_plot("heatmap_all_samples")
total_vst %>% plot_heat_map(
  SAMPLE_DATA %>%
    tibble::rownames_to_column("tmp_sample_name") %>%
    tidyr::unite(col='sample_info', c(tmp_sample_name, FEATURES_FOR_ALL_SAMPLES_HEATMAP),
                 sep = ":", remove = FALSE) %>%
    extract2("sample_info"))
end_plot()

start_plot("count_distribution_norm")
plot_count_distribution(total_dds_data, norm=T)
end_plot()

start_plot("count_distribution")
plot_count_distribution(total_dds_data, norm=F)
end_plot()

#### Differential gene expression ####

gene_info <- get_gene_info(SPECIES)
gene_lengths <- get_gene_lengths(SPECIES)

results <- total_dds_data %>% get_count_data()

fpkms <- results %>%
  get_fpkms(gene_lengths, colnames(results) %>% tail(-1), "_fpkm")

results %<>%
  left_join(fpkms) %>%
  left_join(gene_info) %>%
  left_join(gene_lengths)

# Generate plots of the FPKMs of marker genes in the samples
check_cell_type(results, fpkm_check_cutoff = 5, print_check_log = TRUE, print_fpkm_table = FALSE)

# Run all get_res() functions in parallel.
# For debugging, it may be worth calling stop_parallel(), because the mclapply has a problem printing 
# out stdout in rstudio; see:
# http://dept.stat.lsa.umich.edu/~jerrick/courses/stat701/notes/parallel.html#forking-with-mclapply


comparisons_results <- COMPARISON_TABLE %>% pull(comparison) %>%
                  parallelManager(job_name='DESeq2',parallelParam='MulticoreParam',FUN='run_deseq',
                  results_tbl=results,fpkms=fpkms,species=SPECIES,qSVA=qSVA)

if (exists(x = 'all_comparison_pvalue_distribution')) {
  rm(all_comparison_pvalue_distribution)
}

lapply(comparisons_results, function(cmp) {
  # merge the cmp result table into global results table
  get_global("results") %>%
    left_join(cmp$results_tb %>% dplyr::select(gene,contains('.'))) %>%
    set_global("results")

  # merge the cmp summary table into global SUMMARY_TABLE
  get_global("SUMMARY_TB") %>%
    rbind(cmp$summary_tb %>% as.data.frame(stringsAsFactors=FALSE)) %>%
    set_global("SUMMARY_TB")

  # merge the p value plots
  add_to_patchwork(cmp$p_plot, plot_var_name = 'all_comparison_pvalue_distribution')

  # plot top DE gene for each comparison
  start_plot(str_c('top_de_genes_fpkm_',cmp$comparison_name))
  plot_genes_fpkm(results,results %>% arrange(across(str_c(cmp$comparison_name,'.padj',sep = ''))) %>% head(4) %>% pull(gene),print_fpkm_table = F)
  end_plot()

  # export the res and dds
  # commnet out for now.
  # see https://github.com/sidbdri/cookiecutter-de_analysis_skeleton/issues/83
  # cmp$res %>% set_global(cmp$comparison %>% str_c('res', sep = '_'))
  cmp$dds %>% set_global(cmp$comparison %>% str_c('dds', sep = '_'))
  cmp$misassignment_percentage %>% set_global(cmp$comparison %>% str_c('misassignment_percentage', sep = '_'))
  'success'
})

start_plot("all_comparison_pvalue_distribution")
all_comparison_pvalue_distribution
end_plot()

plot_scatter_fpkm(results)


mitochondrial_terms <- c('GO:0005739',ontology_find_all_children_terms('GO:0005739',as.list(GO.db::GOCCCHILDREN))) %>% unique()
nuclear_encoded_mitochondrial_genes<-topGO::annFUN.org("CC", mapping = switch(SPECIES, mouse = "org.Mm.eg.db", rat = "org.Rn.eg.db", human = "org.Hs.eg.db"), ID = "ensembl") %>%
  extract(mitochondrial_terms) %>% unlist() %>% unique()
mitochondrial_genes <- gene_info %>% filter(chromosome=='MT') %>% pull(gene)
start_plot('mitochondrial_gene_count')
plot_gene_percentage(counts(total_dds_data), list(MT_mito=mitochondrial_genes, NUC_mito=nuclear_encoded_mitochondrial_genes),use_percentage = FALSE) %>% print
end_plot()
start_plot('mitochondrial_gene_percentage')
plot_gene_percentage(counts(total_dds_data), list(MT_mito=mitochondrial_genes, NUC_mito=nuclear_encoded_mitochondrial_genes),use_percentage = TRUE) %>% print
end_plot()

# save results
if (COMPARISON_TABLE %>% pull(group) %>% unique() %>% length() > 1) {
  save_results_by_group(results)
}

results %>%
  dplyr::select(
    gene, gene_name, chromosome, description, entrez_id, gene_type,
    gene_length, max_transcript_length,
    everything(), -dplyr::contains("_fpkm"), -dplyr::ends_with(".stat")) %>%
  write_csv(file.path(OUTPUT_DIR, "de_gene", str_c("deseq2_results_count_", SPECIES, ".csv")), na="")

results %>%
  dplyr::select(
    gene, gene_name, chromosome, description, entrez_id, gene_type,
    gene_length, max_transcript_length,
    dplyr::contains("_fpkm"),
    COMPARISON_TABLE %>%
      pull(comparison) %>%
      sapply(FUN = function(x) results %>% colnames() %>% str_which(str_c("^", x, sep =''))) %>%
      unlist() %>%
      as.vector() %>%
      unique(),
    -dplyr::ends_with(".stat")) %>%
  write_csv(file.path(OUTPUT_DIR, "de_gene", str_c("deseq2_results_fpkm_", SPECIES, ".csv")), na="")

SUMMARY_TB %>%
  write_csv(file.path(OUTPUT_DIR, "de_gene", str_c("de_summary_", SPECIES, ".csv")), na="")

#### GO enrichment analysis ####
# For each comparison:
#   - for the GO/Reactome analyses, we are using all/up/down regulated genes,
#   - for GSA, we are using four gene set categories: "CURATED", "MOTIF", "GO" and "CELL_TYPE"
expressed_genes <- get_total_dds(SAMPLE_DATA, SPECIES, filter_low_counts = TRUE) %>%
  get_count_data()


GO_results <- COMPARISON_TABLE%>%
  pull(comparison) %>% expand.grid(c('all','up','down'),c("BP", "MF", "CC")) %>% 
  tidyr::unite(col='job_string',sep = ';') %>% pull(job_string) %>%
  parallelManager(parallelParam='SnowParam',FUN='run_topgo',
                  result_tbl=results,
                  p.adj.cutoff=P.ADJ.CUTOFF,
                  expressed_genes=expressed_genes,
                  species=SPECIES,
                  out_dir=file.path(OUTPUT_DIR,"enrichment_tests"))
#### Reactome pathway analysis
Reactome_results <- COMPARISON_TABLE%>%
  pull(comparison) %>% expand.grid(c('all','up','down')) %>% 
  tidyr::unite(col='job_string',sep = ';') %>% pull(job_string) %>%
  parallelManager(job_name='Reactome',parallelParam='SnowParam',FUN='run_reactome',
                  result_tbl=results,
                  p.adj.cutoff=P.ADJ.CUTOFF,
                  expressed_genes=expressed_genes,
                  species=SPECIES,
                  out_dir=file.path(OUTPUT_DIR,"reactome")) 

#### Gene set enrichment analysis ####
gene_set_categories <- list("CURATED", "MOTIF", "GO", "CELL_TYPE")

list_of_gene_sets <- gene_set_categories %>%
  set_names(.) %>% lapply(function(category, ...) get_gene_sets(SPECIES, category))

GS_results <- COMPARISON_TABLE %>%
  pull(comparison) %>% expand.grid(gene_set_categories) %>% 
  tidyr::unite(col='job_string',sep = ';') %>% pull(job_string) %>%
  parallelManager(job_name='GS',parallelParam='SnowParam',FUN='run_gs',
                  result_tbl=results,
                  dds_list=str_c(COMPARISON_TABLE%>% pull(comparison),'dds',sep = '_') %>% set_names(.) %>% lapply(get_global),
                  list_of_gene_sets=list_of_gene_sets,
                  gene_info=gene_info,
                  species=SPECIES,
                  out_dir=file.path(OUTPUT_DIR,"gene_set_tests")) 

## plot gs heatmap
GSheatmap <- GS_results %>% names %>%
  parallelManager(job_name='GSheatmap',parallelParam='MulticoreParam',FUN='plot_significant_set_heatmap',
                  results_tbl=results,
                  gs_result=GS_results,
                  list_of_gene_set=list_of_gene_sets,
                  comparison_tbl=COMPARISON_TABLE,
                  species=SPECIES,
                  out_dir=file.path(OUTPUT_DIR,"gene_set_tests")) 

#### Saving and loading the workspace ####

# Save the objects in the workspace for future analysis

rws <- "results/Rworkspace/"
if (!dir.exists(rws)) {
  dir.create(rws,recursive=TRUE)
}

save(list = ls() %>% grep(x = ., pattern='comparisons_results', value = T,invert = T),
     file = str_c(rws,"diff_expr_",SPECIES,".RData"))

# Load the save workspace to get all objects back for analysis

# load_rs_data()
