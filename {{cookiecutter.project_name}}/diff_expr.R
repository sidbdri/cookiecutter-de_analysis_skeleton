source("meta_data.R")

SPECIES <- "{{cookiecutter.species}}"

{% if cookiecutter.qSVA !="no" %}
qSVA <- TRUE
{% else %}
qSVA <- FALSE
{% endif %}
PLOT_TO_FILE <- TRUE

OUTPUT_DIR <- 'results/differential_expression/de_gene/'
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive=TRUE)

GRAPHS_DIR <- 'results/differential_expression/graphs/'
if (!dir.exists(GRAPHS_DIR)) dir.create(GRAPHS_DIR, recursive=TRUE)

#####

total_dds_data <- get_total_dds(SAMPLE_DATA, SPECIES, qSVA=qSVA)
total_vst <- total_dds_data %>% varianceStabilizingTransformation

start_plot("pca_all")
total_vst %>% plot_pca_with_labels(intgroup=PCA_FEATURE)
end_plot()


start_plot("pca_features")
#This is to plot individually every feature defined in the SAMPLE_DATA table
SAMPLE_DATA %>% dplyr::select(-species,-sample_name) %>% colnames() %>%
  walk(function(feature){
    total_vst %>% plot_pca(intgroup=c(feature),FALSE) %>%
                  add_to_patchwork(plot_var_name='pathworkplot')
})
pathworkplot
end_plot()

start_plot("heatmap_all")
total_vst %>% plot_heat_map(
  SAMPLE_DATA %>% 
    tibble::rownames_to_column("tmp_sample_name") %>%
    tidyr::unite(col='sample_info', c(tmp_sample_name, HEAT_MAP_FEATURE), 
                 sep = ":", remove = FALSE) %>% 
    extract2("sample_info"))
end_plot()

start_plot("count_distribution_norm")
plot_count_distribution(total_dds_data, norm=T)
end_plot()

start_plot("count_distribution")
plot_count_distribution(total_dds_data, norm=F)
end_plot()

#####

gene_info <- get_gene_info(SPECIES)
gene_lengths <- get_gene_lengths(SPECIES)

results <- total_dds_data %>% get_count_data() 

fpkms <- results %>% 
  get_fpkms(gene_lengths, colnames(results) %>% tail(-1), "_fpkm")

results %<>% 
  left_join(fpkms) %>%
  left_join(gene_info) %>%
  left_join(gene_lengths)

# run all get_res() functions and add to main "results" object
COMPARISON_TABLE %>% pull(comparison) %>% walk (
  function(comparison_name) {
    res <- get_res(comparison_name, fpkms, SPECIES, qSVA=qSVA)
    
    results <- get("results", envir = .GlobalEnv) %>% 
      left_join(res[[1]], by="gene") %>%
      dplyr::rename(!!str_c(comparison_name, '.l2fc') := log2FoldChange,
                    !!str_c(comparison_name, '.raw_l2fc') := raw_l2fc,
                    !!str_c(comparison_name, '.stat') := stat,
                    !!str_c(comparison_name, '.pval') := pvalue,
                    !!str_c(comparison_name, '.padj') := padj)
    
    p_plot<-plot_pvalue_distribution(results, str_c(comparison,'.pval'))

    add_to_patchwork(p_plot,plot_var_name='all_comparison_pvalue_distribution')

    assign("results", results,envir = .GlobalEnv)
    
    comparison_name %>% str_c('res', sep = '_') %>% assign(res, envir = .GlobalEnv)
  }
)

start_plot("all_comparison_pvalue_distribution")
all_comparison_pvalue_distribution
end_plot()

###########
## workout sargasso error ratio
## Calculate proportion of reads incorrectly assigned to other species in pure
## samples for each 1-to-1 orthologous gene
# RAT_ONLY_SAMPLES <- c('a1','a2','a3')
# rat_only_mouse_counts <- RAT_ONLY_SAMPLES %>% 
#   get_single_species_only_counts("mouse") %>% 
#   rename(mouse_gene = gene, rat_only_mouse_count=total)

# rat_only_rat_counts <- RAT_ONLY_SAMPLES %>% 
#   get_single_species_only_counts("rat") %>% 
#   rename(rat_gene = gene, rat_only_rat_count=total)

# ortholog_info <- get_mouse_rat_ortholog_info()

# rat_only_mapping_info <- ortholog_info %>%
#   left_join(rat_only_mouse_counts) %>%
#   left_join(rat_only_rat_counts) %>%
#   mutate(rat_only_total_count=rat_only_mouse_count + rat_only_rat_count,
#          rat_only_mouse_frac=rat_only_mouse_count/rat_only_total_count)
#
### easier do the join when writing result to csv
# results %<>% left_join(rat_only_mapping_info, by=c("gene" = "rat_gene"))
###########

# save results
results %>% 
  dplyr::select(
    gene, gene_name, chromosome, description, entrez_id, gene_type,
    gene_length, max_transcript_length,
    everything(), -dplyr::contains("_fpkm"), -dplyr::ends_with(".stat")) %>%
  write_csv(str_c(OUTPUT_DIR, "/deseq2_results_count_", SPECIES, ".csv"))

results %>% 
  dplyr::select(
    gene, gene_name, chromosome, description, entrez_id, gene_type,
    gene_length, max_transcript_length,
    dplyr::contains("_fpkm"), 
    COMPARISON_TABLE %>% 
      pull(comparison) %>%
      sapply(FUN = function(x) results %>% colnames() %>% str_which(str_c("^", x, sep =''))) %>%
      as.vector() %>% 
      unique(), 
    -dplyr::ends_with(".stat")) %>%
  write_csv(str_c(OUTPUT_DIR, "/deseq2_results_fpkm_", SPECIES, ".csv"))

SUMMARY_TB %>%
  write_csv(str_c(OUTPUT_DIR, "/de_summary_", SPECIES, ".csv"))

##### GO analyses

expressed_genes <- get_total_dds(SAMPLE_DATA, SPECIES, filter_low_counts=TRUE) %>% 
  get_count_data()

COMPARISON_TABLE %>% pull(comparison) %>% walk(function(comparison_name) {
  p_str <- str_c(comparison_name, '.padj')
  l2fc_str <- str_c(comparison_name, '.l2fc')
  
  results <- get("results",envir = .GlobalEnv)
  
  results %>% 
    filter(get(p_str) < 0.05) %>% 
    perform_go_analyses(expressed_genes, comparison_name, SPECIES)
  
  results %>%
    filter(get(p_str) < 0.05 & get(l2fc_str) > 0) %>% 
    perform_go_analyses(expressed_genes, str_c(comparison_name, '.up'), SPECIES)
  
  results %>%
    filter(get(p_str) < 0.05 & get(l2fc_str) < 0) %>% 
    perform_go_analyses(expressed_genes, str_c(comparison_name, '.down'), SPECIES)
})

##### Gene set enrichment analysis

gene_set_categories <- list("CURATED", "MOTIF", "GO")

list_of_gene_sets <- gene_set_categories %>% 
  map(function(category) get_gene_sets(SPECIES, category))

COMPARISON_TABLE %>% pull(comparison) %>% walk(function(comparison_name) {
  res <- str_c(comparison_name, 'res', sep = '_') %>% get(envir = .GlobalEnv)
  
  camera_results <- list_of_gene_sets %>% 
    map(function(category_gene_sets) {
      get_camera_results(res[[2]], category_gene_sets, gene_info)
    })
  
  assign(str_c(comparison_name, 'camera_results', sep = '_'), 
         camera_results, envir = .GlobalEnv)
  
  for (category in seq(1:length(gene_set_categories))) {
    de_res <- results %>% dplyr::select(
      gene, gene_name, entrez_id, 
      starts_with(str_c(comparison_name, ".")), 
      -starts_with(str_c(comparison_name, ".stat")))  
    write_camera_results(
      gene_set_categories[[category]], list_of_gene_sets[[category]], 
      comparison_name, SPECIES,
      de_res, camera_results[[category]])
  }
})

# results %>% plot_gene_set(list_of_gene_sets[[3]], "GO_<go_term>", "condition.stat")
# results %>% get_gene_set_results(list_of_gene_sets[[3]], "GO_<go_term>", "condition.pval") %>% head
