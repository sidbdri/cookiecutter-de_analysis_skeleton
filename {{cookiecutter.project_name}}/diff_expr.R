SPECIES <- "unknown_species"
META_DATA=str_c('meta_data_',SPECIES,'.R')
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

#####

total_dds_data <- get_total_dds(SAMPLE_DATA, SPECIES, qSVA = qSVA, design_formula = ~1)
total_vst <- total_dds_data %>% varianceStabilizingTransformation(blind = TRUE)

start_plot("pca_all")
total_vst %>% plot_pca(intgroup = PCA_FEATURE) %>% print()
end_plot()

# scale the pdf base on number of features to be plotted
num_features <- SAMPLE_DATA %>% dplyr::select(-contains('species'), -contains('sample_name')) %>% colnames() %>% length()
start_plot("pca_features",num_plots=num_features)

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

# Generate plots of the FPKMs of marker genes in the samples
check_cell_type(results, fpkm_check_cutoff = 5, print_check_log = TRUE, print_fpkm_table = FALSE)

# Run all get_res() functions in parallel.
# For debugging, it may be worth calling stop_parallel(), because the mclapply has a problem printing 
# out stdout in rstudio; see:
# http://dept.stat.lsa.umich.edu/~jerrick/courses/stat701/notes/parallel.html#forking-with-mclapply
comparisons_results <- COMPARISON_TABLE %>% pull(comparison) %>% set_names(.) %>%  lapply_fork(
  function(comparison_name) {
    res <- get_res(comparison_name, fpkms, SPECIES, qSVA = qSVA)
    
    results_tb <- get_global("results") %>% 
      left_join(res[[1]], by = "gene") %>%
      dplyr::rename(!!str_c(comparison_name, '.l2fc') := log2FoldChange,
                    !!str_c(comparison_name, '.raw_l2fc') := raw_l2fc,
                    !!str_c(comparison_name, '.stat') := stat,
                    !!str_c(comparison_name, '.pval') := pvalue,
                    !!str_c(comparison_name, '.padj') := padj)
    
    if (MISASSIGNMENT_PERCENTAGE) {
      P <- get_misassignment_percentages(comparison_name, gene_lengths)
      
      if (!is.na(P$condition_reference_samples)) {
        results_tb %<>% left_join(
          P$P_condition %>% 
            dplyr::select(gene, !!str_c(comparison_name, '.perc.', COMPARISON_TABLE %>% 
                                          filter(comparison == comparison_name) %>% 
                                          pull(condition)) := p))
      }   
      
      if (!is.na(P$condition_base_reference_samples)) {
        results_tb %<>% left_join(
          P$P_condition_base %>% 
            dplyr::select(gene,!!str_c(comparison_name, '.perc.', COMPARISON_TABLE %>% 
                                         filter(comparison == comparison_name) %>% 
                                         pull(condition_base)) := p))
      }   
      
      res$summary_tb_row %<>% 
        mutate(Misassignment_samples_in_comparison_level_condition = 
                 ifelse(Comparison == comparison_name, 
                        P$condition_reference_samples %>% str_c(collapse = ','),
                        Misassignment_samples_in_comparison_level_condition)) %>% 
        mutate(Misassignment_samples_in_base_level_condition = 
                 ifelse(Comparison == comparison_name,
                        P$condition_base_reference_samples %>% str_c(collapse = ','),
                        Misassignment_samples_in_base_level_condition))
      
    }    
    
    p_plot <- plot_pvalue_distribution(results_tb, str_c(comparison_name,'.pval'))
    
    ## return the results and merge them later
    list(comparison_name = comparison_name,
         res = res$res,
         dds = res$dds,
         results_tb = results_tb,
         summary_tb = res$summary_tb_row,
         p_plot = p_plot)
  }
)

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
  
  # export the res and dds
  # commnet out for now.
  # see https://github.com/sidbdri/cookiecutter-de_analysis_skeleton/issues/83
  # cmp$res %>% set_global(cmp$comparison %>% str_c('res', sep = '_'))
  cmp$dds %>% set_global(cmp$comparison %>% str_c('dds', sep = '_'))

  'success'
}) 


start_plot("all_comparison_pvalue_distribution")
all_comparison_pvalue_distribution
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

#####

# For each comparison: 
#   - for the GO/Reactome analyses, we are using all/up/down regulated genes,
#   - for GSA, we are using three gene set categories: "CURATED", "MOTIF" and "GO"
# Thus we need to reduce the number of comparisons we run in parallel to ensure we are not using more cores
# than specified. The total number of cores used after the following line will be 3 * getOption("mc.cores")
if (PARALLEL) {
  adjust_parallel_cores()
}

##### GO analysis

expressed_genes <- get_total_dds(SAMPLE_DATA, SPECIES, filter_low_counts = TRUE) %>% 
  get_count_data()

GO_results <- COMPARISON_TABLE %>% 
  pull(comparison) %>% 
  set_names(.) %>% 
  lapply_socket(X=., function(comparison_name) {
    p_str <- str_c(comparison_name, '.padj')
    l2fc_str <- str_c(comparison_name, '.l2fc')
    
    results <- get_global("results")

    lapply_socket(cores = 3, X = c('', '.up', '.down'), function(cmp) {
        if (cmp == '.up') {
          r <- results %>% filter(get(p_str) < P.ADJ.CUTOFF & get(l2fc_str) > 0)
        } else if((cmp == '.down')) {
          r <- results %>% filter(get(p_str) < P.ADJ.CUTOFF & get(l2fc_str) < 0)
        } else {
          r <- results %>% filter(get(p_str) < P.ADJ.CUTOFF)
        }
        
        perform_go_analyses(r, expressed_genes, comparison_name, cmp, SPECIES, out_dir = file.path(OUTPUT_DIR,"go"))
      }
    ) %>% set_names(str_c(comparison_name,c('.all','.up','.down')))
  }
)

##### Reactome pathway analysis

Reactome_results<- COMPARISON_TABLE %>% 
  pull(comparison) %>% 
  set_names(.) %>% 
  lapply_socket(X = ., function(comparison_name) {
    p_str <- str_c(comparison_name, 'padj', sep = '.')
    l2fc_str <- str_c(comparison_name, 'l2fc', sep = '.')
    
    results <- get_global("results")

    lapply_socket(cores = 3, X = c('', '.up', '.down'), function(cmp) {
        if (cmp=='.up') {
          r <- results %>% filter(get(p_str) < P.ADJ.CUTOFF  & get(l2fc_str) > 0)
        } else if((cmp=='.down')) {
          r <- results %>% filter(get(p_str) < P.ADJ.CUTOFF  & get(l2fc_str) < 0)
        } else {
          r <- results %>% filter(get(p_str) < P.ADJ.CUTOFF)
        }
        perform_pathway_enrichment(r, expressed_genes, comparison_name, cmp, SPECIES, out_dir = file.path(OUTPUT_DIR,"reactome"))
      }
    ) %>% set_names(str_c(comparison_name,c('.all', '.up', '.down')))
  }
)

##### Gene set enrichment analysis

gene_set_categories <- list("CURATED", "MOTIF", "GO")

list_of_gene_sets <- gene_set_categories %>% 
  set_names(.) %>% 
  lapply_fork(cores = length(gene_set_categories), X = ., 
              function(category, ...) get_gene_sets(SPECIES, category))

GS_results <- COMPARISON_TABLE %>% 
  pull(comparison) %>% 
  set_names(.) %>% 
  lapply_fork(X=., function(comparison_name, ...) {
    dds <- str_c(comparison_name, 'dds', sep = '_') %>% get_global()
    
    camera_results <- list_of_gene_sets %>% 
      map(function(category_gene_sets) {
        get_camera_results(dds, category_gene_sets, gene_info)
      })
  
    lapply_fork(cores = 3, X = seq(1:length(gene_set_categories)), function(category,...) {
        de_res <- results %>% dplyr::select(
          gene, gene_name, entrez_id, 
          starts_with(str_c(comparison_name, ".")), 
          -starts_with(str_c(comparison_name, ".stat")))  
        write_camera_results(
          gene_set_categories[[category]], list_of_gene_sets[[category]], 
          comparison_name, SPECIES,
          de_res, camera_results[[category]], out_dir = file.path(OUTPUT_DIR,"gsa"))
      }
    ) 
    
    camera_results
  }
)

# results %>% plot_gene_set(list_of_gene_sets[[3]], "GO_<go_term>", "condition.stat")
# results %>% get_gene_set_results(list_of_gene_sets[[3]], "GO_<go_term>", "condition.pval") %>% head

##### Saving and loading the workspace

# Save the objects in the workspace for future analysis

rws <- "results/Rworkspace/"
if (!dir.exists(rws)) {
  dir.create(rws,recursive=TRUE)
}

save(list = ls() %>% grep(x = ., pattern='comparisons_results', value = T,invert = T),
     file = str_c(rws,"diff_expr.RData"))
  
# Load the save workspace to get all objects back for analysis

# load_rs_data()
