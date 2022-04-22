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

    # for interaction we remove raw_l2fc column
    if(comparison_name %in% INTERACTION_TABLE$comparison){
        results_tb %<>% dplyr::select(-str_c(comparison_name, '.raw_l2fc'))
    }

    P=NULL
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

      res$summary_tb_row %<>% as.data.frame() %>%
        mutate(Misassignment_samples_in_comparison_level_condition = P$condition_reference_samples) %>%
        mutate(Misassignment_samples_in_base_level_condition = P$condition_base_reference_samples)

    }

    p_plot <- plot_pvalue_distribution(results_tb, str_c(comparison_name,'.pval'))

    ## return the results and merge them later
    list(comparison_name = comparison_name,
         res = res$res,
         dds = res$dds,
         results_tb = results_tb,
         summary_tb = res$summary_tb_row,
         p_plot = p_plot,
         misassignment_percentage=P)
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
# Thus we need to reduce the number of comparisons we run in parallel to ensure we are not using more cores
# than specified. The total number of cores used after the following line will be 3 * getOption("mc.cores")
if (PARALLEL) {
  adjust_parallel_cores()
}

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

        perform_go_analyses(r, expressed_genes, comparison_name, cmp, SPECIES, out_dir = file.path(OUTPUT_DIR,"enrichment_tests"))
      }
    ) %>% set_names(str_c(comparison_name,c('.all','.up','.down')))
  }
)

#### Reactome pathway analysis

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

#### Gene set enrichment analysis ####

gene_set_categories <- list("CURATED", "MOTIF", "GO", "CELL_TYPE", "MSIGDB_CELL_TYPE")

list_of_gene_sets <- gene_set_categories %>%
  set_names(.) %>% lapply(function(category, ...) get_gene_sets(SPECIES, category))

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
          de_res, camera_results[[category]], out_dir = file.path(OUTPUT_DIR,"gene_set_tests"))
      }
    )

    camera_results
  }
)

# get entrez IDs from significantly DE gene sets
significant_entrez_ids<-function(gene_set_category, sets_per_comparison){
  # filter by significance
  significant<-sets_per_comparison[[gene_set_category]] %>% filter(FDR < 0.05) %>% rownames()

  # use the globally defined gene set category to get all the entrez IDs in that gene set
  significant_entrez<-sapply(significant, function(x) list_of_gene_sets[[gene_set_category]][[x]])
  return(significant_entrez)
}

# loop through comparisons
for (comp in COMPARISON_TABLE$comparison) {
  # get the gene set results from camera for that comparison
  comparison_list <- GS_results[comp]

  # because of the way it seems to be arranged, the comparison name is used to get all the
  # gene set categories (e.g. GO, CELL TYPE etc.) for that comparison
  gene_sets_per_comparison <- comparison_list[[comp]]

  # for each gene set category, defined elsewhere globally
  # get the entrez ids of the gene sets which are significant e.g. FDR < 0.05
  # replace names using the gene set categories as these seem to be removed during the apply
  significant_entrez <- sapply(gene_set_categories, significant_entrez_ids, gene_sets_per_comparison, simplify = FALSE)
  names(significant_entrez) <- unlist(gene_set_categories)

  #get the comparison criteria for this comparison
  comparison_criteria <- COMPARISON_TABLE %>% filter(comparison == comp) %>% dplyr::select(condition_name, condition, condition_base, filter)
  selected_conditions <- comparison_criteria %>% dplyr::select(condition, condition_base) %>% as.character()

  # pull out samples which have the conditions in this comparison
  # and samples which match the filter specified
  condition_name <- comparison_criteria %>% dplyr::select(condition_name) %>% pull()
  samples_in_comparison <- SAMPLE_DATA %>%
    filter(!!parse_expr(condition_name) %in% selected_conditions &
             !!(parse_expr(comparison_criteria %>% dplyr::select(filter) %>% pull()))) %>%
    dplyr::select(sample_name, !!parse_expr(condition_name)) %>%
    arrange(!!parse_expr(condition_name))

  # loop over significant entrez for gene set categories
  # plot heatmap using the sample names corresponding to the relevant comparison
  sapply(names(significant_entrez), plot_significant_set_heatmap, significant_entrez, comp, samples_in_comparison, SPECIES)
}

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
