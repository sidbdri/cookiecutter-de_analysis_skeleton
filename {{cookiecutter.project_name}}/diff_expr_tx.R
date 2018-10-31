source("meta_data.R")

SPECIES <- "{{cookiecutter.species}}"

TX_LEVEL <- FALSE
QUANT_METHOD <- 'salmon'
USE_TX <- TRUE

{% if cookiecutter.qSVA !="no" %}
qSVA <- TRUE
{% else %}
qSVA <- FALSE
{% endif %}

MISASSIGNMENT_PRECENTAGE = MISASSIGNMENT_SAMPLE_REFERENCE_TABLE %>% nrow() >0


OUTPUT_DIR <- 'results/differential_expression/de_tx'
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR,recursive=TRUE)

GRAPHS_DIR <- 'results/differential_expression/graphs/'
if (!dir.exists(GRAPHS_DIR)) dir.create(GRAPHS_DIR, recursive=TRUE)

#####

total_dds_data <- get_total_dds_tximport(SAMPLE_DATA, QUANT_METHOD, TX_LEVEL)
total_vst <- total_dds_data %>% varianceStabilizingTransformation

start_plot("pca_all_tx")
total_vst %>% plot_pca_with_labels(intgroup=c("condition"))
end_plot()

start_plot("pca_features_tx")
#This is to plot individually every feature defined in the SAMPLE_DATA table
SAMPLE_DATA %>% dplyr::select(-species,-sample_name) %>% colnames() %>%
  walk(function(feature){
    total_vst %>% plot_pca(intgroup=c(feature),FALSE) %>%
                  add_to_patchwork(plot_var_name='pathworkplot')
})
pathworkplot
end_plot()

start_plot("heatmap_all_tx")
total_vst %>% plot_heat_map(SAMPLE_DATA %>% 
                            mutate(sample_info=str_c(condition, ..., sep=":")) %>% 
                            extract2("sample_info"))
end_plot()

start_plot("count_distribution_norm_tx")
plot_count_distribution(total_dds_data, norm=F)
end_plot()

start_plot("count_distribution_tx")
plot_count_distribution(total_dds_data, norm=T)
end_plot()

#####

gene_info <- get_gene_info(SPECIES)
gene_lengths <- get_gene_lengths(SPECIES)

results <- total_dds_data %>% get_count_data()

txi <- get_tximport(SAMPLE_DATA, QUANT_METHOD, TX_LEVEL)
tpms <- txi$abundance %>% 
  as.data.frame() %>%
  rename_all(funs(str_c(., "_tpm")))

if (TX_LEVEL) {
    tpms %<>% tibble::rownames_to_column('transcript')
    tpms %<>% inner_join(get_transcripts_to_genes(SPECIES)) %>%
        left_join(get_avg_tpm(tpms, TX_LEVEL))

    results %<>% rename(transcript=gene) %>%
        left_join(txi$Length %>% 
                    rename(transcript_length=length) %>% 
                    tibble::rownames_to_column('transcript')) %>%
        left_join(get_transcripts_to_genes(SPECIES)) %>%
        left_join(tpms) %>%
        left_join(gene_info) %>%
        left_join(gene_lengths)

    # number of transcripts
    results %<>% left_join(
      results %>%
        dplyr::select(gene,transcript) %>%
        group_by(gene) %>% summarise(number_of_transcript = n()))
} else {
    tpms %<>% tibble::rownames_to_column('gene')
    tpms %<>% inner_join(get_transcripts_to_genes(SPECIES)) %>%
        left_join(get_avg_tpm(tpms,TX_LEVEL))

    results %<>%
        left_join(tpms) %>%
        left_join(gene_info) %>%
        left_join(gene_lengths)
}

# #workout avg tpm
# results %<>% dplyr::select(gene,dplyr::contains("_tpm")) %>% group_by(gene) %>%
#                 summarise_all(.funs = sum) %>%
#                 mutate(sum=rowSums(dplyr::select(., dplyr::contains("_tpm"))),
#                 n=ncol(dplyr::select(., dplyr::contains("_tpm")))) %>%
#                 mutate(avg_gene_tpm=sum/n) %>%
#                 dplyr::select(gene,avg_gene_tpm) %>% right_join(results)

# run all get_res() functions and add to main "results" object
COMPARISON_TABLE %>% pull(comparison) %>% walk (
  function(comparison_name) {
    res <- get_res(comparison_name, tpms, use_tx=USE_TX, 
                   quant_method=QUANT_METHOD, tx_level=TX_LEVEL)

    results <- get("results", envir = .GlobalEnv) %>%
      left_join(res[[1]]) %>%
      dplyr::rename(!!str_c(comparison_name, '.l2fc') := log2FoldChange,
                    !!str_c(comparison_name, '.raw_l2fc') := raw_l2fc,
                    !!str_c(comparison_name, '.stat') := stat,
                    !!str_c(comparison_name, '.pval') := pvalue,
                    !!str_c(comparison_name, '.padj') := padj)

    ##work out misassigned precentage
    if(MISASSIGNMENT_PRECENTAGE){
        P<-get_misassigned_precentage(comparison_name)

        results_sargasso %<>% left_join(P$P_condition %>% dplyr::select(gene,!!str_c(comparison_name, '.perc.',COMPARISON_TABLE %>% filter(comparison==comparison_name) %>% pull(condition)) := p))
        results_sargasso %<>% left_join(P$P_condition_base %>% dplyr::select(gene,!!str_c(comparison_name, '.perc.',COMPARISON_TABLE %>% filter(comparison==comparison_name) %>% pull(condition_base)) := p))

        SUMMARY_TB <- get("SUMMARY_TB", envir = .GlobalEnv) %>%
            mutate(Misassignment_samples_in_comparison_level_condition = ifelse(Comparison == comparison_name,P$condition_reference_samples %>% str_c(collapse = ','),Misassignment_samples_in_comparison_level_condition)) %>%
            mutate(Misassignment_samples_in_base_level_condition = ifelse(Comparison == comparison_name,P$condition_base_reference_samples %>% str_c(collapse = ','),Misassignment_samples_in_base_level_condition))

        assign("SUMMARY_TB", SUMMARY_TB,envir = .GlobalEnv)
    }

    p_plot<-plot_pvalue_distribution(results, str_c(comparison_name,'.pval'))

    add_to_patchwork(p_plot,plot_var_name='all_comparison_pvalue_distribution')

    assign("results", results, envir = .GlobalEnv)

    comparison_name %>% str_c('res', sep = '_') %>% assign(res, envir = .GlobalEnv)
  }
)

start_plot("all_comparison_pvalue_distribution")
all_comparison_pvalue_distribution
end_plot()

# save results

if (TX_LEVEL) {
  columns_included <- c('transcript', 'transcript_length', 'gene', 'number_of_transcript')
  tx_level_str <- "transcript"
} else {
  columns_included <- c('gene')
  tx_level_str <- "gene"
}

results %>%
  dplyr::select(
    !!columns_included, gene_name, chromosome, description, entrez_id, gene_type,
    everything(), -gene_length, -max_transcript_length, 
    -dplyr::contains("_tpm"), -dplyr::ends_with(".stat")) %>%
  write_csv(str_c(OUTPUT_DIR, "/deseq2_results_count_", SPECIES, "_tx_",
            tx_level_str, "_", QUANT_METHOD, ".csv"))

if (TX_LEVEL) {
  results %<>% select(-gene_length, -max_transcript_length)
}

results %>% 
  dplyr::select(
    !!columns_included, gene_name, chromosome, description, entrez_id, gene_type,
    dplyr::contains("_tpm"), 
    COMPARISON_TABLE %>% 
      pull(comparison) %>%
      sapply(FUN = function(x) results %>% colnames() %>% str_which(str_c("^",x,sep =''))) %>%
      as.vector() %>% 
      unique(), 
    -dplyr::ends_with(".stat")) %>% 
  write_csv(str_c(OUTPUT_DIR, "/deseq2_results_tpm_", SPECIES, "_tx_",
            tx_level_str, "_", QUANT_METHOD, ".csv"))

SUMMARY_TB %>%
  write_csv(str_c(OUTPUT_DIR, "/de_summary_", SPECIES, "_tx_",
            tx_level_str, "_", QUANT_METHOD, ".csv"))

# library (knitr)
# sink(str_c("results/differential_expression/de_summary_",
#             SPECIES,"_tx_",tx_level_str,"_",QUANT_METHOD,".md"))
# kable(SUMMARY_TB, format = 'markdown')
# sink()

##### GO and GSA analyses

if (!USE_TX | !TX_LEVEL) {
  expressed_genes <- total_dds_data %>% get_count_data()

   COMPARISON_TABLE %>% pull(comparison) %>% walk(function(comparison_name) {
     p_str <- str_c(comparison_name, '.padj')
     l2fc_str <- str_c(comparison_name ,'.l2fc')

     results <- get("results", envir = .GlobalEnv)

     results %>%
       filter(get(p_str) < 0.05) %>%
       perform_go_analyses(expressed_genes, comparison_name, SPECIES)

     results %>%
      filter(get(p_str) < 0.05  & get(l2fc_str) > 0) %>%
      perform_go_analyses(expressed_genes, str_c(comparison_name, '.up'), SPECIES)

    results %>%
      filter(get(p_str) < 0.05  & get(l2fc_str) < 0) %>%
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
      }
    )

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
}

# results %>% plot_gene_set(list_of_gene_sets[[3]], "GO_<go_term>", "condition.stat")
# results %>% get_gene_set_results(list_of_gene_sets[[3]], "GO_<go_term>", "condition.pval") %>% head
