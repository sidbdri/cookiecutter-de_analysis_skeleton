SPECIES <- "unknown_species"
META_DATA=str_c('meta_data_',SPECIES,'.R')
source(META_DATA)

# Note that when comparisons are run in parallel in RStudio, the output is silent and the R session 
# will be hung until all sub-processes finish or are terminated. When running on the command line, 
# we will see the output but in a random order from each core. Thus, we might want to turn off 
# parallel when debugging in RStudio.
start_parallel(nrow(COMPARISON_TABLE))
#stop_parallel()

TX_LEVEL <- FALSE
ID_COLUMN <- ifelse(TX_LEVEL,'transcript','gene')
QUANT_METHOD <- 'salmon'
USE_TX <- TRUE

{% if cookiecutter.qSVA !="no" %}
qSVA <- TRUE
{% else %}
qSVA <- FALSE
{% endif %}

MISASSIGNMENT_PRECENTAGE = MISASSIGNMENT_SAMPLE_REFERENCE_TABLE %>% nrow() >0

PLOT_TO_FILE <- TRUE

OUTPUT_DIR <- 'results/differential_expression/de_tx'
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR,recursive=TRUE)

GRAPHS_DIR <- 'results/differential_expression/graphs_tx/'
if (!dir.exists(GRAPHS_DIR)) dir.create(GRAPHS_DIR, recursive=TRUE)

#####

total_dds_data <- get_total_dds_tximport(SAMPLE_DATA, QUANT_METHOD, TX_LEVEL)
total_vst <- total_dds_data %>% varianceStabilizingTransformation

start_plot("pca_all_tx")
total_vst %>% plot_pca(intgroup = PCA_FEATURE) %>% print()
end_plot()


num_features <- SAMPLE_DATA %>% dplyr::select(-species,-sample_name) %>% colnames() %>% length()
start_plot("pca_features_tx")
if(exists(x = 'patchworkplot')) rm(patchworkplot)
#This is to plot individually every feature defined in the SAMPLE_DATA table
SAMPLE_DATA %>% dplyr::select(-species,-sample_name) %>% colnames() %>%
    walk(function(feature){
        total_vst %>% plot_pca(intgroup=c(feature),FALSE) %>%
        add_to_patchwork(plot_var_name='patchworkplot')
})
patchworkplot
end_plot()


start_plot("heatmap_all_tx")
total_vst %>% plot_heat_map(
  SAMPLE_DATA %>% 
    tibble::rownames_to_column("tmp_sample_name") %>%
    tidyr::unite(col='sample_info', c(tmp_sample_name, HEAT_MAP_FEATURE), 
                 sep = ":", remove = FALSE) %>% 
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

## Result table contains transcript id, thus now suitable to use this for function check_cell_type
txi <- get_tximport(SAMPLE_DATA, QUANT_METHOD, TX_LEVEL)
tpms <- txi$abundance %>% 
  as.data.frame() %>%
  rename_all(funs(str_c(., "_tpm")))

if (TX_LEVEL) {
  #The gene column here is acturally transcript.
  tpms %<>% tibble::rownames_to_column('gene') 
  tpms %<>% left_join(get_avg_tpm(., TX_LEVEL)) %>%
    rename(transcript=gene)
  
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
    tpms %<>% tibble::rownames_to_column('gene') %>%
        inner_join(get_transcripts_to_genes(SPECIES)) %>%
        left_join(get_avg_tpm(tpms,TX_LEVEL))

    results %<>%
        left_join(tpms) %>%
        left_join(gene_info) %>%
        left_join(gene_lengths)
}



# run all get_res() functions and add to main "results" object
if(exists(x = 'all_comparison_pvalue_distribution')) rm(all_comparison_pvalue_distribution)
comparisons_results<-COMPARISON_TABLE %>% pull(comparison) %>% lapply_fork (
  function(comparison_name) {
    res <- get_res(comparison_name, tpms, use_tx=USE_TX, 
                   quant_method=QUANT_METHOD, tx_level=TX_LEVEL)
    
    results_tb <- get_global("results") %>%
      left_join(res[[1]]) %>%
      dplyr::rename(!!str_c(comparison_name, '.l2fc') := log2FoldChange,
                    !!str_c(comparison_name, '.raw_l2fc') := raw_l2fc,
                    !!str_c(comparison_name, '.stat') := stat,
                    !!str_c(comparison_name, '.pval') := pvalue,
                    !!str_c(comparison_name, '.padj') := padj)
    
    ##work out misassigned percentage
    if(MISASSIGNMENT_PRECENTAGE){
      # not tested 
      # P<-get_misassigned_percentage(comparison_name)
      # 
      # results_sargasso %<>% left_join(P$P_condition %>% dplyr::select(gene,!!str_c(comparison_name, '.perc.',COMPARISON_TABLE %>% filter(comparison==comparison_name) %>% pull(condition)) := p))
      # results_sargasso %<>% left_join(P$P_condition_base %>% dplyr::select(gene,!!str_c(comparison_name, '.perc.',COMPARISON_TABLE %>% filter(comparison==comparison_name) %>% pull(condition_base)) := p))
      # 
      # SUMMARY_TB <- get_global("SUMMARY_TB") %>%
      #   mutate(Misassignment_samples_in_comparison_level_condition = ifelse(Comparison == comparison_name,P$condition_reference_samples %>% str_c(collapse = ','),Misassignment_samples_in_comparison_level_condition)) %>%
      #   mutate(Misassignment_samples_in_base_level_condition = ifelse(Comparison == comparison_name,P$condition_base_reference_samples %>% str_c(collapse = ','),Misassignment_samples_in_base_level_condition)) %>% 
      #   set_global("SUMMARY_TB")
    }
    
    p_plot<-plot_pvalue_distribution(results_tb, str_c(comparison_name,'.pval'))
    
    ## return the results and merge them later
    list(comparison_name = comparison_name,
         res = res$res,
         dds = res$dds,
         results_tb = results_tb %>% dplyr::select(ID_COLUMN,contains('.')),
         summary_tb = res$summary_tb_row,
         p_plot = p_plot)
  }
)

if(exists(x = 'all_comparison_pvalue_distribution')) rm(all_comparison_pvalue_distribution)

lapply(comparisons_results,function(cmp){
  ## merge the cmp result table into global results table
  get_global("results") %>%
    left_join(cmp$results_tb) %>% 
    set_global("results")
  
  ## merge the cmp summary table into global SUMMARY_TABLE
  get_global("SUMMARY_TB") %>%
    rbind(cmp$summary_tb %>% as.data.frame(stringsAsFactors=FALSE)) %>%
    set_global("SUMMARY_TB")
  
  ## merge the p value plots
  add_to_patchwork(cmp$p_plot,plot_var_name='all_comparison_pvalue_distribution')
  
  ## export the res
  # cmp$res %>% set_global(cmp$comparison %>% str_c('res', sep = '_'))
  cmp$dds %>% set_global(cmp$comparison %>% str_c('dds', sep = '_'))
  
  'success'
}) %>% invisible()  ##so lappy will not print out useless merging message

start_plot("all_comparison_pvalue_distribution")
all_comparison_pvalue_distribution
end_plot()

# save results
columns_included <- ifelse(TX_LEVEL,
                            c('transcript', 'transcript_length', 'gene', 'number_of_transcript'),
                            c('gene','gene_length', 'max_transcript_length'))

if (COMPARISON_TABLE %>% pull(group) %>% unique() %>% length() > 1) {
  save_results_by_group(results,USE_TX)
}

results %>%
  dplyr::select(
    columns_included, gene_name, chromosome, description, entrez_id, gene_type, everything(),
    -dplyr::contains("_tpm"), -dplyr::ends_with(".stat")) %>%
  write_csv(str_c(OUTPUT_DIR, "/deseq2_results_count_", SPECIES, "_tx_", ID_COLUMN, "_", QUANT_METHOD, ".csv"),na = "")

results %>%
  dplyr::select(
    columns_included, gene_name, chromosome, description, entrez_id, gene_type,
    dplyr::contains("_tpm"),
    COMPARISON_TABLE %>%
      pull(comparison) %>%
      sapply(FUN = function(x) results %>% colnames() %>% str_which(str_c("^", x, sep =''))) %>% 
      unlist() %>%
      as.vector() %>%
      unique(),
    -dplyr::ends_with(".stat")
  ) %>%
  write_csv(str_c(OUTPUT_DIR, "/deseq2_results_tpm_", SPECIES, "_tx_",ID_COLUMN, "_", QUANT_METHOD, ".csv"),na = "")

SUMMARY_TB %>%
  write_csv(str_c(OUTPUT_DIR, "/de_summary_", SPECIES, "_tx_",ID_COLUMN, "_", QUANT_METHOD, ".csv"),na = "")



# library (knitr)
# sink(str_c("results/differential_expression/de_summary_",
#             SPECIES,"_tx_",tx_level_str,"_",QUANT_METHOD,".md"))
# kable(SUMMARY_TB, format = 'markdown')
# sink()

#####

# For each comparison: 
#   - for the GO/Reactome analyses, we are using all/up/down regulated genes,
#   - for GSA, we are using four gene set categories: "CURATED", "MOTIF", "GO" and "CELL_TYPE"
# Thus we need to reduce the number of comparisons we run in parallel to ensure we are not using more cores
# than specified. The total number of cores used after the following line will be 3 * getOption("mc.cores")
if(PARALLEL) adjust_parallel_cores()

##### GO and GSA analyses

if (!USE_TX | !TX_LEVEL) {
  expressed_genes <- get_total_dds(SAMPLE_DATA, SPECIES, filter_low_counts=TRUE) %>% 
    get_count_data()
  
  GO_results <- COMPARISON_TABLE %>% 
    pull(comparison) %>% 
    set_names(.) %>% 
    lapply_socket(X=., function(comparison_name) {
      p_str <- str_c(comparison_name, '.padj')
      l2fc_str <- str_c(comparison_name, '.l2fc')
    
    results <- get_global("results")
    
    lapply_socket(cores=3,X=c('','.up','.down'),function(cmp){
      if(cmp=='.up'){
        r <- results %>% filter(get(p_str) < P.ADJ.CUTOFF & get(l2fc_str) > 0)
      }else if((cmp=='.down')){
        r <- results %>% filter(get(p_str) < P.ADJ.CUTOFF & get(l2fc_str) < 0)
      }else{
        r <- results %>% filter(get(p_str) < P.ADJ.CUTOFF)
      }
      perform_go_analyses(r, expressed_genes, comparison_name, cmp, SPECIES, out_dir="results/differential_expression/go_tx/")
    }) %>% set_names(str_c(comparison_name,c('.all','.up','.down')))
  })
  
  ##### Gene set enrichment analysis
  
  gene_set_categories <- list("CURATED", "MOTIF", "GO", "CELL_TYPE")
  
  list_of_gene_sets <- gene_set_categories %>% lapply_fork(cores=length(gene_set_categories), X=., function(category,...) get_gene_sets(SPECIES, category))
  
  GS_results <- COMPARISON_TABLE %>% pull(comparison) %>% lapply_fork(X=., function(comparison_name,...) {
    dds <- str_c(comparison_name, 'dds', sep = '_') %>% get_global()
    
    camera_results <- list_of_gene_sets %>% 
      map(function(category_gene_sets) {
        get_camera_results(dds, category_gene_sets, gene_info)
      })
    
    camera_results %>% set_global(str_c(comparison_name, 'camera_results', sep = '_'))
    
    
    lapply_fork(cores = 3, X = seq(1:length(gene_set_categories)), function(category,...) {
      de_res <- results %>% dplyr::select(
        gene, gene_name, entrez_id, 
        starts_with(str_c(comparison_name, ".")), 
        -starts_with(str_c(comparison_name, ".stat")))  
      write_camera_results(
        gene_set_categories[[category]], list_of_gene_sets[[category]], 
        comparison_name, SPECIES,
        de_res, camera_results[[category]], out_dir="results/differential_expression/gsa_tx/")
    }) 
    camera_results
  })
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
        perform_pathway_enrichment(r, expressed_genes, comparison_name, cmp, SPECIES, out_dir="results/differential_expression/reactome_tx/")
      }
      ) %>% set_names(str_c(comparison_name,c('.all', '.up', '.down')))
    }
    )
}

# results %>% plot_gene_set(list_of_gene_sets[[3]], "GO_<go_term>", "condition.stat")
# results %>% get_gene_set_results(list_of_gene_sets[[3]], "GO_<go_term>", "condition.pval") %>% head

##### Saving and loading the workspace

# Save the objects in the workspace for future analysis

rws <- "results/Rworkspace/"
if (!dir.exists(rws)) {
  dir.create(rws,recursive=TRUE)
}

save(list = ls() %>% grep(x = ., pattern='comparisons_results', value = T,invert = T),
     file = str_c(rws,"diff_expr_tx.RData"))

# Load the save workspace to get all objects back for analysis

# load_rs_data()
