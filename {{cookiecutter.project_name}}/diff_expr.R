source("meta_data.R")

SPECIES <- "{{cookiecutter.species}}"

# Note that when comparisons are run in parallel in R studio, the output are slient and the Rsession is hang till all sub-proccess finishs or terminaled.
# When running in command line, we will see the output but in a random order from each core. Thus, we might want to turn off parallel when debugging in Rstudio.
start_parallel(nrow(COMPARISON_TABLE))
#stop_parallel()

{% if cookiecutter.qSVA !="no" %}
qSVA <- TRUE
{% else %}
qSVA <- FALSE
{% endif %}
PLOT_TO_FILE <- TRUE

MISASSIGNMENT_PERCENTAGE <- MISASSIGNMENT_SAMPLE_REFERENCE_TABLE %>% nrow() > 0

OUTPUT_DIR <- 'results/differential_expression/de_gene/'
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive=TRUE)

GRAPHS_DIR <- 'results/differential_expression/graphs/'
if (!dir.exists(GRAPHS_DIR)) dir.create(GRAPHS_DIR, recursive=TRUE)

#####

total_dds_data <- get_total_dds(SAMPLE_DATA, SPECIES, qSVA=qSVA)
total_vst <- total_dds_data %>% varianceStabilizingTransformation

start_plot("pca_all")
total_vst %>% plot_pca_with_labels(intgroup=PCA_FEATURE) %>% print()
end_plot()

num_features <- SAMPLE_DATA %>% dplyr::select(-species,-sample_name) %>% colnames() %>% length()
pdf_scale_factor <- 6
start_plot("pca_features")
if(exists(x = 'patchworkplot',where = .GlobalEnv)) rm(patchworkplot,envir=.GlobalEnv)
# This is to plot individually every feature defined in the SAMPLE_DATA table
SAMPLE_DATA %>% dplyr::select(-species,-sample_name) %>% colnames() %>%
  walk(function(feature){
    total_vst %>% plot_pca(intgroup=c(feature),FALSE) %>%
                  add_to_patchwork(plot_var_name='patchworkplot')
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

# We want to generate plots of the fpkm of the marker genes in the samples
check_cell_type(results, fpkm_check_cutoff=5, print_check_log=TRUE, print_fpkm_table=FALSE)

# run all get_res() functions in parallel
# for debugging, it may be worth calling stop_parallel(), because the mclapply has problem printing out stdout in rstudio.
# see http://dept.stat.lsa.umich.edu/~jerrick/courses/stat701/notes/parallel.html#forking-with-mclapply
comparisons_results<-COMPARISON_TABLE %>% pull(comparison) %>% lapplyFunc.Fork (
  function(comparison_name) {
    res <- get_res(comparison_name, fpkms, SPECIES, qSVA=qSVA)
    
    results <- get("results", envir = .GlobalEnv) %>% 
      left_join(res[[1]], by="gene") %>%
      dplyr::rename(!!str_c(comparison_name, '.l2fc') := log2FoldChange,
                    !!str_c(comparison_name, '.raw_l2fc') := raw_l2fc,
                    !!str_c(comparison_name, '.stat') := stat,
                    !!str_c(comparison_name, '.pval') := pvalue,
                    !!str_c(comparison_name, '.padj') := padj)

    if (MISASSIGNMENT_PERCENTAGE) {
      P <- get_misassignment_percentages(comparison_name, gene_lengths)
    
      if (!is.na(P$condition_reference_samples)) {
        results %<>% left_join(
          P$P_condition %>% 
            dplyr::select(gene, !!str_c(comparison_name, '.perc.',COMPARISON_TABLE %>% filter(comparison==comparison_name) %>% pull(condition)) := p))
      }   
    
      if (!is.na(P$condition_base_reference_samples)) {
        results %<>% left_join(
          P$P_condition_base %>% 
            dplyr::select(gene,!!str_c(comparison_name, '.perc.',COMPARISON_TABLE %>% filter(comparison==comparison_name) %>% pull(condition_base)) := p))
      }   
    
      SUMMARY_TB <- get("SUMMARY_TB", envir = .GlobalEnv) %>% 
        mutate(Misassignment_samples_in_comparison_level_condition = 
                 ifelse(Comparison == comparison_name, 
                        P$condition_reference_samples %>% str_c(collapse = ','),
                        Misassignment_samples_in_comparison_level_condition)) %>% 
        mutate(Misassignment_samples_in_base_level_condition = 
                 ifelse(Comparison == comparison_name,
                        P$condition_base_reference_samples %>% str_c(collapse = ','),
                        Misassignment_samples_in_base_level_condition))
                            
      assign("SUMMARY_TB", SUMMARY_TB,envir = .GlobalEnv)
    }    

    p_plot<-plot_pvalue_distribution(results, str_c(comparison_name,'.pval'))
    
    add_to_patchwork(p_plot,plot_var_name='all_comparison_pvalue_distribution')

    assign("results", results,envir = .GlobalEnv)
    
    comparison_name %>% str_c('res', sep = '_') %>% assign(res, envir = .GlobalEnv)

    ##we return the results and merge them later
    list(comparison_name=comparison_name,
         res=res,
         results_tb=get("results", envir = .GlobalEnv),
         summary_tb=get("SUMMARY_TB", envir = .GlobalEnv),
         p_plot=p_plot)

  }
)


if(exists(x = 'all_comparison_pvalue_distribution')) rm(all_comparison_pvalue_distribution)
lapply(comparisons_results,function(cmp){
  # merge the cmp result table into global results table
  assign("results",
         get("results", envir = .GlobalEnv) %>% left_join(cmp$results_tb %>% dplyr::select(gene,contains('.'))),
         envir = .GlobalEnv)
  
  # merge the cmp summary table into global SUMMARY_TABLE
  assign("SUMMARY_TB",
         get("SUMMARY_TB", envir = .GlobalEnv) %>% rbind(cmp$summary_tb),
         envir = .GlobalEnv)
  
  # merge the p value plots
  add_to_patchwork(cmp$p_plot,plot_var_name='all_comparison_pvalue_distribution')
  
  # export the res
  cmp$comparison %>% str_c('res', sep = '_') %>% assign(cmp$res, envir = .GlobalEnv)

  'success'
})


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

if(COMPARISON_TABLE %>% pull(group) %>% unique()%>% length() > 1) save_results_by_group(results)

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
      sapply(FUN = function(x) results %>% colnames() %>% str_which(str_c("^", x, sep =''))) %>% unlist() %>%
      as.vector() %>% 
      unique(), 
    -dplyr::ends_with(".stat")) %>%
  write_csv(str_c(OUTPUT_DIR, "/deseq2_results_fpkm_", SPECIES, ".csv"))

SUMMARY_TB %>%
  write_csv(str_c(OUTPUT_DIR, "/de_summary_", SPECIES, ".csv"))


#####
## For each comparison, 
##   for the GO/reactome analysis, we are running all/up/down regulated genes,
##   for the GSEA, we are running three categories ("CURATED", "MOTIF", "GO")
## Thus we need to reduce the number of comparison we analysis in parallel to ensure we are not using more cores than specified.
## The total number of cores used after the following line will be 3 * getOption("mc.cores")
if(PARALLEL) adjust_parallel_cores()

##### GO analyses
expressed_genes <- get_total_dds(SAMPLE_DATA, SPECIES, filter_low_counts=TRUE) %>% 
  get_count_data()

COMPARISON_TABLE %>% pull(comparison) %>% lapplyFunc.Socket(X=.,function(comparison_name) {
  p_str <- str_c(comparison_name, '.padj')
  l2fc_str <- str_c(comparison_name, '.l2fc')
  
  results <- get("results",envir = .GlobalEnv)
  
  lapplyFunc.Socket(cores=3,X=c('','.up','.down'),function(cmp){
    if(cmp=='.up'){
      r <- results %>% filter(get(p_str) < P.ADJ.CUTOFF & get(l2fc_str) > 0)
    }else if((cmp=='.down')){
      r <- results %>% filter(get(p_str) < P.ADJ.CUTOFF & get(l2fc_str) < 0)
    }else{
      r <- results %>% filter(get(p_str) < P.ADJ.CUTOFF)
    }
    perform_go_analyses(r, expressed_genes, comparison_name, cmp, SPECIES)
  })
})

##### Reactome pathway analysis

COMPARISON_TABLE %>% pull(comparison) %>% lapplyFunc.Socket(X=.,function(comparison_name) {
  p_str=str_c(comparison_name, 'padj', sep = '.')
  l2fc_str=str_c(comparison_name, 'l2fc', sep = '.')

  results <- get("results",envir = .GlobalEnv)

  lapplyFunc.Socket(cores=3,X=c('','.up','.down'),function(cmp){
    if(cmp=='.up'){
      r <- results %>% filter(get(p_str) < P.ADJ.CUTOFF  & get(l2fc_str) > 0)
    }else if((cmp=='.down')){
      r <- results %>% filter(get(p_str) < P.ADJ.CUTOFF  & get(l2fc_str) < 0)
    }else{
      r <- results %>% filter(get(p_str) < P.ADJ.CUTOFF)
    }
    perform_pathway_enrichment(r, expressed_genes, comparison_name, cmp, SPECIES)
  })
})

##### Gene set enrichment analysis

gene_set_categories <- list("CURATED", "MOTIF", "GO")

list_of_gene_sets <- gene_set_categories %>% lapplyFunc.Fork(cores=length(gene_set_categories), X=., function(category,...) get_gene_sets(SPECIES, category))

COMPARISON_TABLE %>% pull(comparison) %>% lapplyFunc.Fork(X=., function(comparison_name,...) {
  res <- str_c(comparison_name, 'res', sep = '_') %>% get(envir = .GlobalEnv)
  
  camera_results <- list_of_gene_sets %>% 
    map(function(category_gene_sets) {
      get_camera_results(res[[2]], category_gene_sets, gene_info)
    })
  
  assign(str_c(comparison_name, 'camera_results', sep = '_'), 
         camera_results, envir = .GlobalEnv)
  
  lapplyFunc.Fork(cores=3,X=seq(1:length(gene_set_categories)),function(category,...){
    de_res <- results %>% dplyr::select(
      gene, gene_name, entrez_id, 
      starts_with(str_c(comparison_name, ".")), 
      -starts_with(str_c(comparison_name, ".stat")))  
    write_camera_results(
      gene_set_categories[[category]], list_of_gene_sets[[category]], 
      comparison_name, SPECIES,
      de_res, camera_results[[category]])
  })
  'success'
})

# results %>% plot_gene_set(list_of_gene_sets[[3]], "GO_<go_term>", "condition.stat")
# results %>% get_gene_set_results(list_of_gene_sets[[3]], "GO_<go_term>", "condition.pval") %>% head


save.image(file="diff_expr.RData")