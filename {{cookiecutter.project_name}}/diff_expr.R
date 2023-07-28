source("load_packages.R")
source("utility_functions.R")
source("parallel.R")
source("qc.R")
source("deseq.R")
source("enrichment_tests.R")
source("gene_set_tests.R")
source("reactome.R")

SPECIES <- "unknown_species"
source(str_c('meta_data_', SPECIES, '.R'))

# to run jobs in serial, for debugging
#set_num_cores(1) 

P.ADJ.CUTOFF <- 0.05
PLOT_TO_FILE <- TRUE
MISASSIGNMENT_PERCENTAGE <- MISASSIGNMENT_SAMPLE_REFERENCE_TABLE %>% nrow() > 0

{% if cookiecutter.qSVA !="no" %}
qSVA <- TRUE
{% else %}
qSVA <- FALSE
{% endif %}

OUTPUT_DIR <- file.path('results','differential_expression')
DE_OUT_DIR=file.path(OUTPUT_DIR, "de_gene")
GRAPHS_DIR <- file.path(OUTPUT_DIR, "graphs")
LOG_DIR=file.path('results','logs','R','diff_expr')
dir.create(DE_OUT_DIR, recursive = TRUE)
dir.create(GRAPHS_DIR, recursive = TRUE)
dir.create(LOG_DIR, recursive = TRUE)

# https://github.com/sidbdri/cookiecutter-de_analysis_skeleton/issues/210
# add mitocondiral gene count percentage in SAMPLE_DATA for PCR plot later
gene_info <- get_gene_info(SPECIES)
mitochondrial_genes <- gene_info %>% filter(chromosome == 'MT') %>% pull(gene)
nuclear_encoded_mitochondrial_genes <- get_nuclear_encoded_mitochondrial_genes()
gene_list <- c(mitochondrial_genes, nuclear_encoded_mitochondrial_genes)
SAMPLE_DATA$mt <- read_mitochondrial_gene_percent(SAMPLE_DATA$sample_name, gene_list)

### Save sample data ###
SAMPLE_DATA %>% write_csv(file.path(OUTPUT_DIR, str_c("sample_data", ".csv")), na = "")

#### Differential gene expression ####

## Check read distribution for sample region
## For a large region or large number of samples, increase the bin_width to 10, 100 or 1000
## For multiple samples, the plot will order by the input samples order

# g <- check_sample_bam(samples=SAMPLE_DATA  %>% arrange(CultureType) %>% pull(sample_name),
#                       species=SPECIES,chr='2',start= 75505857,end=75534985,bin_width = 1000)
# plot(g)
# # g + scale_y_continuous(trans='log10')

## Create dds and vst objects containing all samples 

total_dds_data <- get_total_dds(SAMPLE_DATA, SPECIES, qSVA = qSVA, design_formula = ~1)
total_vst <- total_dds_data %>% varianceStabilizingTransformation(blind = TRUE)

## Main PCA plot of all samples

start_plot("pca_all_samples")
total_vst %>% plot_pca(intgroup = FEATURES_FOR_ALL_SAMPLES_PCA, 
                       output_data_table_path = file.path(GRAPHS_DIR, str_c('pca_all_samples_', SPECIES, '.csv'))) %>% print()
end_plot()

## PCA plots showing, individually, every feature defined in the SAMPLE_DATA table

# scale the pdf base on number of features to be plotted
num_features <- SAMPLE_DATA %>% 
  dplyr::select(-contains('species'), -contains('sample_name')) %>% 
  colnames() %>% 
  length()

start_plot("pca_all_features_all_samples",num_plots = num_features)

if (global_exists('patchworkplot')) {
  rm_global('patchworkplot')
}

SAMPLE_DATA %>% dplyr::select(-species,-sample_name) %>% colnames() %>%
  walk(function(feature) {
    total_vst %>%
      plot_pca(intgroup = c(feature), FALSE) %>%
      add_to_patchwork(plot_var_name = 'patchworkplot')
})

patchworkplot %>% print()
end_plot()

## Heatmap of all samples

start_plot("heatmap_all_samples")
total_vst %>% plot_heat_map(
  SAMPLE_DATA %>%
    tibble::rownames_to_column("tmp_sample_name") %>%
    tidyr::unite(col = 'sample_info', c(tmp_sample_name, FEATURES_FOR_ALL_SAMPLES_HEATMAP),
                 sep = ":", remove = FALSE) %>%
    extract2("sample_info"))
end_plot()

## Plots of distributions of normalised and raw counts

start_plot("count_distribution_norm")
plot_count_distribution(total_dds_data, norm = T)
end_plot()

start_plot("count_distribution")
plot_count_distribution(total_dds_data, norm = F)
end_plot()

## Plot percentage of mitochondrial RNA in each sample

start_plot('mitochondrial_gene_count')
plot_gene_percentage(counts(total_dds_data), 
                     list(MT_mito = mitochondrial_genes, NUC_mito = nuclear_encoded_mitochondrial_genes),
                     use_percentage = FALSE) %>% print
end_plot()

start_plot('mitochondrial_gene_percentage')
plot_gene_percentage(counts(total_dds_data), 
                     list(MT_mito = mitochondrial_genes, NUC_mito = nuclear_encoded_mitochondrial_genes),
                     use_percentage = TRUE) %>% print
end_plot()

## Initialise a D.E. results object containing gene information, counts and FPKMs

gene_lengths <- get_gene_lengths(SPECIES)
results <- total_dds_data %>% get_count_data()

fpkms <- results %>%
  get_fpkms(gene_lengths, colnames(results) %>% tail(-1), "_fpkm")

results %<>%
  left_join(fpkms) %>%
  left_join(gene_info) %>%
  left_join(gene_lengths)

## Generate plots of the FPKMs of marker genes in the samples

check_cell_type(results, fpkm_check_cutoff = 5, print_check_log = TRUE, print_fpkm_table = FALSE)

## Perform gene-level differential expression

comparisons_results <- COMPARISON_TABLE %>% 
  pull(comparison) %>%
  run_jobs_with_shared_memory(job_function = 'run_deseq',
                              job_source = 'deseq.R',
                              results_tbl = results,
                              fpkms = fpkms,
                              species = SPECIES,
                              qSVA = qSVA)

## Merge all D.E. results into the main results table

if (global_exists('all_comparison_pvalue_distribution')) {
  rm(all_comparison_pvalue_distribution)
}

lapply(comparisons_results, function(cmp) {
  # merge the comparison results tables into global results table
  get_global("results") %>%
    left_join(cmp$results_tb %>% dplyr::select(gene,contains('.'))) %>%
    set_global("results")

  # merge the comparison summary tables into global SUMMARY_TABLE
  get_global("SUMMARY_TB") %>%
    rbind(cmp$summary_tb %>% as.data.frame(stringsAsFactors = FALSE)) %>%
    set_global("SUMMARY_TB")

  # merge the p value plots
  add_to_patchwork(cmp$p_plot, plot_var_name = 'all_comparison_pvalue_distribution')

  # plot top DE gene for each comparison
  start_plot(str_c('top_de_genes_fpkm_',cmp$comparison_name))
  plot_genes_fpkm(results,results %>% arrange(across(str_c(cmp$comparison_name,'.padj',sep = ''))) %>% head(4) %>% pull(gene),print_fpkm_table = F)
  end_plot()
  
  # volcano plot for each comparison - top 5 up and down regulated genes are labelled
  start_plot(str_c('volcano_plot_', cmp$comparison_name))
  plot_volcano(cmp$results_tb, cmp$comparison_name)
  end_plot()

  # export the res and dds
  # commnet out for now.
  # see https://github.com/sidbdri/cookiecutter-de_analysis_skeleton/issues/83
  # cmp$res %>% set_global(cmp$comparison %>% str_c('res', sep = '_'))
  # cmp$dds %>% set_global(cmp$comparison %>% str_c('dds', sep = '_'))
  cmp$vst %>% set_global(cmp$comparison %>% str_c('vst', sep = '_'))
  cmp$misassignment_percentage %>% set_global(cmp$comparison %>% str_c('misassignment_percentage', sep = '_'))
  'success'
})

## Plot distributions of p-values for each comparison

start_plot("all_comparison_pvalue_distribution")
all_comparison_pvalue_distribution
end_plot()

## Create scatter plots of average FPKMs for each comparison

plot_scatter_fpkm(results)


## Plot gene expression heatmap for each comparison

for(comparison_name in COMPARISON_TABLE %>% pull(comparison)) {
  start_plot(str_c('gene_expression_heatmap_',comparison_name,sep = ''))
  plot_expression_heatmap(comparison_name = comparison_name,top = 50) %>%
    ComplexHeatmap::draw(column_title=comparison_name,column_title_gp=grid::gpar(fontsize=16))
  end_plot()
}


## Save results to file

if (COMPARISON_TABLE %>% pull(group) %>% unique() %>% length() > 1) {
  save_results_by_group(results)
}

results %>%
  dplyr::select(
    gene, gene_name, chromosome, description, entrez_id, gene_type,
    gene_length, max_transcript_length,
    everything(), -dplyr::contains("_fpkm"), -dplyr::ends_with(".stat")) %>%
  write_csv(file.path(DE_OUT_DIR, str_c("deseq2_results_count_", SPECIES, ".csv")), na = "")

results %>%
  dplyr::select(
    gene, gene_name, chromosome, description, entrez_id, gene_type,
    gene_length, max_transcript_length,
    dplyr::contains("_fpkm"),
    COMPARISON_TABLE %>%
      pull(comparison) %>%
      sapply(FUN = function(x) results %>% colnames() %>% str_which(str_c("^", x, sep = ''))) %>%
      unlist() %>%
      as.vector() %>%
      unique(),
    -dplyr::ends_with(".stat")) %>%
  write_csv(file.path(DE_OUT_DIR, str_c("deseq2_results_fpkm_", SPECIES, ".csv")), na = "")

SUMMARY_TB %>%
  write_csv(file.path(DE_OUT_DIR, str_c("de_summary_", SPECIES, ".csv")), na = "")

#### GO enrichment analysis ####

expressed_genes <- 
  get_total_dds(SAMPLE_DATA, SPECIES, filter_low_counts = TRUE) %>%
  get_count_data()

GO_results <- COMPARISON_TABLE %>%
  get_run_topgo_job_strings() %>% 
  run_jobs_with_separate_memory(job_function = 'run_topgo',
                                job_source = "enrichment_tests.R",
                                result_tbl = results,
                                p.adj.cutoff = P.ADJ.CUTOFF,
                                expressed_genes = expressed_genes,
                                species = SPECIES,
                                out_dir = file.path(OUTPUT_DIR, "enrichment_tests"))

#### Reactome pathway analysis ####

Reactome_results <- COMPARISON_TABLE %>%
  get_run_reactome_job_strings() %>%
  run_jobs_with_separate_memory(job_function = 'run_reactome',
                                job_source = "reactome.R",
                                result_tbl = results,
                                p.adj.cutoff = P.ADJ.CUTOFF,
                                expressed_genes = expressed_genes,
                                species = SPECIES,
                                out_dir = file.path(OUTPUT_DIR, "reactome"),
                                log_dir=LOG_DIR)

#### Gene set enrichment analysis ####

gene_set_categories <- list("CURATED", "MOTIF", "GO", "CELL_TYPE", "MSIGDB_CELL_TYPE")

list_of_gene_sets <- gene_set_categories %>%
  set_names(.) %>% 
  lapply(function(category, ...) get_gene_sets(SPECIES, category))

gene_set_analysis_results <- COMPARISON_TABLE %>%
  get_run_gene_set_analysis_job_strings(gene_set_categories) %>%
  run_jobs_with_separate_memory(job_function = 'run_gene_set_analysis',
                                job_source = 'gene_set_tests.R',
                                result_tbl = results,
                                comparison_table = COMPARISON_TABLE,
                                vst_list = str_c(COMPARISON_TABLE %>% pull(comparison), 'vst', sep = '_') %>%
                                  set_names(.) %>% 
                                  lapply(get_global),
                                list_of_gene_sets = list_of_gene_sets,
                                gene_info = gene_info,
                                species = SPECIES,
                                out_dir = file.path(OUTPUT_DIR, "gene_set_tests"),
                                log_dir = LOG_DIR)

gene_set_analysis_heatmaps <- gene_set_analysis_results %>% 
  names %>%
  run_jobs_with_shared_memory(job_function = 'plot_significant_set_heatmap',
                              results_tbl = results,
                              gs_result = gene_set_analysis_results,
                              list_of_gene_set = list_of_gene_sets,
                              comparison_tbl = COMPARISON_TABLE,
                              species = SPECIES,
                              out_dir = file.path(OUTPUT_DIR, "gene_set_tests")) 

#### Saving and loading the workspace ####

## Save the objects in the workspace for future analysis

rws <- "results/Rworkspace/"
if (!dir.exists(rws)) {
  dir.create(rws, recursive = TRUE)
}

  save(list = ls() %>% grep(x = ., pattern='comparisons_results', value = T,invert = T),
     file = str_c(rws, "diff_expr_", SPECIES,".RData"))

## Load the save workspace to get all objects back for analysis

# load_rs_data()
