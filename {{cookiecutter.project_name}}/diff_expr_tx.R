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

USE_TX <- TRUE
TX_LEVEL <- TRUE
QUANT_METHOD <- 'salmon'

# to run jobs in serial, for debugging
# set_num_cores(1)

P.ADJ.CUTOFF <- 0.05
PLOT_TO_FILE <- TRUE
MISASSIGNMENT_PERCENTAGE <- MISASSIGNMENT_SAMPLE_REFERENCE_TABLE %>% nrow() > 0

{% if cookiecutter.qSVA !="no" %}
qSVA <- TRUE
{% else %}
qSVA <- FALSE
{% endif %}

OUTPUT_DIR <- file.path('results','differential_expression_tx')
DE_OUT_DIR=file.path(OUTPUT_DIR, "de_tx")
GRAPHS_DIR <- file.path(OUTPUT_DIR, "graphs")
dir.create(DE_OUT_DIR, recursive = TRUE)
dir.create(GRAPHS_DIR, recursive = TRUE)

#### Differential gene expression ####

## Check read distribution for sample region
## For a large region or large number of samples, increase the bin_width to 10, 100 or 1000
## For multiple samples, the plot will order by the input samples order

# g <- check_sample_bam(samples=SAMPLE_DATA  %>% arrange(CultureType) %>% pull(sample_name),
#                       species=SPECIES,chr='2',start= 75505857,end=75534985,bin_width = 1000)
# plot(g)
# # g + scale_y_continuous(trans='log10')

## Create dds and vst objects containing all samples
total_dds_data <- get_total_dds_tximport(SAMPLE_DATA, SPECIES, QUANT_METHOD, TX_LEVEL,design_formula = ~1)
total_vst <- total_dds_data %>% varianceStabilizingTransformation()

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

gene_info <- get_gene_info(SPECIES)
mitochondrial_genes <- gene_info %>% filter(chromosome == 'MT') %>% pull(gene)
nuclear_encoded_mitochondrial_genes <- get_nuclear_encoded_mitochondrial_genes()

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

txi <- get_tximport(SAMPLE_DATA, SPECIES, QUANT_METHOD, TX_LEVEL)
tpms <- txi$abundance %>% as.data.frame() %>%
  rename_all(.funs=list(~str_c(., "_fpkm"))) %>%
  tibble::rownames_to_column('gene') %>%
  left_join(get_avg_tpm(.),by='gene')

if (TX_LEVEL) {
  #The gene column here is acturally transcript.
  tpms %<>% rename(transcript=gene)

  results %<>% rename(transcript=gene) %>%
    left_join(txi$Length %>% rename(transcript_length=length) %>%
                tibble::rownames_to_column('transcript'),by='transcript') %>%
    left_join(get_transcripts_to_genes(SPECIES),by='transcript') %>%
    left_join(tpms,by='transcript') %>%
    left_join(gene_info,by='gene') %>%
    left_join(gene_lengths,by='gene')

  # number of transcripts
  results %<>% left_join(
    results %>% dplyr::select(gene,transcript) %>%
      group_by(gene) %>% summarise(number_of_transcript = n()),by='gene')
} else {
  results %<>%
    left_join(tpms,by='gene') %>%
    left_join(gene_info,by='gene') %>%
    left_join(gene_lengths,by='gene')
}

## Generate plots of the FPKMs of marker genes in the samples
check_cell_type(results, fpkm_check_cutoff = 5, print_check_log = TRUE, print_fpkm_table = F)

## Perform gene-level differential expression
comparisons_results <- COMPARISON_TABLE %>%
  pull(comparison) %>%
  run_jobs_with_shared_memory(job_function = 'run_deseq',
                              job_source = 'deseq.R',
                              results_tbl = results,
                              fpkms = tpms,
                              species = SPECIES,
                              qSVA = qSVA,
                              use_tx=USE_TX,
                              quant_method=QUANT_METHOD,
                              tx_level=TX_LEVEL)

## Merge all D.E. results into the main results table
if (global_exists('all_comparison_pvalue_distribution')) {
  rm(all_comparison_pvalue_distribution)
}

lapply(comparisons_results, function(cmp) {
  # merge the comparison results tables into global results table
  get_global("results") %>%
    left_join(cmp$results_tb %>% dplyr::select(gene,any_of('transcript'),contains('.')),by='gene') %>%
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

  # export the res and dds
  # commnet out for now.
  # see https://github.com/sidbdri/cookiecutter-de_analysis_skeleton/issues/83
  # cmp$res %>% set_global(cmp$comparison %>% str_c('res', sep = '_'))
  cmp$dds %>% set_global(cmp$comparison %>% str_c('dds', sep = '_'))
  cmp$misassignment_percentage %>% set_global(cmp$comparison %>% str_c('misassignment_percentage', sep = '_'))
  'success'
})

## Plot distributions of p-values for each comparison

start_plot("all_comparison_pvalue_distribution")
all_comparison_pvalue_distribution
end_plot()

## Create scatter plots of average FPKMs for each comparison

plot_scatter_fpkm(results)

## Save results to file
## we change the _fpkm back to _tpm
results %<>% rename_at(.vars = vars(contains('_fpkm')),.funs=function(x){gsub(x = x,pattern = '_fpkm',replacement =  "_tpm")})


if (COMPARISON_TABLE %>% pull(group) %>% unique() %>% length() > 1) {
  save_results_by_group(results)
}

results %>%
  dplyr::select(
    any_of(c('transcript', 'transcript_length', 'gene', 'number_of_transcript','gene_length','max_transcript_length')),
    gene_name, chromosome, description, entrez_id, gene_type,
    everything(), -dplyr::contains("_tpm"), -dplyr::ends_with(".stat")) %>%
  write_csv(file.path(DE_OUT_DIR, str_c("deseq2_results_count_", SPECIES, ".csv")), na = "")


results %>%
  dplyr::select(
    any_of(c('transcript', 'transcript_length', 'gene', 'number_of_transcript','gene_length','max_transcript_length')),
    gene_name, chromosome, description, entrez_id, gene_type,
    everything(), -dplyr::contains("_tpm"), -dplyr::ends_with(".stat")) %>%
  write_csv(file.path(DE_OUT_DIR, str_c("deseq2_results_count_", SPECIES, ".csv")), na = "")


if(!TX_LEVEL){
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
                                  out_dir = file.path(OUTPUT_DIR, "reactome"))

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
                                  dds_list = str_c(COMPARISON_TABLE %>% pull(comparison), 'dds', sep = '_') %>%
                                    set_names(.) %>%
                                    lapply(get_global),
                                  list_of_gene_sets = list_of_gene_sets,
                                  gene_info = gene_info,
                                  species = SPECIES,
                                  out_dir = file.path(OUTPUT_DIR, "gene_set_tests"))

  gene_set_analysis_heatmaps <- gene_set_analysis_results %>%
    names %>%
    run_jobs_with_shared_memory(job_function = 'plot_significant_set_heatmap',
                                # we change the name back to _fpkm as the plot_significant_set_heatmap expect the result table to have fpkm column
                                results_tbl = results %>% rename_at(.vars = vars(contains('_tpm')),.funs=function(x){gsub(x = x,pattern = '_tpm',replacement =  "_fpkm")}),
                                gs_result = gene_set_analysis_results,
                                list_of_gene_set = list_of_gene_sets,
                                comparison_tbl = COMPARISON_TABLE,
                                species = SPECIES,
                                out_dir = file.path(OUTPUT_DIR, "gene_set_tests"))
}


#### Saving and loading the workspace ####

## Save the objects in the workspace for future analysis

rws <- "results/Rworkspace/"
if (!dir.exists(rws)) {
  dir.create(rws, recursive = TRUE)
}

save(list = ls() %>% grep(x = ., pattern='comparisons_results', value = T,invert = T),
     file = str_c(rws, "diff_expr_", SPECIES,"_transcript.RData"))

## Load the save workspace to get all objects back for analysis

# load_rs_data()
