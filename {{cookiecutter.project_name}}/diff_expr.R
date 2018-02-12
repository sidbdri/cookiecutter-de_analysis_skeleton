source("load_packages.R")
source("common_functions.R")

SAMPLE_NAMES <- c(condition1, condition2, etc) %>%
  outer(c(rep1, rep2, etc), str_c, sep="todo") %>%
  t %>%
  as.vector

SAMPLE_DATA <- data.frame(
  condition=...,
  row.names=SAMPLE_NAMES
)

get_total_dds <- function(filter_low_counts=FALSE) {
  # Collate count data
  total_count_data <- SAMPLE_DATA %>% 
    row.names() %>% 
    map(read_counts) %>%
    purrr::reduce(inner_join) %>%
    remove_gene_column()
    
  get_deseq2_dataset(
    total_count_data, SAMPLE_DATA,
    filter_low_counts=filter_low_counts, design_formula=~condition)
}

get_total_dds_tximport <- function(quant_method) {
  total_sample_data <- data.frame(
    condition=c(),
    sample=c(),
    filename=c(),
    row.names=c("<SAMPLE1>", "<SAMPLE2>", etc)
  )

  quant_file <- get_transcript_quant_file(quant_method)
  quant_files <- str_c("results/", quant_method, "_quant/", total_sample_data$filename, "/", quant_file)
  names(quant_files) <- total_sample_data$filename

  txi <- tximport(quant_files, type=quant_method, tx2gene=get_transcripts_to_genes(), reader=read_tsv)

  total_dds <- DESeqDataSetFromTximport(txi, total_sample_data, ~condition)
  total_dds <- DESeq(total_dds)

  return(total_dds)
}

get_condition_res <- function() {
  sample_data <- SAMPLE_DATA %>% 
    filter_with_rownames(...)

  dds <- sample_data %>% 
    row.names() %>% 
    map(read_counts) %>% 
    purrr::reduce(inner_join) %>%
    remove_gene_column() %>% 
    get_deseq2_dataset(sample_data, design_formula=~condition)

  res <- dds %>% 
    get_deseq2_results("<condition>", "<cond2>", "<cond1>") %>% 
    left_join(dds %>% get_raw_l2fc(sample_data, condition="<cond2>"))

  list(res, dds %<>% varianceStabilizingTransformation)
}
  
get_condition_res_tximport <- function(quant_method) {
  sample_data <- data.frame(
    condition=c(),
    sample=c(),
    filename=c(),
    row.names=c("<SAMPLE1>", "<SAMPLE2>", etc)
  )

  quant_file <- get_transcript_quant_file(quant_method)
  quant_files <- str_c("results/", quant_method, "_quant/", sample_data$filename, "/", quant_file)
  names(quant_files) <- sample_data$filename

  txi <- tximport(quant_files, type=quant_method, tx2gene=get_transcripts_to_genes(), reader=read_tsv)

  dds <- DESeqDataSetFromTximport(txi, sample_data, ~condition)
  dds <- dds[rowSums(counts(dds)) > 1, ]
  dds <- DESeq(dds)

  results <- dds %>% get_deseq2_results("condition", "<cond2>", "<cond1>")
  
  l2fc <- get_count_data(dds) %>%
    mutate(l2fc=log2(() / ()) %>%
    dplyr::select(gene, l2fc)
    
  results %<>% left_join(l2fc)
  
  return(results)
}

#####

total_dds_data <- get_total_dds()
total_vst <- total_dds_data %>% varianceStabilizingTransformation
total_vst %>% plot_pca_with_labels(intgroup=c("condition"))
total_vst %>% plot_heat_map(SAMPLE_DATA %>% 
                              mutate(sample_info=str_c(condition, ..., sep=":")) %>% 
                              extract2("sample_info"))

plot_count_distribution(total_dds_data, norm=F)
plot_count_distribution(total_dds_data, norm=T)

#####

gene_info <- get_gene_info("{{cookiecutter.species}}")
gene_lengths <- read_csv("results/gene_lengths.csv")

results <- total_dds_data %>% get_count_data() 

fpkms <- results %>% 
  get_fpkms(gene_lengths, colnames(results) %>% tail(-1), "_fpkm")

fpkms %<>% mutate(
  <CONDITION1>_fpkm_avg = (<SAMPLE1>_fpkm + etc)/n,
  etc.
)

results %<>% 
  inner_join(fpkms) %>%
  inner_join(gene_info) %>% 
  inner_join(gene_lengths)

condition_res <- get_condition_res()

results %<>% 
  left_join(condition_res[[1]], by="gene") %>%
  dplyr::rename(condition.l2fc=log2FoldChange,
         condition.raw_l2fc=raw_l2fc,
         condition.stat=stat,
         condition.pval=pvalue,
         condition.padj=padj)

plot_pvalue_distribution(results, "condition.pval")

results %>% 
  dplyr::select(gene, gene_name, chromosome, description, entrez_id, gene_type,
                gene_length, max_transcript_length,
                everything(), -dplyr::contains("_fpkm"), -dplyr::ends_with(".stat")) %>%
  write_csv("results/differential_expression/deseq2_results.csv")

results %>% 
  dplyr::select(gene, gene_name, chromosome, description, entrez_id, gene_type,
                gene_length, max_transcript_length,
         dplyr::contains("_fpkm"), 
         starts_with(comparison), etc., 
         -dplyr::ends_with(".stat")) %>% 
  write_csv("results/differential_expression/deseq2_results_fpkm.csv")

##### GO analyses

expressed_genes <- get_total_dds(TRUE) %>% get_count_data()

results %>% 
  filter(padj < 0.1) %>%
  perform_go_analyses(expressed_genes, "<condition_comparison>")

##### Gene set enrichment analysis

gene_set_categories <- list("CURATED", "MOTIF", "GO")

gene_sets <- gene_set_categories %>% 
  map(function(x) get_gene_sets("{{cookiecutter.species}}", x))

condition_camera_results <- gene_sets %>% 
  map(function(x) get_camera_results(condition_res[[2]], x, gene_info, ~condition))

for (category in seq(1:length(gene_set_categories))) {
  write_camera_results(
    gene_set_categories[[category]], gene_sets[[category]], "condition",
    results, condition_camera_results[[category]])
}

# results %>% plot_gene_set(gene_sets[[3]], "GO_<go_term>", "condition.stat")
# results %>% get_gene_set_results(gene_sets[[3]], "GO_<go_term>", "condition.pval") %>% head

##### Salmon/tximport analysis

results_salmon <- get_total_dds_tximport("salmon") %>% 
  get_count_data()

salmon_tpms <- SAMPLE_NAMES %>%
  map(get_salmon_tpms) %>%
  purrr::reduce(inner_join)

results_salmon %<>% 
  inner_join(fpkms) %>%
  inner_join(gene_info) %>% 
  inner_join(gene_lengths)

results_salmon %<>% 
  left_join(get_condition_res_tximport("salmon", "quant.sf"), by="gene") %>%
  dplyr::rename(condition.l2fc=log2FoldChange,
                condition.raw_l2fc=l2fc,
                condition.pval=pvalue,
                condition.padj=padj)

results_salmon %>% 
  dplyr::select(gene, gene_name, chromosome, description, entrez_id, gene_type,
                gene_length, max_transcript_length,
                everything(), -dplyr::contains("_tpm")) %>%
  write_csv("results/differential_expression/deseq2_salmon_results.csv")
           
results_salmon %>% 
  dplyr::select(gene, gene_name, chromosome, description, entrez_id, gene_type,
                gene_length, max_transcript_length,
         dplyr::contains("_tpm"), 
         starts_with(condition), etc.) %>% 
  write_csv("results/differential_expression/deseq2_salmon_results_fpkm.csv")

##### Kallisto/tximport analysis

results_kallisto <- get_total_dds_tximport("kallisto") %>% 
  get_count_data() 

kallisto_tpms <- SAMPLE_NAMES %>%
  map(get_kallisto_tpms) %>%
  purrr::reduce(inner_join)

results_kallisto %<>% 
  inner_join(kallisto_tpms) %>%
  inner_join(gene_info) %>% 
  inner_join(gene_lengths)

results_kallisto %<>% 
  left_join(get_condition_res_tximport("kallisto", "abundance.tsv"), by="gene") %>%
  dplyr::rename(condition.l2fc=log2FoldChange,
                condition.raw_l2fc=l2fc,
                condition.pval=pvalue,
                condition.padj=padj)

results_kallisto %>% 
  dplyr::select(gene, gene_name, chromosome, description, entrez_id, gene_type,
                gene_length, max_transcript_length,
                everything(), -dplyr::contains("_tpm")) %>%
  write_csv("results/differential_expression/deseq2_kallisto_results.csv")
           
results_kallisto %>% 
  dplyr::select(gene, gene_name, chromosome, description, entrez_id, gene_type,
                gene_length, max_transcript_length,
         dplyr::contains("_tpm"), 
         starts_with(condition), etc.) %>% 
  write_csv("results/differential_expression/deseq2_kallisto_results_fpkm.csv")
