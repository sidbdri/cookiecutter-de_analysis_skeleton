library(DESeq2)
library(dplyr)
library(gplots)
library(ggplot2)
library(magrittr)
library(RColorBrewer)
library(readr)
library(stringr)
library(tximport)

##### FUNCTIONS

read_counts <- function(sample) {
  counts_file_name <- str_c("results/read_counts/", sample, ".counts")
  counts <- read_tsv(counts_file_name, col_names=c("gene", str_c(sample)))
  return(counts)
}

remove_gene_column <- function(count_data) {
  row.names(count_data) <- count_data$gene
  count_data %<>% select(-gene)
  
  return(count_data)
}

get_deseq2_dataset <- function(count_data, sample_data, filter_low_counts=TRUE, 
                               design_formula=~condition) {
  
  dds <- DESeqDataSetFromMatrix(countData=count_data, colData=sample_data, design=design_formula)
  
  if (filter_low_counts) {
    dds <- dds[rowSums(counts(dds)) > 1, ]
  }
  
  return(DESeq(dds))
}

get_deseq2_results <- function(dds, comparison, condition, condition_base) {
  res <- results(dds, c(comparison, condition, condition_base))
  print(summary(res))
  
  res %<>% as.data.frame() %>%
    tibble::rownames_to_column(var="gene") %>%
    select(-baseMean, -lfcSE, -stat)
  
  return(res)
}

get_deseq2_results_name <- function(dds, name) {
  res <- results(dds, name=name)
  print(summary(res))
  
  res %<>% as.data.frame() %>%
    tibble::rownames_to_column(var="gene") %>%
    select(-baseMean, -lfcSE, -stat)
  
  return(res)
}

get_count_data <- function(dds, norm=T) {
  counts <- dds %>%
    counts(normalized=norm) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var="gene")
  
  return(counts)
}

plot_heat_map <- function(rld, sample_data) {
  distsRL <- dist(t(assay(rld)))
  mat <- as.matrix(distsRL)
  rownames(mat) <- colnames(mat) <- row.names(sample_data)
  hc <- hclust(distsRL)
  hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
  heatmap.2(mat, Rowv=as.dendrogram(hc), 
            symm=TRUE, trace="none",
            col = rev(hmcol), margin=c(10, 10))
}

get_gene_info <- function() {
  gene_info <- read_tsv(str_c("data/{{cookiecutter.species}}_ensembl-{{cookiecutter.ensembl_version}}/genes.tsv"),
                        col_names = c("gene", "description", "chromosome", "gene_name"),
                        col_types = list(chromosome = col_character()))
  
  return(gene_info)
}

get_fpkms <- function(all_counts, gene_lengths, samples, col_suffix) {
  all_counts %<>% inner_join(gene_lengths)
  
  for (sample in samples) {
    mmr <- sum(all_counts[[sample]]) / 1000000
    
    all_counts[,str_c(sample, col_suffix)] <- 
      map2_dbl(all_counts[[sample]], 
               all_counts[["max_transcript_length"]], 
               function(x, y) x / y / mmr * 1000)
  }
  
  all_counts %>% dplyr::select(gene, dplyr::contains(col_suffix))
}

#####

get_total_dds <- function() {
  # Collate count data
  total_count_data <- 
    read_counts("<SAMPLE1>") %>%
    inner_join(read_counts("<SAMPLE2>") %>%
    etc. %>%
    remove_gene_column()
    
  total_sample_data <- data.frame(
    condition=c(),
    sample=c(),
    row.names=c("<SAMPLE1>", "<SAMPLE2>", etc)
  )
  
  #total_sample_data$genotype %<>% relevel("WT")
  
  total_dds <- get_deseq2_dataset(
    total_count_data, total_sample_data,
    filter_low_counts=FALSE, design_formula=~genotype)
  
  return(total_dds)
}

get_total_dds_tximport <- function(quant_method, quant_file) {
  tx2gene <- read_delim("results/tx2gene.tsv", " ", col_names=c("transcript", "gene"))

  sample_data <- data.frame(
    condition=c(),
    sample=c(),
    filename=c(),
    row.names=c("<SAMPLE1>", "<SAMPLE2>", etc)
  )

  quant_files <- str_c("results/", quant_method, "_quant/", sample_data$filename, "/", quant_file)
  names(quant_files) <- sample_data$filename

  txi <- tximport(quant_files, type=quant_method, tx2gene=tx2gene, reader=read_tsv)

  total_dds <- DESeqDataSetFromTxImport(txi, sample_data, ~genotype)
  total_dds <- DESeq(dds)

  return(total_dds)
}

get_condition_res <- function() {
  count_data <- 
    read_counts("<SAMPLE1>") %>%
    inner_join(read_counts("<SAMPLE2>") %>%
    etc. %>%
    remove_gene_column()
    
  sample_data <- data.frame(
    condition=c(),
    sample=c(),
    row.names=c("<SAMPLE1>", "<SAMPLE2>", etc)
  )
  
  #sample_data$genotype %<>% relevel("WT")
  
  dds <- get_deseq2_dataset(
    count_data, sample_data, design_formula=~genotype)
  
  rld <- rlog(dds)
  
  plotPCA(rld, intgroup=c("genotype"))
  plot_heat_map(rld, sample_data)
  
  results <- get_deseq2_results(dds, "genotype", "KO", "WT")
  
  l2fc <- get_count_data(dds) %>%
    mutate(l2fc=log2(() / ()) %>%
    select(gene, l2fc)
    
  results %<>% left_join(l2fc)
  
  return(results)
}
  
get_condition_res_tximport <- function(quant_method, quant_file) {
  tx2gene <- read_delim("results/tx2gene.tsv", " ", col_names=c("transcript", "gene"))

  sample_data <- data.frame(
    condition=c(),
    sample=c(),
    filename=c(),
    row.names=c("<SAMPLE1>", "<SAMPLE2>", etc)
  )

  quant_files <- str_c("results/", quant_method, "_quant/", sample_data$filename, "/", quant_file)
  names(quant_files) <- sample_data$filename

  txi <- tximport(quant_files, type=quant_method, tx2gene=tx2gene, reader=read_tsv)

  dds <- DESeqDataSetFromTxImport(txi, sample_data, ~genotype)
  dds <- dds[rowSums(counts(dds)) > 1, ]
  dds <- DESeq(dds)

  results <- get_deseq2_results(dds, "genotype", "KO", "WT")
  
  l2fc <- get_count_data(dds) %>%
    mutate(l2fc=log2(() / ()) %>%
    select(gene, l2fc)
    
  results %<>% left_join(l2fc)
  
  return(results)

}

#####

gene_info <- get_gene_info()
gene_lengths <- read_csv("results/gene_lengths.csv")

results <- get_total_dds() %>% get_count_data() 

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

results %<>% 
  left_join(get_condition_res(), by="gene") %>%
  dplyr::rename(l2fc=log2FoldChange,
         raw_l2fc=l2fc,
         pval=pvalue,
         padj=padj)

results %>% 
  dplyr::select(gene, gene_name, chromosome, description, gene_length, max_transcript_length,
                everything(), -dplyr::contains("_fpkm")) %>%
  write_csv("results/differential_expression/deseq2_results.csv")

results %>% 
  dplyr::select(gene, gene_name, chromosome, description, gene_length, max_transcript_length,
         dplyr::contains("_fpkm"), 
         etc.)
  write_csv("results/differential_expression/deseq2_results_fpkm.csv")

#####

# TODO - look up quant file name
results_salmon <- get_total_dds_tximport("salmon", "quant.genes.sf") %>% 
  get_count_data() 

# TODO - read in accurate TPMs from quant files
#fpkms <- results %>% 
#  get_fpkms(gene_lengths, colnames(results) %>% tail(-1), "_fpkm")

results %<>% 
  inner_join(tpms) %>%
  inner_join(gene_info) %>% 
  inner_join(gene_lengths)

results %<>% 
  left_join(get_condition_res_tximport("salmon", "quant.genes.sf"), by="gene") %>%
  dplyr::rename(l2fc=log2FoldChange,
         raw_l2fc=l2fc,
         pval=pvalue,
         padj=padj)

results %>% 
  dplyr::select(gene, gene_name, chromosome, description, gene_length, max_transcript_length,
                everything(), -dplyr::contains("_tpm")) %>%
  write_csv("results/differential_expression/deseq2_salmon_results.csv")

#####

# TODO - look up quant file name
results_kallisto <- get_total_dds_tximport("kallisto", "abundance.tsv") %>% 
  get_count_data() 

# TODO - read in accurate TPMs from quant files
#fpkms <- results %>% 
#  get_fpkms(gene_lengths, colnames(results) %>% tail(-1), "_fpkm")

results %<>% 
  inner_join(tpms) %>%
  inner_join(gene_info) %>% 
  inner_join(gene_lengths)

results %<>% 
  left_join(get_condition_res_tximport("kallisto", "abundance.tsv"), by="gene") %>%
  dplyr::rename(l2fc=log2FoldChange,
         raw_l2fc=l2fc,
         pval=pvalue,
         padj=padj)

results %>% 
  dplyr::select(gene, gene_name, chromosome, description, gene_length, max_transcript_length,
                everything(), -dplyr::contains("_tpm")) %>%
  write_csv("results/differential_expression/deseq2_kallisto_results.csv")
