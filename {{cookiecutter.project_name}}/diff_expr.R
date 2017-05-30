source("load_packages.R")

KALLISTO <- "kallisto"
SALMON <- "salmon"

##### FUNCTIONS

read_counts <- function(sample) {
  counts_file_name <- str_c("results/read_counts/", sample, ".counts")
  counts_file_name %>% read_tsv(col_names=c("gene", str_c(sample)))
}

remove_gene_column <- function(count_data) {
  row.names(count_data) <- count_data$gene
  count_data %>% dplyr::select(-gene)
}

get_deseq2_dataset <- function(count_data, sample_data, filter_low_counts=TRUE, 
                               design_formula=~condition) {
  
  dds <- DESeqDataSetFromMatrix(countData=count_data, colData=sample_data, design=design_formula)
  
  if (filter_low_counts) {
    dds <- dds[rowSums(counts(dds)) > 1, ]
  }
  
  dds %>% DESeq
}

get_deseq2_results <- function(dds, comparison, condition, condition_base) {
  res <- dds %>% results(c(comparison, condition, condition_base))
  res %>% summary %>% print
  
  res %>% as.data.frame() %>%
    tibble::rownames_to_column(var="gene") %>%
    dplyr::select(-baseMean, -lfcSE, -stat)
}

get_deseq2_results_name <- function(dds, name) {
  res <- dds %>% results(name=name)
  res %>% summary %>% print
  
  res %>% as.data.frame() %>%
    tibble::rownames_to_column(var="gene") %>%
    dplyr::select(-baseMean, -lfcSE, -stat)
}

get_count_data <- function(dds, norm=T) {
  dds %>%
    counts(normalized=norm) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var="gene")
}

plot_heat_map <- function(vst, sample_data) {
  distsRL <- vst %>% assay %>% t %>% dist

  mat <- distsRL %>% as.matrix()
  rownames(mat) <- colnames(mat) <- sample_data %>% row.names

  hc <- distsRL %>% hclust
  hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
  heatmap.2(mat, Rowv=hc %>% as.dendrogram, 
            symm=TRUE, trace="none",
            col = hmcol %>% rev, margin=c(10, 10))
}

plot_count_distribution <- function(dds, norm=T) {
  counts <- dds %>% 
    get_count_data(norm=norm) %>% 
    melt(id.vars=c("gene"), variable.name="sample", value.name="count") 
  
  p <- ggplot(counts, aes(sample, 1 + count)) + 
    geom_violin(aes(fill=sample), scale="width") + 
    geom_boxplot(width=.1, outlier.shape=NA) + 
    coord_trans(y = "log10") + 
    scale_y_continuous(breaks=c(1, 10, 100, 1000, 10000, 100000, 1000000)) +
    guides(fill=FALSE)
  
  print(p)
}

plot_pvalue_distribution <- function(results, pvalue_column) {
  pvals <- results %>% 
    filter_(str_c("!is.na(", pvalue_column, ")")) %>% 
    dplyr::select_(pvalue_column)
  
  p <- ggplot(pvals, aes_string(pvalue_column)) + 
    geom_histogram(binwidth=0.025) 
  
  print(p)
}

get_gene_info <- function() {
  "data/{{cookiecutter.species}}_ensembl_{{cookiecutter.ensembl_version}}/genes.tsv" %>% 
    read_tsv(col_names = c("gene", "description", "chromosome", "gene_name"),
             col_types = list(chromosome = col_character())) %>% 
    group_by(gene) %>% 
    filter(row_number()==1) %>% 
    ungroup
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

get_transcripts_to_genes <- function() {
  read_delim("results/tx2gene.tsv", " ", col_names=c("transcript", "gene"))
}

get_transcript_quant_file <- function(quant_method) {
  if (quant_method == SALMON) {
    "quant.sf"
  } else if (quant_method == KALLISTO) {
    "abundance.tsv" 
  } 
}

get_salmon_tpms <- function(sample) {
  sample %>% 
    str_c("results/salmon_quant/", ., "/quant.genes.sf") %>% 
    read_tsv %>% 
    dplyr::select(Name, TPM) %>% 
    rename_(.dots = setNames(names(.), c("gene", str_c(sample, "_tpm"))))
}

get_kallisto_tpms <- function(sample) {
  transcript_tpms <- sample %>% 
    str_c("results/kallisto_quant/", ., "/abundance.tsv") %>% 
    read_tsv %>% 
    dplyr::select(target_id, tpm)
  
  transcript_tpms %<>% inner_join(get_transcripts_to_genes(), by=c("target_id"="transcript")) 
  
  transcript_tpms %>% 
    group_by(gene) %>% 
    summarise(tpm = sum(tpm)) %>%
    rename_(.dots = setNames(names(.), c("gene", str_c(sample, "_tpm"))))
}

get_significant_genes <- function(term, GOdata, gene_info) {
  genes_for_term <- GOdata %>% genesInTerm(term) %>% extract2(1)
  significant_genes <- GOdata %>% sigGenes
  significant_genes_for_term <- genes_for_term %>% intersect(significant_genes)
  
  gene_info %>% 
    filter(gene %in% significant_genes_for_term) %>% 
    dplyr::select(gene_name) %>% 
    extract2(1) %>% 
    paste(collapse=", ")
}

perform_go_analysis <- function(gene_universe, significant_genes, ontology="BP") {
  gene_list <- (gene_universe$gene %in% significant_genes$gene) %>% as.integer %>% factor
  names(gene_list) <- gene_universe$gene
  
  mapping <- switch("{{cookiecutter.species}}",
         mouse = "org.Mm.eg.db",
         rat = "org.Rn.eg.db",
         human = "org.Hs.eg.db")
  
  go_data <- new("topGOdata", ontology=ontology, allGenes=gene_list,
                 annot=annFUN.org, mapping=mapping, ID="Ensembl")
  
  result_fisher <- go_data %>% runTest(algorithm="weight01", statistic="fisher")
  result_fisher %>% print
  
  go_results <- go_data %>% GenTable(weight_fisher=result_fisher, orderBy="weight_fisher", topNodes=150)
  
  go_results$Genes <- sapply(go_results[,c('GO.ID')], 
                             function(x) get_significant_genes(x, go_data, get_gene_info()))
  
  go_results
}

perform_go_analyses <- function(significant_genes, expressed_genes, file_prefix) {
  c("BP", "MF", "CC") %>% walk(
    function(x) {
      perform_go_analysis(expressed_genes, significant_genes, x) %>%
        write_csv(str_c("results/differential_expression/go/{{cookiecutter.species}}_", file_prefix, "_go_", x %>% tolower, ".csv"))      
    }
  )
}

#####

SAMPLE_NAMES <- c(condition1, condition2, etc) %>%
  outer(c(rep1, rep2, etc), str_c, sep="todo") %>%
  t %>%
  as.vector

get_total_dds <- function() {
  # Collate count data
  total_count_data <- SAMPLE_NAMES %>%
    map(read_counts) %>%
    reduce(inner_join) %>%
    remove_gene_column()
    
  total_sample_data <- data.frame(
    condition=c(),
    sample=c(),
    row.names=c("<SAMPLE1>", "<SAMPLE2>", etc)
  )
  
  #total_sample_data$genotype %<>% relevel("WT")
  
  total_dds <- get_deseq2_dataset(
    total_count_data, total_sample_data,
    filter_low_counts=FALSE, design_formula=~condition)
  
  list(total_sample_data, total_dds)
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
  count_data <- SAMPLE_NAMES %>%
    map(read_counts) %>%
    reduce(inner_join) %>%
    remove_gene_column()
    
  sample_data <- data.frame(
    condition=c(),
    sample=c(),
    row.names=c("<SAMPLE1>", "<SAMPLE2>", etc)
  )
  
  #sample_data$genotype %<>% relevel("WT")
  
  dds <- get_deseq2_dataset(
    count_data, sample_data, design_formula=~condition)
  
  vst <- dds %>% varianceStabilizingTransformation
  vst %>% plotPCA(intgroup=c("condition")) %>% print
  vst %>% plot_heat_map(sample_data)
  
  results <- dds %>% get_deseq2_results("condition", "<cond2>", "<cond1>")
  
  l2fc <- dds %>% get_count_data %>%
    mutate(l2fc=log2(() / ())) %>%
    dplyr::select(gene, l2fc)
    
  results %>% left_join(l2fc)
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
total_vst <- total_dds_data %>% extract2(2) %>% varianceStabilizingTransformation
total_vst %>% plotPCA(intgroup=c("condition"))
total_vst %>% plot_heat_map(total_dds_data %>% extract2(1))

plot_count_distribution(total_dds_data %>% extract2(2), norm=F)
plot_count_distribution(total_dds_data %>% extract2(2), norm=T)

#####

gene_info <- get_gene_info()
gene_lengths <- read_csv("results/gene_lengths.csv")

results <- total_dds_data %>% extract2(2) %>% get_count_data() 

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
  dplyr::rename(condition.l2fc=log2FoldChange,
         condition.raw_l2fc=l2fc,
         condition.pval=pvalue,
         condition.padj=padj)

plot_pvalue_distribution(results, "condition.pval")

results %>% 
  dplyr::select(gene, gene_name, chromosome, description, gene_length, max_transcript_length,
                everything(), -dplyr::contains("_fpkm")) %>%
  write_csv("results/differential_expression/deseq2_results.csv")

results %>% 
  dplyr::select(gene, gene_name, chromosome, description, gene_length, max_transcript_length,
         dplyr::contains("_fpkm"), 
         starts_with(condition), etc.) %>% 
  write_csv("results/differential_expression/deseq2_results_fpkm.csv")

##### GO analyses

results %>% 
  filter(padj < 0.1) %>%
  perform_go_analyses(expressed_genes, "<condition_comparison>")

##### Salmon/tximport analysis

results_salmon <- get_total_dds_tximport("salmon") %>% 
  get_count_data()

salmon_tpms <- SAMPLE_NAMES %>%
  map(get_salmon_tpms) %>%
  reduce(inner_join)

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
  dplyr::select(gene, gene_name, chromosome, description, gene_length, max_transcript_length,
                everything(), -dplyr::contains("_tpm")) %>%
  write_csv("results/differential_expression/deseq2_salmon_results.csv")
           
results_salmon %>% 
  dplyr::select(gene, gene_name, chromosome, description, gene_length, max_transcript_length,
         dplyr::contains("_tpm"), 
         starts_with(condition), etc.) %>% 
  write_csv("results/differential_expression/deseq2_salmon_results_fpkm.csv")

##### Kallisto/tximport analysis

results_kallisto <- get_total_dds_tximport("kallisto") %>% 
  get_count_data() 

kallisto_tpms <- SAMPLE_NAMES %>%
  map(get_kallisto_tpms) %>%
  reduce(inner_join)

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
  dplyr::select(gene, gene_name, chromosome, description, gene_length, max_transcript_length,
                everything(), -dplyr::contains("_tpm")) %>%
  write_csv("results/differential_expression/deseq2_kallisto_results.csv")
           
results_kallisto %>% 
  dplyr::select(gene, gene_name, chromosome, description, gene_length, max_transcript_length,
         dplyr::contains("_tpm"), 
         starts_with(condition), etc.) %>% 
  write_csv("results/differential_expression/deseq2_kallisto_results_fpkm.csv")
