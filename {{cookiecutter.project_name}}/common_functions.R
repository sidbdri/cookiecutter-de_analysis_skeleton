KALLISTO <- "kallisto"
SALMON <- "salmon"

##### FUNCTIONS

filter_with_rownames <- function(.data, ...) {
  .data %>%
    tibble::rownames_to_column(var = "tmp_row_names") %>%
    filter_(.dots = lazyeval::lazy_dots(...)) %>%
    tibble::column_to_rownames(var = "tmp_row_names")
}

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
  
  dds %>% DESeq(betaPrior=TRUE)
}

get_deseq2_results <- function(dds, comparison, condition, condition_base) {
  res <- dds %>% results(c(comparison, condition, condition_base))
  res %>% summary %>% print
  
  res %>% as.data.frame() %>%
    tibble::rownames_to_column(var="gene") %>%
    dplyr::select(-baseMean, -lfcSE, -stat)
}

get_raw_l2fc <- function(dds, sample_data, ...) {
  sample_data %<>% tibble::rownames_to_column(var="tmp_sample_name")
  
  all_samples <- sample_data %>% extract2("tmp_sample_name") %>% as.vector
  
  comparison_sample_data <- filter_(sample_data, .dots=lazyeval::lazy_dots(...))
  comparison_samples <- comparison_sample_data %>% extract2("tmp_sample_name") %>% as.vector
  
  base_samples = all_samples %>% setdiff(comparison_samples)
  
  dds %>% 
    get_count_data %>% 
    mutate(comparison_mean=rowMeans(.[comparison_samples]),
           base_mean=rowMeans(.[base_samples]),
           raw_l2fc=log2(comparison_mean / base_mean)) %>%
    dplyr::select(gene, raw_l2fc)
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
    guides(fill=FALSE) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
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
    read_tsv(col_names = c("gene", "description", "chromosome", "gene_name", "entrez_id"),
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
  
  gene_info <- get_gene_info()
  go_results$Genes <- sapply(go_results[,c('GO.ID')], 
                             function(x) get_significant_genes(x, go_data, gene_info))
  
  go_results
}

perform_go_analyses <- function(significant_genes, expressed_genes, file_prefix) {
  if (significant_genes %>% nrow == 0) {
    stop("No significant genes supplied.")
  }
  
  c("BP", "MF", "CC") %>% walk(
    function(x) {
      perform_go_analysis(expressed_genes, significant_genes, x) %>%
        write_csv(str_c("results/differential_expression/go/{{cookiecutter.species}}_", file_prefix, "_go_", x %>% tolower, ".csv"))      
    }
  )
}

read_de_results <- function(filename, num_samples, num_conditions, num_comparisons, extra_columns="") {
  col_types_string <- str_c(
    "ccccii",
    strrep("d", num_samples + num_conditions + num_comparisons * 4),
    extra_columns)
  print(nchar(col_types_string))
  
  read_csv(filename, col_types=col_types_string)
}