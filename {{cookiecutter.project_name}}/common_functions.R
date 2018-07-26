KALLISTO <- "kallisto"
SALMON <- "salmon"

##### FUNCTIONS

filter_with_rownames <- function(.data, ...) {
  .data %>%
    tibble::rownames_to_column(var = "tmp_row_names") %>%
    filter_(.dots = lazyeval::lazy_dots(...)) %>%
    tibble::column_to_rownames(var = "tmp_row_names")
}

read_counts <- function(sample, species) {
  counts_file_name <- str_c("results/read_counts/", sample, ".", species,".counts")
  counts_file_name %>% read_tsv(col_names=c("gene", str_c(sample)))
}

remove_gene_column <- function(count_data) {
  row.names(count_data) <- count_data$gene
  count_data %>% dplyr::select(-gene)
}

get_deseq2_dataset <- function(count_data, sample_data, filter_low_counts=TRUE, 
                               design_formula=~condition,qSVA=FALSE) {
  
  dds <- DESeqDataSetFromMatrix(countData=count_data, colData=sample_data, design=design_formula)
  
  if (filter_low_counts) {
    dds <- dds[rowSums(counts(dds)) > 1, ]
  }

  if(qSVA) {
    dds %<>% get_qsva_dds()
  }
  
  dds %>% DESeq(betaPrior=TRUE)
}

get_deseq2_results <- function(dds, comparison, condition, condition_base, alpha=0.05) {
  res <- dds %>% results(c(comparison, condition, condition_base), alpha=alpha)
  res %>% summary %>% print
  
  res %>% as.data.frame() %>%
    tibble::rownames_to_column(var="gene") %>%
    dplyr::select(-baseMean, -lfcSE)
}

get_raw_l2fc <- function(dds, sample_data, comparison_samples_filter) {
  sample_data %<>% tibble::rownames_to_column(var="tmp_sample_name")
  
  all_samples <- sample_data %>% extract2("tmp_sample_name") %>% as.vector

  comparison_sample_data <- filter(sample_data, !!comparison_samples_filter)
  comparison_samples <- comparison_sample_data %>% extract2("tmp_sample_name") %>% as.vector

  base_samples = all_samples %>% setdiff(comparison_samples)
  
  dds %>% 
    get_count_data %>% 
    mutate(comparison_mean=rowMeans(.[comparison_samples]),
           base_mean=rowMeans(.[base_samples]),
           raw_l2fc=log2(comparison_mean / base_mean)) %>%
    dplyr::select(gene, raw_l2fc)
}

get_deseq2_results_name <- function(dds, name, alpha=0.05) {
  res <- dds %>% results(name=name, alpha=alpha)
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

plot_pca_with_labels <- function(vst, intgroup=c("condition")) {
  pca_data <- vst %>% plotPCA(intgroup=intgroup, returnData=TRUE)
  
  percent_var <- round(100 * attr(pca_data, "percentVar"))
  
  intgroup.df <- as.data.frame(colData(vst)[, intgroup, drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(vst)[[intgroup]]
  } 
  
  pca_data %>% 
    ggplot(aes(PC1, PC2, color=group)) +
    geom_point(size=3) +
    geom_text(aes(label = name), 
              colour="darkgrey", 
              position=position_nudge(y = 1), size=3) + 
    xlab(str_c("PC1: ", percent_var[1], "% variance")) +
    ylab(str_c("PC2: ", percent_var[2], "% variance")) + 
    theme(legend.position="none")
}

plot_heat_map <- function(vst, sample_names) {
  distsRL <- vst %>% assay %>% t %>% dist
  
  mat <- distsRL %>% as.matrix()
  rownames(mat) <- colnames(mat) <- sample_names
  
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

get_gene_info <- function(species) {
  species %>%
    str_c("data/", ., "_ensembl_{{cookiecutter.ensembl_version}}/genes.tsv") %>% 
    read_tsv(col_names = c("gene", "description", "chromosome", "gene_name", "entrez_id", "gene_type"),
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

  fpkms<-all_counts %>% dplyr::select(gene, dplyr::contains(col_suffix))
  # workout average fpkm
  fpkms %>% left_join(get_avg_fpkm(fpkms))
}

read_de_results <- function(filename, num_samples, num_conditions, num_comparisons, extra_columns="") {
  col_types_string <- str_c(
    "ccccii",
    strrep("d", num_samples + num_conditions + num_comparisons * 4),
    extra_columns)
  print(nchar(col_types_string))
  
  read_csv(filename, col_types=col_types_string)
}

##### Transcript-level D. E. analyses

get_transcripts_to_genes <- function(species='human') {
  str_c("data/", species , "_ensembl_{{cookiecutter.ensembl_version}}/tx2gene.tsv") %>%
    read_delim(delim = " ", col_names=c("transcript", "gene"))
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

##### Gene ontology enrichment analyses

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

perform_go_analysis <- function(gene_universe, significant_genes, ontology="BP", species) {
  gene_list <- (gene_universe$gene %in% significant_genes$gene) %>% as.integer %>% factor
  names(gene_list) <- gene_universe$gene
  
  mapping <- switch(species,
                    mouse = "org.Mm.eg.db",
                    rat = "org.Rn.eg.db",
                    human = "org.Hs.eg.db")
  
  go_data <- new("topGOdata", ontology=ontology, allGenes=gene_list,
                 annot=annFUN.org, mapping=mapping, ID="Ensembl")
  
  result_fisher <- go_data %>% runTest(algorithm="weight01", statistic="fisher")
  result_fisher %>% print
  
  go_results <- go_data %>% GenTable(weight_fisher=result_fisher, orderBy="weight_fisher", topNodes=150)
  
  gene_info <- get_gene_info(species)
  go_results$Genes <- sapply(go_results[,c('GO.ID')], 
                             function(x) get_significant_genes(x, go_data, gene_info))
  
  go_results
}

perform_go_analyses <- function(significant_genes, expressed_genes, file_prefix, species) {
  if (significant_genes %>% nrow == 0) {
    message("No significant genes supplied.")
    return()
  }
  
  top_dir<-str_c("results/differential_expression/go/",species,sep = '')
  if (!dir.exists(top_dir)) {
    dir.create(top_dir,recursive=TRUE)
  }

  c("BP", "MF", "CC") %>% walk(
    function(x) {
    perform_go_analysis(expressed_genes, significant_genes, x, species) %>%
    write_csv(str_c("results/differential_expression/go/",species,"/" ,file_prefix, "_go_", x %>% tolower, ".csv"))
    }
  )
}

##### Camera gene set enrichment analysis

get_human_vs_species_ortholog_info <- function(species) {
  orthologs <- species %>% 
    str_c("data/", ., "_ensembl_{{cookiecutter.ensembl_version}}/human_orthologs.tsv") %>% 
    read_tsv(col_names=c("species_gene", "human_gene", "type"))
  
  species_genes <- get_gene_info(species)
  human_genes <- get_gene_info("human")
  
  orthologs_entrez <- orthologs %>% 
    left_join(species_genes %>% 
                dplyr::select(gene, entrez_id) %>% 
                dplyr::rename(species_gene=gene, species_entrez_id=entrez_id)) %>%
    left_join(human_genes %>% 
                dplyr::select(gene, entrez_id) %>% 
                dplyr::rename(human_gene=gene, human_entrez_id=entrez_id)) %>%
    dplyr::select(species_entrez_id, human_entrez_id) %>% 
    filter(!is.na(species_entrez_id) & !is.na(human_entrez_id)) %>% 
    distinct()
  
  same_name_entrez <- (species_genes %>% 
                         mutate(gene_name=tolower(gene_name)) %>% 
                         dplyr::select(gene_name, entrez_id) %>% 
                         filter(!is.na(entrez_id)) %>% 
                         dplyr::rename(species_entrez_id=entrez_id)) %>%
    inner_join(human_genes %>% 
                 mutate(gene_name=tolower(gene_name)) %>% 
                 dplyr::select(gene_name, entrez_id) %>% 
                 filter(!is.na(entrez_id)) %>% 
                 dplyr::rename(human_entrez_id=entrez_id)) %>%
    dplyr::select(-gene_name) %>% 
    distinct()
  
  human_species_entrez_mappings <- orthologs_entrez %>% 
    rbind(same_name_entrez) %>% 
    distinct
}

get_gene_sets <- function(species, gene_set_name) {
  msigdb_data <- str_c("data/msigdb/v5.2/", gene_set_name, ".all.v5.2.entrez.gmt.txt") %>% 
    read_lines()
  
  gene_set_names <- map_chr(msigdb_data, function(x) {str_split(x, "\t")[[1]][1]})
  gene_sets <- map(msigdb_data, function(x) {str_split(str_split(x, "\t", n=3)[[1]][3], "\t")[[1]]})
  
  names(gene_sets) <- gene_set_names
  
  if (species == "human") {
    return(gene_sets)
  }
  
  pb <- txtProgressBar(max=length(gene_sets), style=3)
  count <- 0
  
  ortholog_info <- get_human_vs_species_ortholog_info(species)
  
  ret <- gene_sets %>% map(function(y) {
    count <<- count + 1
    setTxtProgressBar(pb, count)
    ortholog_info[which(ortholog_info$human_entrez_id %in% y),]$species_entrez_id %>% unique
  }) 
  
  close(pb)
  
  ret
}

get_camera_results <- function(dds, gene_sets, gene_info) {
  vst <- dds %>% varianceStabilizingTransformation
  design_formula <- dds %>% design()
  expression_data <- vst %>% assay
  
  ids <- expression_data %>% 
    as.data.frame %>% 
    tibble::rownames_to_column(var="gene") %>% 
    inner_join(gene_info) 
  
  idx <- gene_sets %>% ids2indices(id=ids$entrez_id)
  
  design_matrix <- model.matrix(design_formula, vst %>% colData)
  
  expression_data %>% camera(idx, design_matrix)
}

plot_gene_set <- function(results, gene_sets, gene_set_name, prefix) {
  idx <- gene_sets %>% 
    ids2indices(id=results$entrez_id) %>%
    extract2(gene_set_name)
  
  pval <- prefix %>% str_c(".pval") %>% rlang::sym()
  l2fc <- prefix %>% str_c(".l2fc") %>% rlang::sym()
  
  results %>% 
    mutate(signed_p = -log10(!!pval) * sign(!!l2fc)) %>%
    pull(signed_p) %>% 
    barcodeplot(index=idx, quantiles = c(-1,1)*(-log10(0.05)))
}

get_gene_set_results <- function(results, gene_sets, gene_set_name, pvalue) {
  idx <- gene_sets %>% 
    extract2(gene_set_name) %>% 
    match(results$entrez_id) %>% 
    na.omit
  
    results %>% magrittr::extract(idx, ) %>% arrange_(pvalue)
}

write_camera_results <- function(
gene_set_collection_name, gene_set_collection, comparison_name, species, de_results, camera_results,
  barcodeplots=FALSE) {
  
  camera_results %<>% 
    tibble::rownames_to_column(var="GeneSet") %>% 
    filter(FDR < 0.1)
  
  if ((camera_results %>% nrow) == 0) {
    return()
  }
  
  top_dir <- str_c("results/differential_expression/gsa/",species, "/", comparison_name)
  if (!dir.exists(top_dir)) {
    dir.create(top_dir,recursive=TRUE)
  }
  
  camera_results %>% write_csv(str_c(top_dir, "/", gene_set_collection_name, "_enriched_sets.csv"))
  
  sub_dir <- str_c(top_dir, "/", gene_set_collection_name)
  if (!dir.exists(sub_dir)) {
    dir.create(sub_dir,recursive=TRUE)
  }
  
  camera_results %>% 
    extract2("GeneSet") %>% 
    walk(function(x) {
      gene_set_results <- de_results %>% get_gene_set_results(gene_set_collection, x, str_c(comparison_name, ".pval"))
      gene_set_results %>% write_csv(str_c(sub_dir, "/", x, ".csv"))
      
      if (barcodeplots) {
        pdf(str_c(sub_dir, "/", x, ".pdf"))
        plot_gene_set(de_results, gene_set_collection, x, comparison_name)
        dev.off()
      }
    })
}

get_avg_fpkm <- function(fpkms){
  sample_data = SAMPLE_DATA
  sample_data %<>% group_by(.dots=AVG_FPKM_GROUP) %>%
    summarise(samples=str_c(sample_name,'_fpkm',sep = '',collapse = ',')) %>%
    tidyr::unite('avg_name',AVG_FPKM_GROUP,sep='_')

  for (avg in sample_data %>% pull(avg_name)){
    samples <- sample_data %>%
      filter(avg_name == avg) %>%
      pull(samples) %>%
      str_split(',') %>% unlist()

    avg_fpkm <- fpkms %>% dplyr::select(one_of(samples)) %>%
      mutate(avg_fpkm=rowMeans(.)) %>%
      dplyr::pull(avg_fpkm)

    fpkms %<>% mutate(!!str_c(avg,'_avg_fpkm',sep='') := avg_fpkm)
  }

  fpkms %>% dplyr::select(gene,contains('avg'))
}



get_avg_tpm<-function(sample_data,tpms,comparision){
  sample_data = SAMPLE_DATA
  sample_data %<>% group_by(.dots=AVG_FPKM_GROUP) %>%
    summarise(samples=str_c(sample_name,'_tpm',sep = '',collapse = ',')) %>%
    tidyr::unite('avg_name',AVG_FPKM_GROUP,sep='_')

  for (avg in sample_data %>% pull(avg_name)){
    samples <- sample_data %>%
      filter(avg_name == avg) %>%
      pull(samples) %>%
      str_split(',') %>% unlist()

    avg_tpm <- tpms %>% dplyr::select(one_of(samples)) %>%
      mutate(avg_tpm=rowMeans(.)) %>%
      dplyr::pull(avg_tpm)

    tpms %<>% mutate(!!str_c(avg,'_avg_tpm',sep='') := avg_tpm)
  }

  id_column=ifelse(tx_level,'transcript','gene')


  tpms %>% dplyr::select(!!id_column,contains('avg'))
}

get_quality_surrogate_variables <- function(dds) {
  design_formula = dds %>% design()
  sample_data<- dds %>% colData()
  sample_names <- sample_data %>% rownames()

  quality_surrogate_variables <- sample_names %>%
    str_c("results/read_counts/", ., ".dm.tsv") %>%
    read.degradation.matrix(
      sampleNames=sample_names,
      totalMapped=dds %>% counts %>% colSums,
      readLength=75,
      type="region_matrix_single",
      BPPARAM=SerialParam()) %>%
    qsva(mod=design_formula %>% model.matrix(sample_data))

  quality_surrogate_variables %>% print

  quality_surrogate_variables
}

get_qsva_dds <- function(dds) {

  design_formula = dds %>% design()
  sample_data<- dds %>% colData()

  quality_surrogate_variables <- get_quality_surrogate_variables(dds)


  if (is.vector(quality_surrogate_variables)) {
   quality_surrogate_variables %<>% data.frame(qSVA=.)
  }

  colData(dds) %<>% cbind(quality_surrogate_variables)

  design(dds) <- design_formula %>%
    terms() %>% attr("term.labels") %>%
    c(quality_surrogate_variables %>% colnames, .) %>% 
    reformulate
  
  dds
}


get_total_dds <- function(sample_data, species, filter_low_counts=FALSE, qSVA=FALSE) {
  # Collate count data
  total_count_data <- sample_data %>%
    row.names() %>%
    map(read_counts,species) %>%
    purrr::reduce(inner_join) %>%
    remove_gene_column()

  get_deseq2_dataset(
  total_count_data, sample_data,
  filter_low_counts=filter_low_counts, design_formula=~1,qSVA=qSVA)
}



#get res for given condition name
get_res <- function(comparison_name,sample_data,comparison_table, tpms, species,qSVA=FALSE,use_tx=FALSE,quant_method='salmon',tx_level=FALSE) {
  x=comparison_table %>% filter(comparison==comparison_name)
  sample_data %<>%
    tibble::rownames_to_column(var = "tmp_row_names") %>%
    mutate(!!x$condition_name:= factor(!!parse_expr(x$condition_name))) %>%
    filter(!!parse_expr(x$filter)) %>%
    tibble::column_to_rownames(var = "tmp_row_names")
  
  ##Ensure that conditions to be used in GSA comparisons are factors with
  # the correct base level set.
  sample_data[,x$condition_name] %<>% relevel(x$condition_base)

  if(use_tx){
    txi<-get_tximport(sample_data,quant_method,tx_level)
    dds <- DESeqDataSetFromTximport(txi, sample_data, x$formula %>% as.formula())
    dds <- dds[rowSums(counts(dds)) > 1, ]
    dds <- DESeq(dds, betaPrior = TRUE)
  }else{
    dds <- sample_data %>%
      row.names() %>%
      map(read_counts,species) %>%
      purrr::reduce(inner_join) %>%
      remove_gene_column() %>%
      get_deseq2_dataset(sample_data, design_formula = x$formula %>% as.formula(),qSVA=qSVA)
  }
  

  #res contains transcript, rather than gene
  if(use_tx & tx_level){
    res <- dds %>%
        get_deseq2_results(x$condition_name, x$condition, x$condition_base) %>%
        left_join(dds %>% get_raw_l2fc(sample_data, expr(!!sym(x$condition_name) == !!(x$condition)))) %>%
        dplyr::rename(transcript=gene)
  }else{
    res <- dds %>%
        get_deseq2_results(x$condition_name, x$condition, x$condition_base) %>%
        left_join(dds %>% get_raw_l2fc(sample_data, expr(!!sym(x$condition_name) == !!(x$condition))))
  }
  
  #fill summary table
  SUMMARY_TB <- get("SUMMARY_TB", envir = .GlobalEnv) %>% 
    add_row(Comparison = x$comparison, DESeq_model_formula = design(dds) %>% format(),
            Condition_tested = x$condition_name,
            Total_number_of_samples_data=sample_data %>% nrow(),
            Base_level_condition=x$condition_base,
            Number_of_samples_in_base_level_condition=sample_data %>% filter(!!parse_expr(x$condition_name)==x$condition_base)%>% nrow(),
            Sample_names_in_base_level_condition=sample_data %>% filter(!!parse_expr(x$condition_name)==x$condition_base)%>% pull(sample_name) %>% str_c(collapse = ','),
            Comparison_level_condition=x$condition,
            Number_of_samples_in_comparison_level_condition=sample_data %>% filter(!!parse_expr(x$condition_name)==x$condition)%>% nrow(),
            Sample_names_in_comparison_level_condition=sample_data %>% filter(!!parse_expr(x$condition_name)==x$condition)%>% pull(sample_name) %>% str_c(collapse = ','),
            p.adj.cutoff=0.05,
            Up_regulated=res %>% filter( padj < 0.05 & log2FoldChange > 0 ) %>% nrow(),
            Down_regulated=res %>% filter( padj < 0.05 & log2FoldChange < 0 ) %>% nrow(),
            D.E.total=res %>% filter( padj < 0.05) %>% nrow())
  
  assign("SUMMARY_TB", SUMMARY_TB,envir = .GlobalEnv)
  
  list(res, dds)
}


get_tximport<-function(sample_data,quant_method='salmon',tx_level=TRUE){
  quant_file <- get_transcript_quant_file(quant_method)
  quant_dirs <- sample_data %>%
    tibble::rownames_to_column(var="tmp") %>%
    pull("tmp")

  quant_files <- str_c("results/",quant_method,"_quant/", quant_dirs, "/", quant_file)
  names(quant_files) <- quant_dirs

  txi <- tximport(quant_files, type=quant_method, txOut = tx_level, tx2gene=get_transcripts_to_genes(), dropInfReps = TRUE)
  txi$Length <- read.csv(quant_files[1],sep = '\t',stringsAsFactors = FALSE) %>%
    dplyr::select(id=1,length=2) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(var='id')
  txi
}

get_total_dds_tximport <- function(sample_data,quant_method='salmon',tx_level=TRUE) {

  txi <- get_tximport(sample_data,quant_method,tx_level)

  total_dds <- DESeqDataSetFromTximport(txi, sample_data, ~1)
  total_dds <- DESeq(total_dds)

  total_dds
}

get_res_tx <- function(comparison_name,sample_data,comparison_table,quant_method='salmon',tx_level=FALSE) {
  x=comparison_table %>% filter(comparison==comparison_name)
  sample_data %<>%
  tibble::rownames_to_column(var = "tmp_row_names") %>%
    mutate(!!x$condition_name:= factor(!!parse_expr(x$condition_name))) %>%
    filter(!!parse_expr(x$filter)) %>%
    tibble::column_to_rownames(var = "tmp_row_names")

  ##Ensure that conditions to be used in GSA comparisons are factors with
  # the correct base level set.
  sample_data[,x$condition_name] %<>% relevel(x$condition_base)

  txi<-get_tximport(sample_data,quant_method,tx_level)

  dds <- DESeqDataSetFromTximport(txi, sample_data, x$formula %>% as.formula())
  dds <- dds[rowSums(counts(dds)) > 1, ]
  dds <- DESeq(dds)

  res <- dds %>%
    get_deseq2_results(x$condition_name, x$condition, x$condition_base) %>%
    left_join(dds %>% get_raw_l2fc_tx(sample_data, expr(!!sym(x$condition_name) == !!(x$condition))))

  #fill summary table
  SUMMARY_TB <- get("SUMMARY_TB", envir = .GlobalEnv) %>%
  add_row(Comparison = x$comparison, DESeq_model_formula = design(dds) %>% format(),
  Condition_tested = x$condition_name,
  Total_number_of_samples_data=sample_data %>% nrow(),
  Base_level_condition=x$condition_base,
  Number_of_samples_in_base_level_condition=sample_data %>% filter(!!parse_expr(x$condition_name)==x$condition_base)%>% nrow(),
  Sample_names_in_base_level_condition=sample_data %>% filter(!!parse_expr(x$condition_name)==x$condition_base)%>% pull(sample_name) %>% str_c(collapse = ','),
  Comparison_level_condition=x$condition,
  Number_of_samples_in_comparison_level_condition=sample_data %>% filter(!!parse_expr(x$condition_name)==x$condition)%>% nrow(),
  Sample_names_in_comparison_level_condition=sample_data %>% filter(!!parse_expr(x$condition_name)==x$condition)%>% pull(sample_name) %>% str_c(collapse = ','),
  p.adj.cutoff=0.05,
  Up_regulated=res %>% filter( padj < 0.05 & log2FoldChange > 0 ) %>% nrow(),
  Down_regulated=res %>% filter( padj < 0.05 & log2FoldChange < 0 ) %>% nrow(),
  D.E.total=res %>% filter( padj < 0.05) %>% nrow())

  assign("SUMMARY_TB", SUMMARY_TB,envir = .GlobalEnv)

  list(res, dds)
}



perform_rmats <- function(sample_data,comparison){

  top_dir<-str_c("results/rMATS/",species,"/",comparison)
  if (!dir.exists(top_dir)) {
     dir.create(top_dir,recursive=TRUE)
  }

  file.create(str_c(top_dir,"/b1.txt"))
  file.create(str_c(top_dir,"/b2.txt"))

  b1<-str_c(top_dir,"/b1.txt") %>% normalizePath()
  b2<-str_c(top_dir,"/b2.txt") %>% normalizePath()


  reps<-sample_data %>%
    group_by(!!parse_expr(x$condition_name)) %>%
    mutate(bam_file=str_c("results/final_bams/",sample_name,".",species,".bam",sep = '') %>% normalizePath) %>%
    summarise(replicates=str_c(bam_file,collapse = ','))

  reps %>%
    filter(!!parse_expr(x$condition_name) == x$condition_base) %>%
    pull(replicates) %>%
    write(b1)

  reps %>%
    filter(!!parse_expr(x$condition_name) == x$condition) %>%
    pull(replicates) %>%
    write(b2)

  cmd <- str_c("cd /opt/rMATS.4.0.2 && python rmats.py",
                "--b1", b1,
                "--b2", b2,
                "--gtf", str_c("data/",dir(path = "data/", pattern = str_c(species,"_ensembl_*"))) %>%
                list.files(pattern = '.gtf',full.names = T) %>% normalizePath,
                "--od",top_dir %>% normalizePath,
                "-t paired",
                "--nthread 12",
                "--tstat 12",
                "--readLength 150",
                "--cstat 0.05",
                "--libType fr-unstranded",sep = " "
  )

#2hours
  system(cmd)
}

#######
# comment out for now, needs project information to fill 
# this function in order to source the function.
# get_condition_res_tximport <- function(quant_method) {
#   sample_data <- data.frame(
#     condition=c(),
#     sample=c(),
#     filename=c(),
#     row.names=c("<SAMPLE1>", "<SAMPLE2>", etc)
#   )
# 
#   quant_file <- get_transcript_quant_file(quant_method)
#   quant_files <- str_c("results/", quant_method, "_quant/", sample_data$filename, "/", quant_file)
#   names(quant_files) <- sample_data$filename
# 
#   txi <- tximport(quant_files, type=quant_method, tx2gene=get_transcripts_to_genes(), reader=read_tsv)
# 
#   dds <- DESeqDataSetFromTximport(txi, sample_data, ~condition)
#   dds <- dds[rowSums(counts(dds)) > 1, ]
#   dds <- DESeq(dds)
# 
#   results <- dds %>% get_deseq2_results("condition", "<cond2>", "<cond1>")
# 
#   l2fc <- get_count_data(dds) %>%
#     mutate(l2fc=log2(() / ()) %>%
#              dplyr::select(gene, l2fc)
# 
#            results %<>% left_join(l2fc)
# 
#            return(results)
# }
