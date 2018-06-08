source("load_packages.R")
source("common_functions.R")


SAMPLE_NAMES <- c(condition1, condition2, etc) %>%
  outer(c(rep1, rep2, etc), str_c, sep="todo") %>%
  t %>%
  as.vector

# n.b. Ensure that conditions to be used in GSA comparisons are factors with
# the correct base level set.
SAMPLE_DATA <- data.frame(
  condition=...,
  row.names=SAMPLE_NAMES
)

#example can be found https://github.com/sidbdri/cookiecutter-sargasso-de_analysis_skeleton
comparison_table<-tribble(
~comparision, ~fomular, ~condition_name, ~condition, ~condition_base, ~filter,...,
#"P10_Ctx_KO_vs_WT", "~genotype", "genotype", "KO", "WT", "age=='P10' & region=='Ctx'",...,
)

SUMMARY_TB<-setNames(data.frame(matrix(ncol = 14, nrow = 0)),
                     c("Comparision", "DESeq_model_formula", "Condition_tested",
                       "Base_level_condition","Total_number_of_samples_data","Number_of_samples_in_base_level_condition",
                       "Sample_names_in_base_level_condition",
                       "Comparison_level_condition","Number_of_samples_in_comparison_level_condition",
                       "Sample_names_in_comparison_level_condition",
                       "p.adj.cutoff","Up_regulated","Down_regulated","D.E.total"))


get_total_dds <- function(filter_low_counts=FALSE) {
  # Collate count data
  total_count_data <- SAMPLE_DATA %>% 
    row.names() %>% 
    map(read_counts) %>%
    purrr::reduce(inner_join) %>%
    remove_gene_column()
    
  get_deseq2_dataset(
    total_count_data, SAMPLE_DATA,
    filter_low_counts=filter_low_counts, design_formula=~1)
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

#get res for given condition name
get_res <- function(comparision_name) {
  x=comparison_table %>% filter(comparision==comparision_name)
  sample_data <- SAMPLE_DATA %>%
    tibble::rownames_to_column(var = "tmp_row_names") %>%
    mutate(!!x$condition_name:= factor(!!parse_expr(x$condition_name))) %>%
    filter(!!parse_expr(x$filter)) %>%
    tibble::column_to_rownames(var = "tmp_row_names")

  ##Ensure that conditions to be used in GSA comparisons are factors with
  # the correct base level set.
  sample_data[,x$condition_name] %<>% relevel(x$condition_base)


  dds <- sample_data %>%
    row.names() %>%
    map(read_counts) %>%
    purrr::reduce(inner_join) %>%
    remove_gene_column() %>%
    get_deseq2_dataset(sample_data, design_formula = x$fomular %>% as.formula() )

  res <- dds %>%
    get_deseq2_results(x$condition_name, x$condition, x$condition_base) %>%
    left_join(dds %>% get_raw_l2fc(sample_data, expr(!!sym(x$condition_name) == !!(x$condition))))

  #fill sumary table
  SUMMARY_TB <- get("SUMMARY_TB", envir = .GlobalEnv) %>% 
    add_row(Comparision = x$comparision, DESeq_model_formula = x$fomular, 
            Condition_tested = x$condition_name,
            Base_level_condition=x$condition_base,
            Total_number_of_samples_data=sample_data %>% nrow(),
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

  list(res, dds %<>% varianceStabilizingTransformation)
}

# get_condition_res <- function() {
#   sample_data <- SAMPLE_DATA %>% 
#     filter_with_rownames(...)
# 
#   dds <- sample_data %>% 
#     row.names() %>% 
#     map(read_counts) %>% 
#     purrr::reduce(inner_join) %>%
#     remove_gene_column() %>% 
#     get_deseq2_dataset(sample_data, design_formula=~condition)
# 
#   #cond2 is to be compared against the baseline, cond1.
#   res <- dds %>% 
#     get_deseq2_results("<condition>", "<cond2>", "<cond1>") %>%
#     left_join(dds %>% get_raw_l2fc(sample_data, condition="<cond2>"))
# 
# 
#   list(res, dds %<>% varianceStabilizingTransformation)
# }
  
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


# fpkms %<>% mutate(
#   P10_Ctx_KO_fpkm_avg = !!get_avg_fpkm(filter="age=='P10' & genotype=='KO' & region=='Ctx'"),
#   P10_Piri_KO_fpkm_avg = !!get_avg_fpkm(filter="age=='P10' & genotype=='KO' & region=='Piri'")
# ) 
fpkms %<>% mutate(
  <CONDITION1>_fpkm_avg = !!get_avg_fpkm(filter="condition1=='' & condition2==''"),
  etc.
)



results %<>% 
  inner_join(fpkms) %>%
  inner_join(gene_info) %>% 
  inner_join(gene_lengths)


##run all get_res functions and add to results object
comparison_table %>% pull(comparision) %>% walk ( function(x){
  res_name<-str_c(x,'res',sep = '_')
  assign(str_c(x,'res',sep = '_'), get_res(x),envir = .GlobalEnv)
  
  res <-get(res_name, envir = .GlobalEnv)
  results<-get("results",envir = .GlobalEnv) %>% 
    left_join(res[[1]], by="gene") %>%
    dplyr::rename(!!str_c(x,'l2fc',sep = '.'):=log2FoldChange,
                  !!str_c(x,'raw_l2fc',sep = '.'):=raw_l2fc,
                  !!str_c(x,'stat',sep = '.'):=stat,
                  !!str_c(x,'pval',sep = '.'):=pvalue,
                  !!str_c(x,'padj',sep = '.'):=padj)
  plot_pvalue_distribution(results, str_c(x,'pval',sep = '.'))
  assign("results", results,envir = .GlobalEnv)
}) 




#save results
results %>% 
  dplyr::select(gene, gene_name, chromosome, description, entrez_id, gene_type,
                gene_length, max_transcript_length,
                everything(), -dplyr::contains("_fpkm"), -dplyr::ends_with(".stat")) %>%
  write_csv("results/differential_expression/deseq2_results.csv")

results %>% 
  dplyr::select(gene, gene_name, chromosome, description, entrez_id, gene_type,
                gene_length, max_transcript_length,
         dplyr::contains("_fpkm"), 
         comparison_table %>% pull(comparision) %>%
           sapply(FUN = function(x) results %>% colnames() %>% str_which(str_c("^",x,sep =''))) %>%
           as.vector() %>% unique(), 
         -dplyr::ends_with(".stat")) %>% 
  write_csv("results/differential_expression/deseq2_results_fpkm.csv")

SUMMARY_TB %>%
write_csv("results/differential_expression/de_summary.csv")

##### GO analyses

expressed_genes <- get_total_dds(TRUE) %>% get_count_data()

comparison_table %>% pull(comparision) %>% walk( function(x){
  p_str=str_c(x,'padj',sep = '.')
  l2fc_str=str_c(x,'l2fc',sep = '.')
  
  get("results",envir = .GlobalEnv) %>% 
    filter( get(p_str) < 0.05 ) %>% 
    perform_go_analyses(expressed_genes, x)
  
  get("results",envir = .GlobalEnv) %>%
    filter( get(p_str) < 0.05  & get(l2fc_str) > 0 ) %>% 
    perform_go_analyses(expressed_genes, str_c(x,'up',sep = '.'))
  
  get("results",envir = .GlobalEnv) %>%
    filter( get(p_str) < 0.05  & get(l2fc_str) < 0 ) %>% 
    perform_go_analyses(expressed_genes, str_c(x,'down',sep = '.'))
})



##### Gene set enrichment analysis

gene_set_categories <- list("CURATED", "MOTIF", "GO")

gene_sets <- gene_set_categories %>% 
  map(function(x) get_gene_sets("{{cookiecutter.species}}",x))


comparison_table %>% pull(comparision) %>% walk( function(x){
  x=comparison_table %>% filter(x==comparision)
  res<-str_c(x$comparision,'res',sep = '_') %>% get(envir = .GlobalEnv)
  
  camera_results <- get('gene_sets',envir = .GlobalEnv) %>% 
    map(function(y) get_camera_results( res[[2]],y,gene_info, x$fomular %>% as.formula()))
  
  assign(str_c(x$comparision,'camera_results',sep = '_'),camera_results,envir = .GlobalEnv)
  
  
  for (category in seq(1:length(gene_set_categories))) {
    write_camera_results(
      gene_set_categories[[category]], gene_sets[[category]], x$comparision,
      results,camera_results[[category]])
  }
  
})

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
