source("load_packages.R")
source("common_functions.R")

species="{{cookiecutter.species}}"

{% if cookiecutter.qSVA !="no" %}
qSVA<-TRUE
{% else %}
qSVA<-FALSE
{% endif %}

SAMPLE_NAMES <- c(condition1, condition2, etc) %>%
  outer(c(rep1, rep2, etc), str_c, sep="todo") %>%
  t %>%
  as.vector

# n.b. Ensure that conditions to be used in GSA comparisons are factors with
# the correct base level set.
SAMPLE_DATA <- data.frame(
  condition=...,
  species=species,
  row.names=SAMPLE_NAMES%>%str_c(species,sep = '.')
)

SAMPLE_DATA %<>% tibble::rownames_to_column(var = "tmp_row_names") %>%
                 filter(species==!!species) %>%
                 tibble::column_to_rownames(var = "tmp_row_names")

#example can be found https://github.com/sidbdri/cookiecutter-sargasso-de_analysis_skeleton
COMPARISON_TABLE<-tribble(
~comparison, ~formula, ~condition_name, ~condition, ~condition_base, ~filter,
#"P10_Ctx_KO_vs_WT", "~genotype", "genotype", "KO", "WT", "age=='P10' & region=='Ctx'",
)

SUMMARY_TB<-setNames(data.frame(matrix(ncol = 14, nrow = 0)),
                     c("Comparison", "DESeq_model_formula", "Condition_tested", "Total_number_of_samples_data",
                       "Base_level_condition", "Number_of_samples_in_base_level_condition",
                       "Sample_names_in_base_level_condition",
                       "Comparison_level_condition", "Number_of_samples_in_comparison_level_condition",
                       "Sample_names_in_comparison_level_condition",
                       "p.adj.cutoff", "Up_regulated", "Down_regulated", "D.E.total"))


#####

total_dds_data <- get_total_dds(SAMPLE_DATA,qSVA=qSVA)
total_vst <- total_dds_data %>% varianceStabilizingTransformation
total_vst %>% plot_pca_with_labels(intgroup=c("condition"))
total_vst %>% plot_heat_map(SAMPLE_DATA %>% 
                              mutate(sample_info=str_c(condition, ..., sep=":")) %>% 
                              extract2("sample_info"))

plot_count_distribution(total_dds_data, norm=F)
plot_count_distribution(total_dds_data, norm=T)

#####

gene_info <- get_gene_info(species)
gene_lengths <- read_csv(str_c("data/",species,"_ensembl_",{{cookiecutter.ensembl_version}},"/gene_lengths.csv", sep=""))

results <- total_dds_data %>% get_count_data() 

fpkms <- results %>% 
  get_fpkms(gene_lengths, colnames(results) %>% tail(-1), "_fpkm")


# fpkms %<>% mutate(
#   P10_Ctx_KO_fpkm_avg = !!get_avg_fpkm(filter="age=='P10' & genotype=='KO' & region=='Ctx'"),
#   P10_Piri_KO_fpkm_avg = !!get_avg_fpkm(filter="age=='P10' & genotype=='KO' & region=='Piri'")
# ) 
fpkms %<>% mutate(
  <CONDITION1>_fpkm_avg = !!get_avg_fpkm(SAMPLE_DATAfilter="condition1=='' & condition2==''"),
  etc.
)

results %<>% 
  inner_join(fpkms) %>%
  inner_join(gene_info) %>% 
  inner_join(gene_lengths)

##run all get_res functions and add to results object
COMPARISON_TABLE %>% pull(comparison) %>% walk ( function(x){
  res_name<-str_c(x,'res',sep = '_')
  assign(str_c(x,'res',sep = '_'), get_res(x,SAMPLE_DATA,COMPARISON_TABLE,qSVA=qSVA),envir = .GlobalEnv)
  
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
  write_csv(str_c("results/differential_expression/deseq2_results_",species,".csv"))

results %>% 
  dplyr::select(gene, gene_name, chromosome, description, entrez_id, gene_type,
                gene_length, max_transcript_length,
         dplyr::contains("_fpkm"), 
         COMPARISON_TABLE %>% pull(comparison) %>%
           sapply(FUN = function(x) results %>% colnames() %>% str_which(str_c("^",x,sep =''))) %>%
           as.vector() %>% unique(), 
         -dplyr::ends_with(".stat")) %>% 
  write_csv(str_c("results/differential_expression/deseq2_results_fpkm_",species,".csv"))

SUMMARY_TB %>%
  write_csv(str_c("results/differential_expression/de_summary_",species,".csv"))

##### GO analyses

expressed_genes <- get_total_dds(SAMPLE_SAMPLE_DATA,filter_low_counts=TRUE) %>% get_count_data()

COMPARISON_TABLE %>% pull(comparison) %>% walk( function(x){
  p_str=str_c(x,'padj',sep = '.')
  l2fc_str=str_c(x,'l2fc',sep = '.')
  
  get("results",envir = .GlobalEnv) %>% 
    filter( get(p_str) < 0.05 ) %>% 
    perform_go_analyses(expressed_genes, x, species)
  
  get("results",envir = .GlobalEnv) %>%
    filter( get(p_str) < 0.05  & get(l2fc_str) > 0 ) %>% 
    perform_go_analyses(expressed_genes, str_c(x,'up',sep = '.'),species)
  
  get("results",envir = .GlobalEnv) %>%
    filter( get(p_str) < 0.05  & get(l2fc_str) < 0 ) %>% 
    perform_go_analyses(expressed_genes, str_c(x,'down',sep = '.'),species)
})

##### Gene set enrichment analysis

gene_set_categories <- list("CURATED", "MOTIF", "GO")

gene_sets <- gene_set_categories %>% 
  map(function(x) get_gene_sets(species, x))


COMPARISON_TABLE %>% pull(comparison) %>% walk( function(x){
  x <- COMPARISON_TABLE %>% filter(x==comparison)
  res<-str_c(x$comparison,'res',sep = '_') %>% get(envir = .GlobalEnv)
  
  camera_results <- get('gene_sets',envir = .GlobalEnv) %>% 
    map(function(y) get_camera_results( res[[2]],y,gene_info, x$formula %>% as.formula()))
  
  assign(str_c(x$comparison,'camera_results',sep = '_'),camera_results,envir = .GlobalEnv)
  
  for (category in seq(1:length(gene_set_categories))) {
    de_res <- results %>% dplyr::select(
      gene, gene_name, entrez_id, 
      starts_with(str_c(x$comparison, ".")), 
      -starts_with(str_c(x$comparison, ".stat")))  
    write_camera_results(
      gene_set_categories[[category]], gene_sets[[category]], x$comparison, species,
      de_res, camera_results[[category]])
  }
})

# results %>% plot_gene_set(gene_sets[[3]], "GO_<go_term>", "condition.stat")
# results %>% get_gene_set_results(gene_sets[[3]], "GO_<go_term>", "condition.pval") %>% head

##### Salmon/tximport analysis

# results_salmon <- get_total_dds_tximport("salmon") %>% 
#   get_count_data()
# 
# salmon_tpms <- SAMPLE_NAMES %>%
#   map(get_salmon_tpms) %>%
#   purrr::reduce(inner_join)
# 
# results_salmon %<>% 
#   inner_join(fpkms) %>%
#   inner_join(gene_info) %>% 
#   inner_join(gene_lengths)
# 
# results_salmon %<>% 
#   left_join(get_condition_res_tximport("salmon", "quant.sf"), by="gene") %>%
#   dplyr::rename(condition.l2fc=log2FoldChange,
#                 condition.raw_l2fc=l2fc,
#                 condition.pval=pvalue,
#                 condition.padj=padj)
# 
# results_salmon %>% 
#   dplyr::select(gene, gene_name, chromosome, description, entrez_id, gene_type,
#                 gene_length, max_transcript_length,
#                 everything(), -dplyr::contains("_tpm")) %>%
#   write_csv("results/differential_expression/deseq2_salmon_results.csv")
#            
# results_salmon %>% 
#   dplyr::select(gene, gene_name, chromosome, description, entrez_id, gene_type,
#                 gene_length, max_transcript_length,
#          dplyr::contains("_tpm"), 
#          starts_with(condition), etc.) %>% 
#   write_csv("results/differential_expression/deseq2_salmon_results_fpkm.csv")
# 
# ##### Kallisto/tximport analysis
# 
# results_kallisto <- get_total_dds_tximport("kallisto") %>% 
#   get_count_data() 
# 
# kallisto_tpms <- SAMPLE_NAMES %>%
#   map(get_kallisto_tpms) %>%
#   purrr::reduce(inner_join)
# 
# results_kallisto %<>% 
#   inner_join(kallisto_tpms) %>%
#   inner_join(gene_info) %>% 
#   inner_join(gene_lengths)
# 
# results_kallisto %<>% 
#   left_join(get_condition_res_tximport("kallisto", "abundance.tsv"), by="gene") %>%
#   dplyr::rename(condition.l2fc=log2FoldChange,
#                 condition.raw_l2fc=l2fc,
#                 condition.pval=pvalue,
#                 condition.padj=padj)
# 
# results_kallisto %>% 
#   dplyr::select(gene, gene_name, chromosome, description, entrez_id, gene_type,
#                 gene_length, max_transcript_length,
#                 everything(), -dplyr::contains("_tpm")) %>%
#   write_csv("results/differential_expression/deseq2_kallisto_results.csv")
#            
# results_kallisto %>% 
#   dplyr::select(gene, gene_name, chromosome, description, entrez_id, gene_type,
#                 gene_length, max_transcript_length,
#          dplyr::contains("_tpm"), 
#          starts_with(condition), etc.) %>% 
#   write_csv("results/differential_expression/deseq2_kallisto_results_fpkm.csv")
