source("meta_data.R")


SPECIES="{{cookiecutter.species}}"

TX_LEVEL=FALSE
QUANT_METHOD='salmon'
USE_TX=TRUE

{% if cookiecutter.qSVA !="no" %}
qSVA<-TRUE
{% else %}
qSVA<-FALSE
{% endif %}


output_folder = 'results/differential_expression/de_tx'
if (!dir.exists(output_folder)) dir.create(output_folder,recursive=TRUE)

#####

total_dds_data <- get_total_dds_tximport(SAMPLE_DATA,QUANT_METHOD,TX_LEVEL)

total_vst <- total_dds_data %>% varianceStabilizingTransformation
total_vst %>% plot_pca_with_labels(intgroup=c("condition"))
total_vst %>% plot_heat_map(SAMPLE_DATA %>% 
                              mutate(sample_info=str_c(condition, ..., sep=":")) %>% 
                              extract2("sample_info"))

plot_count_distribution(total_dds_data, norm=F)
plot_count_distribution(total_dds_data, norm=T)

#####

gene_info <- get_gene_info(SPECIES)
gene_lengths <- read_csv(str_c("data/",SPECIES,"_ensembl_",{{cookiecutter.ensembl_version}},"/gene_lengths.csv", sep=""))

results <- total_dds_data %>% get_count_data()

txi<-get_tximport(SAMPLE_DATA,QUANT_METHOD,TX_LEVEL)
tpms<-txi$abundance %>% as.data.frame() %>%
  rename_all(funs(str_c(.,"_tpm")))


# fpkms <- results %>%
#   get_fpkms(gene_lengths, colnames(results) %>% tail(-1), "_fpkm")


# fpkms %<>% mutate(
#   P10_Ctx_KO_fpkm_avg = !!get_avg_fpkm(filter="age=='P10' & genotype=='KO' & region=='Ctx'"),
#   P10_Piri_KO_fpkm_avg = !!get_avg_fpkm(filter="age=='P10' & genotype=='KO' & region=='Piri'")
# ) 
# fpkms %<>% mutate(
#   <CONDITION1>_fpkm_avg = !!get_avg_fpkm(SAMPLE_DATAfilter="condition1=='' & condition2==''"),
#   etc.
# )

if(TX_LEVEL){
    tpms %<>% tibble::rownames_to_column('transcript')
    tpms %<>%inner_join(get_transcripts_to_genes(SPECIES))  %>%
        left_join(get_avg_tpm(tpms,TX_LEVEL))

    results %<>% rename(transcript=gene) %>%
        left_join(txi$Length %>% rename(transcript_length=length) %>% tibble::rownames_to_column('transcript')) %>%
        left_join(get_transcripts_to_genes(SPECIES)) %>%
        left_join(tpms) %>%
        left_join(gene_info) %>%
        left_join(gene_lengths)

    #number of transcript
    results %<>% left_join(
    results %>%
        dplyr::select(gene,transcript) %>%
        group_by(gene) %>% summarise(number_of_transcript = n()))
}else{
    tpms %<>% tibble::rownames_to_column('gene')
    tpms %<>%inner_join(get_transcripts_to_genes(SPECIES))  %>%
        left_join(get_avg_tpm(tpms,TX_LEVEL))

    results %<>%
        left_join(tpms) %>%
        left_join(gene_info) %>%
        left_join(gene_lengths)

}

# #workout avg tpm
# results %<>% dplyr::select(gene,dplyr::contains("_tpm")) %>% group_by(gene) %>%
#                 summarise_all(.funs = sum) %>%
#                 mutate(sum=rowSums(dplyr::select(., dplyr::contains("_tpm"))),
#                 n=ncol(dplyr::select(., dplyr::contains("_tpm")))) %>%
#                 mutate(avg_gene_tpm=sum/n) %>%
#                 dplyr::select(gene,avg_gene_tpm) %>% right_join(results)

##run all get_res functions and add to results object
COMPARISON_TABLE %>% pull(comparison) %>% walk ( function(x){
  res_name<-str_c(x,'res',sep = '_')
  assign(str_c(x,'res',sep = '_'), get_res(x,SAMPLE_DATA,COMPARISON_TABLE,tpms,use_tx=USE_TX,quant_method=QUANT_METHOD,tx_level=TX_LEVEL),envir = .GlobalEnv)

  res <-get(res_name, envir = .GlobalEnv)
  results<-get("results",envir = .GlobalEnv) %>%
    left_join(res[[1]]) %>%
    dplyr::rename(!!str_c(x,'l2fc',sep = '.'):=log2FoldChange,
                  !!str_c(x,'raw_l2fc',sep = '.'):=raw_l2fc,
                  !!str_c(x,'stat',sep = '.'):=stat,
                  !!str_c(x,'pval',sep = '.'):=pvalue,
                  !!str_c(x,'padj',sep = '.'):=padj)
  plot_pvalue_distribution(results, str_c(x,'pval',sep = '.'))
  assign("results", results,envir = .GlobalEnv)
})

#save results
column_inclued<-c('gene')
if(TX_LEVEL){
    column_inclued<-c('transcript','transcript_length','gene','number_of_transcript')
}

#count
results %>%
    dplyr::select(!!column_inclued, gene_name, chromosome, description, entrez_id, gene_type,
    everything(),  -gene_length, -max_transcript_length, -dplyr::contains("_tpm"), -dplyr::ends_with(".stat")) %>%
         write_csv(str_c(output_folder,"/deseq2_results_count_",SPECIES,"_tx_",
                    ifelse(TX_LEVEL,"transcript","gene"),"_",QUANT_METHOD,".csv"))

if(TX_LEVEL){
  results %>% 
    dplyr::select(!!column_inclued, gene_name, chromosome, description, entrez_id, gene_type,
                  -gene_length, -max_transcript_length, dplyr::contains("_tpm"), 
                  COMPARISON_TABLE %>% pull(comparison) %>%
                    sapply(FUN = function(x) results %>% colnames() %>% str_which(str_c("^",x,sep =''))) %>%
                    as.vector() %>% unique(), 
                  -dplyr::ends_with(".stat")) %>% 
    write_csv(str_c(output_folder,"/deseq2_results_tpm_",SPECIES,"_tx_",
                    ifelse(TX_LEVEL,"transcript","gene"),"_",QUANT_METHOD,".csv"))
}else{
  results %>% 
    dplyr::select(!!column_inclued, gene_name, chromosome, description, entrez_id, gene_type,
                  dplyr::contains("_tpm"), 
                  COMPARISON_TABLE %>% pull(comparison) %>%
                    sapply(FUN = function(x) results %>% colnames() %>% str_which(str_c("^",x,sep =''))) %>%
                    as.vector() %>% unique(), 
                  -dplyr::ends_with(".stat")) %>% 
    write_csv(str_c(output_folder,"/deseq2_results_tpm_",SPECIES,"_tx_",
                    ifelse(TX_LEVEL,"transcript","gene"),"_",QUANT_METHOD,".csv"))
}
SUMMARY_TB %>%
    write_csv(str_c(output_folder,"/de_summary_",SPECIES,"_tx_",
              ifelse(TX_LEVEL,"transcript","gene"),"_",QUANT_METHOD,".csv"))


# library (knitr)
# sink(str_c("results/differential_expression/de_summary_",
#             SPECIES,"_tx_",ifelse(TX_LEVEL,"transcript","gene"),"_",QUANT_METHOD,".md"))
# kable(SUMMARY_TB, format = 'markdown')
# sink()

##### GO analyses
if( !USE_TX | !TX_LEVEL ){
    expressed_genes <- total_dds_data %>% get_count_data()
    COMPARISON_TABLE %>% pull(comparison) %>% walk( function(x){
      p_str=str_c(x,'padj',sep = '.')
      l2fc_str=str_c(x,'l2fc',sep = '.')

      get("results",envir = .GlobalEnv) %>%
        filter( get(p_str) < 0.05 ) %>%
        perform_go_analyses(expressed_genes, x, SPECIES)

      get("results",envir = .GlobalEnv) %>%
        filter( get(p_str) < 0.05  & get(l2fc_str) > 0 ) %>%
        perform_go_analyses(expressed_genes, str_c(x,'up',sep = '.'),SPECIES)

      get("results",envir = .GlobalEnv) %>%
        filter( get(p_str) < 0.05  & get(l2fc_str) < 0 ) %>%
        perform_go_analyses(expressed_genes, str_c(x,'down',sep = '.'),SPECIES)
    })

    ##### Gene set enrichment analysis

    gene_set_categories <- list("CURATED", "MOTIF", "GO")

    gene_sets <- gene_set_categories %>%
      map(function(x) get_gene_sets(SPECIES, x))


    COMPARISON_TABLE %>% pull(comparison) %>% walk( function(x){
      x <- COMPARISON_TABLE %>% filter(x==comparison)
      res<-str_c(x$comparison,'res',sep = '_') %>% get(envir = .GlobalEnv)

      camera_results <- get('gene_sets',envir = .GlobalEnv) %>%
        map(function(y) get_camera_results( res[[2]],y,gene_info))

      assign(str_c(x$comparison,'camera_results',sep = '_'),camera_results,envir = .GlobalEnv)

      for (category in seq(1:length(gene_set_categories))) {
        de_res <- results %>% dplyr::select(
          gene, gene_name, entrez_id,
          starts_with(str_c(x$comparison, ".")),
          -starts_with(str_c(x$comparison, ".stat")))
        write_camera_results(
          gene_set_categories[[category]], gene_sets[[category]], x$comparison, SPECIES,
          de_res, camera_results[[category]])
      }
    })
}

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
