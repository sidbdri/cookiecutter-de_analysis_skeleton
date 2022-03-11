get_run_reactome_job_strings <- function(comparison_table) {
  comparison_table %>%
    pull(comparison) %>% 
    expand.grid(c('all','up','down')) %>% 
    tidyr::unite(col = 'job_string', sep = ';') %>% 
    pull(job_string)
}

run_reactome <- function(job_string,result_tbl,p.adj.cutoff,expressed_genes,species,out_dir,...){
  
  futile.logger::flog.info(paste0('running Reactome for:',job_string))
  
  c(comparison_name, direction) %<-% (strsplit(job_string, split = ';') %>% unlist())
  c(p_str, l2fc_str) %<-% str_c(comparison_name, c('.padj', '.l2fc'))
  result_tbl %<>% dplyr::select(gene, entrez_id, contains(comparison_name))
  
  if (direction == 'up') {
    result_tbl %<>% 
      filter_at(vars(p_str), any_vars(. < p.adj.cutoff)) %>% filter_at( vars(l2fc_str), any_vars(. > 0))
  } else if ((direction == 'down')) {
    result_tbl %<>% 
      filter_at(vars(p_str), any_vars(. < p.adj.cutoff)) %>% filter_at( vars(l2fc_str), any_vars(. < 0))
  } else {
    result_tbl %<>% filter_at( vars(p_str), any_vars(. < p.adj.cutoff))
  }
  
  if (result_tbl %>% nrow == 0) {
    message("No significant genes supplied.")
    return()
  }
  
  top_dir <-  file.path(out_dir,species,comparison_name)
  if (!dir.exists(top_dir)) {
    dir.create(top_dir, recursive = TRUE)
  }
  ret <- perform_pathway_enrichment(result_tbl, expressed_genes, comparison_name, direction, species, out_dir)
  
  ret
}

perform_pathway_enrichment <- function(significant_genes, expressed_genes, 
                                       comparison_name, file_prefix, species, out_dir="results/differential_expression/reactome/") {
  
  gene_info <- get_gene_info(species)
  
  universe <- expressed_genes %>% 
    inner_join(gene_info) %>%
    filter(!is.na(entrez_id)) %>%
    pull("entrez_id")
  
  gene_list <- significant_genes %>% 
    filter(!is.na(entrez_id)) %>%
    pull("entrez_id")
  
  pathways <- enrichPathway(
    gene = gene_list, organism = species, universe = as.character(universe), 
    pvalueCutoff = 0.1, readable = T) %>%
    as.data.frame()
  
  top_dir <-  file.path(out_dir,species,comparison_name)
  if (!dir.exists(top_dir)) {
    dir.create(top_dir, recursive = TRUE)
  }
  
  ret <- data.frame()
  
  if (pathways %>% nrow() > 0) {
    pathways %>%
      dplyr::select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, geneID) %>%
      write_csv(file.path(top_dir, str_c(comparison_name, file_prefix, "_reactome.csv")), na = "")
    
    ret <- pathways %>% dplyr::select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, geneID)
  }
  
  ret
}