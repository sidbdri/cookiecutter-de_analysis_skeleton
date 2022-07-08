GO_TERMS <- Term(GOTERM) %>% as.data.frame() %>% tibble::rownames_to_column(var="GO.ID")
colnames(GO_TERMS) <- c("GO.ID", "FullTerm")

get_run_topgo_job_strings <- function(comparison_table) {
  comparison_table %>%
    pull(comparison) %>% 
    expand.grid(c('all','up','down'), c("BP", "MF", "CC")) %>% 
    tidyr::unite(col = 'job_string', sep = ';') %>% 
    pull(job_string) 
}

run_topgo <- function(job_string, result_tbl, p.adj.cutoff, expressed_genes, species, out_dir, ...) {

  futile.logger::flog.info(paste0('running topGO for:',job_string))
  
  c(comparison_name, direction, ontology) %<-% (strsplit(job_string, split = ';') %>% unlist())
  c(p_str, l2fc_str) %<-% str_c(comparison_name, c('.padj', '.l2fc'))
  result_tbl %<>% dplyr::select(gene, contains(comparison_name))
  
  if (direction == 'up') {
    result_tbl %<>% 
      filter_at(vars(p_str), any_vars(. < p.adj.cutoff)) %>% filter_at(vars(l2fc_str), any_vars(. > 0))
  } else if ((direction == 'down')) {
    result_tbl %<>% 
      filter_at(vars(p_str), any_vars(. < p.adj.cutoff)) %>% filter_at(vars(l2fc_str), any_vars(. < 0))
  } else {
    result_tbl %<>% filter_at(vars(p_str), any_vars(. < p.adj.cutoff))
  }
  
  if (result_tbl %>% nrow == 0) {
    message("No significant genes supplied.")
    return()
  }
  
  top_dir <-  file.path(out_dir,species,comparison_name)
  if (!dir.exists(top_dir)) {
    dir.create(top_dir, recursive = TRUE)
  }
  
  ret <- perform_go_analysis(expressed_genes, result_tbl, ontology, species, top_dir ,comparison_name)
  message( "Job finished:  ", job_string)
  
  ret %>% 
    extract2('go_results') %>% 
    inner_join(GO_TERMS) %>% 
    dplyr::mutate(Term = FullTerm) %>% 
    dplyr::select(-FullTerm) %>%
    dplyr::rename(annotated_in_background = Annotated, annotated_in_gene_set = Significant,
                  expected_annotated_in_gene_set = Expected, p.value = weight_fisher) %>% # Changing column names in results
    write_csv(file.path(top_dir, str_c(comparison_name, '.' , direction, "_go_", ontology %>% tolower, ".csv")), na = "")
  
  ret
}

perform_go_analysis <- function(gene_universe, significant_genes, ontology="BP", 
                                species, top_dir, comparison_name) {
  
  gene_list <- (gene_universe$gene %in% significant_genes$gene) %>% as.integer %>% factor
  names(gene_list) <- gene_universe$gene
  
  mapping <- switch(species,
                    mouse = "org.Mm.eg.db",
                    rat = "org.Rn.eg.db",
                    human = "org.Hs.eg.db")
  
  go_data <- new("topGOdata", ontology = ontology, allGenes = gene_list,
                 annot = annFUN.org, mapping = mapping, ID = "Ensembl")
  
  result_weight <- go_data %>% runTest('weight01', 'fisher')
  # result_classic <- go_data %>% runTest('classic', 'fisher')
  # result_elim <- go_data %>% runTest('elim', 'fisher')
  result_weight %>% print()

  # this could refactor somehow so the GO/REACTOME use the same function
  # but currently the parallel code only source one file, we need to change that
  # so GO/REACTOME source the same file containing this function.
  .calculate_odds_ratio_go <- function(go_data,go_result,term){
    # @numM: number of gene annotated to the term
    # @total_anotated: total number of gene (gene universe)
    # @numHits: number of gene annotated to the term AND in the significant list
    # @numSig: total number of significant gene
    c(numM,numHits,expected) %<-% termStat(go_data, term)%>%unlist
    c(total_anotated,numSig) %<-% go_result@geneData[c("Annotated","Significant")]

    # method 1
    # (numHits/(numM - numHits)) / ((numSig - numHits)/(total_anotated - numM - numSig + numHits))
    # method 2
    contMat <- cbind(sig = c(numHits, numSig - numHits),
                     notSig = c(numM - numHits, total_anotated - numM - numSig + numHits)) %>%
      set_rownames(c("anno", "notAnno"))
    odd.ratio <- fisher.test(contMat) %>% extract2('estimate') %>% unname()
    odd.ratio
  }
  
  go_results <- go_data %>% GenTable(weight_fisher = result_weight, orderBy = "weight_fisher", topNodes = 150)
  if (go_results %>% nrow() > 0) {
    go_results$odds_ratio <- sapply(go_results[,c('GO.ID')], function(x) .calculate_odds_ratio_go(go_data,result_weight, term=x))
  }

  gene_info <- get_gene_info(species)
  go_results$Genes <- sapply(go_results[,c('GO.ID')], function(x) get_significant_genes(x, go_data, gene_info))
  
  list(go_results = go_results,
       # go_data = go_data,
       # result_classic = result_classic,
       # result_elim = result_elim,
       result_weight = result_weight
  )
}

get_significant_genes <- function(term, GOdata, gene_info) {
  genes_for_term <- GOdata %>% genesInTerm(term) %>% extract2(1)
  significant_genes <- GOdata %>% sigGenes
  significant_genes_for_term <- genes_for_term %>% intersect(significant_genes)
  
  gene_info %>% 
    filter(gene %in% significant_genes_for_term) %>% 
    dplyr::select(gene_name) %>% 
    extract2(1) %>% 
    paste(collapse = ", ")
}


