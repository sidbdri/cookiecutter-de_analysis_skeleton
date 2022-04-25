get_gene_sets <- function(species, gene_set_name) {
  str_c("data/",species,"_ensembl_{{cookiecutter.ensembl_version}}/msigdb/v{{cookiecutter.msigdb_version}}/", gene_set_name, ".all.v{{cookiecutter.msigdb_version}}.entrez.gmt.Rdata") %>% readRDS
}

get_run_gene_set_analysis_job_strings <- function(comparison_table, gene_set_categories) {
  comparison_table %>%
    pull(comparison) %>% 
    expand.grid(gene_set_categories) %>% 
    tidyr::unite(col = 'job_string', sep = ';') %>% 
    pull(job_string)
}

run_gene_set_analysis <- function(job_string, result_tbl, dds_list, list_of_gene_sets, gene_info,species, out_dir, ...) {
  
  futile.logger::flog.info(paste0('running GS for:',job_string))
  
  c(comparison_name, category) %<-% (strsplit(job_string, split = ';') %>% unlist())
  dds <- dds_list[[str_c(comparison_name,'dds',sep = '_')]]
  
  camera_results <- list_of_gene_sets[category] %>%
    map(function(category_gene_sets) {
      get_camera_results(dds, category_gene_sets, gene_info)
    })
  
  de_res <- result_tbl %>% dplyr::select(
    gene, gene_name, entrez_id,
    starts_with(str_c(comparison_name, ".")),
    -starts_with(str_c(comparison_name, ".stat")))
  
  write_camera_results(
    category, list_of_gene_sets[[category]],
    comparison_name, species,
    de_res, camera_results[[category]], out_dir = out_dir)
  
  camera_results[[category]]
}

get_camera_results <- function(dds, gene_sets, gene_info) {
  vst <- dds %>% varianceStabilizingTransformation
  design_formula <- dds %>% design()
  expression_data <- vst %>% assay
  
  ids <- expression_data %>% 
    as.data.frame %>% 
    tibble::rownames_to_column(var = "gene") %>% 
    inner_join(gene_info) 
  
  idx <- gene_sets %>% ids2indices(id = ids$entrez_id)
  
  design_matrix <- model.matrix(design_formula, vst %>% colData)
  
  expression_data %>% camera(idx, design_matrix)
}

plot_gene_set <- function(results, gene_sets, gene_set_name, prefix) {
  idx <- gene_sets %>% 
    ids2indices(id = results$entrez_id) %>%
    extract2(gene_set_name)
  
  pval <- prefix %>% str_c(".pval") %>% rlang::sym()
  l2fc <- prefix %>% str_c(".l2fc") %>% rlang::sym()
  
  results %>% 
    mutate(signed_p = -log10(!!pval) * sign(!!l2fc)) %>%
    pull(signed_p) %>% 
    barcodeplot(index = idx, quantiles = c(-1,1)*(-log10(0.05)))
}

get_gene_set_results_matrix <- function(results, gene_sets, gene_set_name) {
  idx <- gene_sets %>%
    extract2(gene_set_name) %>%
    match(results$entrez_id) %>%
    na.omit
  
  genes_in_set <- results %>%
    extract(idx, ) %>%
    dplyr::select(gene) %>%
    group_by(gene) %>%
    filter(row_number() == 1) %>%
    ungroup
  
  genes_in_set[gene_set_name] = "T"
  
  results %>%
    dplyr::select(gene) %>%
    left_join(genes_in_set) %>%
    dplyr::select(-gene)
}

write_camera_results <- function(
  gene_set_collection_name, gene_set_collection, comparison_name, species, de_results, camera_results,
  barcodeplots = FALSE, fdr_cutoff = 0.1, out_dir="results/differential_expression/gsa/") {
  
  top_dir <- file.path(out_dir,species,comparison_name)
  if (!dir.exists(top_dir)) {
    dir.create(top_dir, recursive = TRUE)
  }
  
  camera_results %>% 
    tibble::rownames_to_column(var = "GeneSet") %>% 
    filter(FDR < fdr_cutoff) %>%
    write_csv(file.path(top_dir, str_c(comparison_name, "-", gene_set_collection_name, "_sets.csv")), na = "")
  
  ret <- list(enriched_sets = camera_results %>% 
                tibble::rownames_to_column(var="GeneSet") %>% 
                filter(FDR < fdr_cutoff))
  
  # Save only the significant gene sets
  camera_results %<>%
    tibble::rownames_to_column(var="GeneSet") %>%
    filter(FDR < fdr_cutoff)
  
  if ((camera_results %>% nrow) == 0) {
    return()
  }
  
  camera_results %>%
    extract2("GeneSet") %>%
    walk(function(x) {
      # We do not save a CSV file for each gene set now...
      # gene_set_results <- de_results %>% get_gene_set_results(gene_set_collection, x, str_c(comparison_name, ".pval"))
      # gene_set_results %>% write_csv(str_c(sub_dir, "/", x, ".csv"))
      
      if (barcodeplots) {
        sub_dir <- str_c(top_dir, "/", gene_set_collection_name)
        if (!dir.exists(sub_dir))
          dir.create(sub_dir, recursive = TRUE)
        
        start_plot(str_c(sub_dir, "/", x))
        plot_gene_set(de_results, gene_set_collection, x, comparison_name)
        end_plot()
      }
    })
  
  # Merge all GSA results into one file
  gene_set_names <- camera_results %>% extract2("GeneSet") %>% sort(method = "radix")
  gene_set_results <- gene_set_names %>%
    map_dfc(function(x) {
      de_results %>% get_gene_set_results_matrix(gene_set_collection, x)
    }) %>%
    setNames(gene_set_names)
  
  de_results %>%
    cbind(gene_set_results) %>%
    write_csv(file.path(top_dir, str_c(comparison_name, "-", gene_set_collection_name, "_genes_in_sets.csv")), na = "")
  
  ret[['genes_in_sets']] <- de_results %>% cbind(gene_set_results)
  
  ret
}

plot_significant_set_heatmap <- function(job_string,gs_result,list_of_gene_set,results_tbl,comparison_tbl,species,out_dir,...){
  
  futile.logger::flog.info(paste0('running heatmap for:',job_string))
  
  c(comparison_name, gs_category) %<-% (strsplit(job_string, split = ';', fixed = T) %>% unlist())

  significant_gs <- gs_result[[job_string]] %>% filter(FDR < 0.05) %>% rownames()
  significant_gs_entrezs <- list_of_gene_set[[gs_category]][significant_gs]
  
  #get the comparison criteria for this comparison
  comparison_criteria <- comparison_tbl %>% filter(comparison == comparison_name) %>% dplyr::select(condition_name, condition, condition_base, filter)
  selected_conditions <- comparison_criteria %>% dplyr::select(condition, condition_base) %>% as.character()
  
  # pull out samples which have the conditions in this comparison
  # and samples which match the filter specified
  condition_name <- comparison_criteria %>% dplyr::select(condition_name) %>% pull()
  samples_in_comparison <- SAMPLE_DATA %>%
    filter(!!parse_expr(condition_name) %in% selected_conditions &
             !!(parse_expr(comparison_criteria %>% dplyr::select(filter) %>% pull()))) %>%
    dplyr::select(sample_name, !!parse_expr(condition_name)) %>%
    arrange(!!parse_expr(condition_name))
  
  # samples_in_comparison %<>% mutate(fpkm_columns = str_c(sample_name, "_fpkm"))
  # fpkm_columns <- samples_in_comparison %>% pull(fpkm_columns)
  # num_samples <- length(fpkm_columns)
  
  # get names of FPKM data columns
  samples_in_comparison %<>% mutate(fpkm_columns = str_c(sample_name, "_fpkm"))
  fpkm_columns <- samples_in_comparison %>% pull(fpkm_columns)
  num_samples <- length(fpkm_columns)
  
  # construct sample annotation
  rownames(samples_in_comparison) <- NULL
  samples_in_comparison %<>% tibble::column_to_rownames(var = "fpkm_columns")
  annot <- samples_in_comparison %<>% dplyr::select(-sample_name)
  
  futile.logger::flog.info(paste0('Number of gene set to be plot: ', significant_gs_entrezs%>%names%>%length))
  lapply(significant_gs_entrezs%>%names,function(gs_name){
    futile.logger::flog.info(paste0('running heatmap for:',job_string,'/',gs_name))
    # get the significant entrez ids
    entrez_ids <- significant_gs_entrezs[[gs_name]]
    # get the name of the column of the log 2 fold changes
    comparison_log2fc <- paste(comparison_name, "l2fc", sep = '.')
    
    # from the global results variable, get rows and columns corresponding to significant entrez IDs
    # and samples from the correct comparison; remove genes with no name
    to_heatmap <- results_tbl %>%
      dplyr::select(gene_name, entrez_id, fpkm_columns, comparison_log2fc) %>%
      filter(entrez_id %in% entrez_ids) %>%
      filter(!is.na(gene_name))
    
    # reorder
    to_heatmap <- to_heatmap %>% arrange(desc(!!sym(comparison_log2fc)))
    
    # declare the path to the heatmaps, based on the gene set category and comparison
    # create the dir if it doesnt exist
    heatmap_path <- file.path(out_dir, species, comparison_name, gs_category )
    ifelse(!dir.exists(file.path(heatmap_path)), dir.create(file.path(heatmap_path), recursive = TRUE), FALSE)
    
    # TODO: currently, the column to rownames call complains about duplicate row names (i.e. gene names)
    # I have removed duplicates - is this how we want to do this? Is there a better way?
    to_heatmap_unique <- distinct(to_heatmap, gene_name, .keep_all = TRUE)
    
    # prepare heatmap data so the row names are the gene IDs and remove the entrez column
    heatmap_data <- to_heatmap_unique %>% 
      tibble::column_to_rownames(var="gene_name") %>% 
      dplyr::select(-entrez_id, -comparison_log2fc)
    
    # divide each row by the mean of that row
    heatmap_data <- t(apply(heatmap_data, 1, function(x) x/mean(x)))
    heatmap_data <- log2(heatmap_data)
    
    heatmap_data[is.infinite(heatmap_data)] <- NA
    heatmap_data[is.nan(heatmap_data)] <- NA
    heatmap_data <- subset(heatmap_data,rowSums(!is.na(heatmap_data)) == nrow(samples_in_comparison))
    
    if(nrow(heatmap_data)!=0){
      max_data <- max(heatmap_data, na.rm = TRUE)
      min_data <- -min(heatmap_data, na.rm = TRUE)
      range <- min(max_data, min_data)
      
      start_plot(prefix = gs_name, path = heatmap_path)
      pheatmap(heatmap_data,
               breaks = seq(-range, range, length.out = 100),
               cluster_rows = FALSE, cluster_cols = FALSE,
               border_color = NA, show_rownames = (heatmap_data %>% nrow()) < 100,
               annotation_col = annot)
      end_plot()
    }
  })
}

#' #' This function reads in all the GSA results and tracks the given gene sets
#' #' @example:
#' #' track_gene_sets(target_terms=c('GO_RIBOSOME','GO_PROTEASOME_COMPLEX','GO_TRANSLATIONAL_INITIATION'), category='GO',comparison_table=COMPARISON_TABLE)
#' track_gene_sets <- function(target_terms = c('GO_RIBOSOME', 'GO_PROTEASOME_COMPLEX', 'GO_TRANSLATIONAL_INITIATION'), 
#'                             category = 'GO',
#'                             gs_results = get_global('GS_results'),
#'                             comparison_table = COMPARISON_TABLE,
#'                             print_table = TRUE, 
#'                             output_table_file = NA,
#'                             left_out_comparison = c(),
#'                             heat_map.fdr.midpoint = 0.05){
#'   
#'   ## create a master table for all comparisons and the gene sets results
#'   gsa_res_tb <- COMPARISON_TABLE %>% 
#'     pull(comparison) %>% 
#'     extract(which(!. %in% left_out_comparison)) %>% 
#'     set_names(.) %>%
#'     sapply(simplify = FALSE, USE.NAMES = TRUE, function(x) {
#'       gs_results %>% 
#'         extract2(x) %>% 
#'         extract2(category) %>%
#'         tibble::rownames_to_column('GeneSet') %>% 
#'         dplyr::mutate(comparison = x)
#'     }) %>% 
#'     reduce(rbind)
#'   
#'   ## we only keep gene sets of interests
#'   gsa_res_tb <- gsa_res_tb %>% 
#'     filter(GeneSet %in% target_terms) %>%
#'     mutate(log10fdr = log10(FDR)) # for plotting
#'   
#'   gsa_res_tb$comparison %<>% factor(levels = COMPARISON_TABLE %>% pull(comparison) %>% rev)
#'   gsa_res_tb$GeneSet %<>% factor(levels = target_terms)
#'   
#'   if (print_table) {
#'     print(gsa_res_tb %>% arrange(GeneSet, FDR))
#'   }
#'   
#'   # we output the FDR table to csv
#'   if (!is.na(output_table_file)) {
#'     gsa_res_tb %>% 
#'       arrange(GeneSet, FDR) %>% 
#'       dplyr::select(GeneSet, FDR, comparison) %>%
#'       reshape(idvar = "comparison", timevar = "GeneSet", direction = "wide") %>%
#'       write.csv(file = output_table_file)
#'   }
#'   
#'   gsa_res_tb %>% 
#'     ggplot(aes(GeneSet, comparison)) +
#'     geom_tile(aes(fill = log10fdr)) + 
#'     scale_fill_gradient2(low = "red", high = "white", mid = "white",
#'                          midpoint = log10(heat_map.fdr.midpoint)) +
#'     geom_text(aes(label = if_else(FDR < heat_map.fdr.midpoint, 
#'                                   ifelse(Direction == 'Up', "UP", "DOWN"),
#'                                   "")), alpha = 0.75, size = 3) +
#'     labs(fill = "log10(FDR)") +
#'     labs(title = str_c('FDR cutoff = ', heat_map.fdr.midpoint),
#'          xlab = "Gene Set", ylab = "Comparison") + 
#'     theme_minimal() + 
#'     theme(axis.text.x = element_text(angle = -90, size=7), 
#'           panel.border = element_blank(), panel.background = element_blank())
#' }

# track_go <- function(target_terms = c('GO:0051492', 'GO:0010811'),
#                      category = 'BP',
#                      gs_results = get_global('GO_results'),
#                      comparison_table = COMPARISON_TABLE,
#                      print_table = TRUE,
#                      output_table_file = NA,
#                      left_out_comparison = c(),
#                      heat_map.p.midpoint = 0.05){
#   
#   ## create a master table for all comparisons and the gene sets results
#   gsa_res_tb <- COMPARISON_TABLE %>%
#     pull(comparison) %>%
#     extract(which(!. %in% left_out_comparison)) %>%
#     set_names(.) %>%
#     sapply(simplify = FALSE, USE.NAMES = TRUE, function(x) {
#       gs_results %>%
#         extract2(x) %>%  extract2(str_c(x,'.all')) %>%
#         extract2(category) %>%  extract2('go_results') %>%
#         dplyr::mutate(comparison = x)
#     }) %>%
#     reduce(rbind)
#   
#   ## we only keep gene sets of interests
#   gsa_res_tb <- gsa_res_tb %>%
#     filter(GO.ID %in% target_terms) %>%
#     mutate(log10p=log10(as.numeric(weight_fisher))) # for plotting
#   
#   gsa_res_tb$comparison %<>% factor(levels = COMPARISON_TABLE %>% pull(comparison) %>% rev)
#   gsa_res_tb$GO.ID %<>% factor(levels = target_terms)
#   
#   if (print_table) {
#     print(gsa_res_tb %>% arrange(GO.ID, log10p))
#   }
#   
#   # we output the FDR table to csv
#   if (!is.na(output_table_file)) {
#     gsa_res_tb %>%
#       arrange(GO.ID, log10p) %>%
#       dplyr::select(GO.ID, weight_fisher, comparison) %>%
#       reshape(idvar = "comparison", timevar = "GO.ID", direction = "wide") %>%
#       write.csv(file = output_table_file)
#   }
#   
#   gsa_res_tb %>%
#     ggplot(aes(GO.ID, comparison)) +
#     geom_tile(aes(fill = log10p)) +
#     scale_fill_gradient2(low = "red", high = "white", mid = "white",
#                          midpoint = log10(heat_map.p.midpoint)) +
#     # geom_text(aes(label = if_else(log10p < log10(heat_map.p.midpoint),
#     #                               ifelse(Direction=='Up', "UP", "DOWN"),
#     #                               "")), alpha = 0.75, size=3) +
#     labs(fill ="log10(p)") +
#     labs(title=str_c('p cutoff = ',heat_map.p.midpoint),
#          xlab="GO.ID", ylab="Comparison") +
#     theme_minimal() +
#     theme(axis.text.x = element_text(angle = -90, size=7),
#           panel.border = element_blank(),panel.background = element_blank())
# }