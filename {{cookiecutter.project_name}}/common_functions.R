KALLISTO <- "kallisto"
SALMON <- "salmon"

GO_TERMS <- Term(GOTERM) %>% as.data.frame() %>% tibble::rownames_to_column(var="GO.ID")
colnames(GO_TERMS) <- c("GO.ID", "FullTerm")

#### Utility functions ####

set_global <- function(value, variable) {
  assign(variable, envir=.GlobalEnv, value)
}

get_global <- function(variable) {
  get(variable, envir = .GlobalEnv)  
}

global_exists <- function(variable) {
  exists(variable, where = .GlobalEnv)
}

rm_global <- function(variable) {
  rm(list=c(variable), envir = .GlobalEnv)
}

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
  count_data %>% tibble::column_to_rownames(var="gene")
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

get_gene_lengths <- function(species) {
  species %>% 
    str_c("data/", ., "_ensembl_{{cookiecutter.ensembl_version}}/gene_lengths.csv") %>% 
    read_csv
}

check_formulas <- function() {
  for (r in COMPARISON_TABLE %>% rownames()) {
    row <- COMPARISON_TABLE[r,]
    f <- row$formula %>% as.formula() %>% terms()
    condition <- row$condition_name

    # if there is an interaction, we check if we have the interaction detail in the INTERACTION_TABLE,
    # otherwise, we check if we have the deciding condition correctly setup
    if(any(f %>% attr( "order") > 1)){
      if(! row$comparison %in% INTERACTION_TABLE$comparison){
        print(row)
        stop("interaction formula was detected but interaction details cannot be found in the INTERACTION_TABLE.")
      }
    }else{
      deciding_condition <- labels(f) %>% tail(1)

      if (condition != deciding_condition) {
        print(row)
        stop("The formula ends with a label which is different to the one specified in the condition_name column.
           This will cause the GSA algorithm to pick up the wrong condition.")
      }
    }
  }
}

check_samples <- function(){
    COMPARISON_TABLE %>% pull(comparison) %>% lapply(function(comparison_name){
        x <- COMPARISON_TABLE %>% filter(comparison==comparison_name)
        sample_data<-SAMPLE_DATA %>% filter(!!parse_expr(x$filter))
        list(
        Comparison = x$comparison,
        Total_samples = sample_data %>% nrow(),
        base_level_condition = x$condition_base,
        num_base = sample_data %>%
            filter(!!parse_expr(x$condition_name) == x$condition_base) %>%
            nrow(),
        base_samples = sample_data %>%
            filter(!!parse_expr(x$condition_name) == x$condition_base) %>%
            pull(sample_name) %>%
            str_c(collapse = ','),
        Comparison_level_condition = x$condition,
        comparison_samples = sample_data %>%
            filter(!!parse_expr(x$condition_name) == x$condition) %>%
            nrow(),
        Sample_names_in_comparison_level_condition = sample_data %>%
            filter(!!parse_expr(x$condition_name) == x$condition) %>%
            pull(sample_name) %>%
            str_c(collapse = ',')
        ) %>% as.data.frame()
    }) %>% reduce(rbind)
}

check_comparison_name <- function(){
    invalid=COMPARISON_TABLE %>% pull(comparison) %>% grepl('^[0-9]',x = .,perl=T)
    if(any(invalid)){
        message('Comparison cannot start with a number. The following comparison names are invalid: ')
        cat(COMPARISON_TABLE %>% pull(comparison) %>% extract(which(not_valid))%>%paste(collapse = '\n'))
    }
}

start_plot <- function(prefix,width=12, height=12, path=GRAPHS_DIR, num_plots=1) {
  .adjust_pdf_size<-function(num_plots){
    num_features <- num_plots
    num_row<-sqrt(num_features) %>% ceiling()
    num_column<-num_features/num_row %>% ceiling()
    c('width'= max(num_row/2,1),'height'=max(num_column/2,1))
  }

  sf <- .adjust_pdf_size(num_plots)

  if (PLOT_TO_FILE) {
    prefix %>%
      str_c(path, ., "_", SPECIES, ".pdf") %>%
      pdf(width=width*sf['width'], height=height*sf['height'])
  }
}

end_plot <- function() {
  if (PLOT_TO_FILE) {
    dev.off() 
  }
}

add_to_patchwork <- function(plot2add, plot_var_name='pathworkplot') {
  if (global_exists(plot_var_name)) {
    p <- get_global(plot_var_name)
    p <- p + plot2add
  } else {
    p <- plot2add
  }
  
  p %>% set_global(plot_var_name)
}

start_parallel <- function(cores=NA) {
  if (is.na(cores)) {
    cores <- get_global('NUM_CORES')
  }
  
  options("mc.cores" = cores)
  set_global(TRUE, "PARALLEL")
}

stop_parallel <- function() {
  options("mc.cores" = 1L)
  set_global(FALSE, "PARALLEL")
}

adjust_parallel_cores <- function() {
  currect_cores <- getOption("mc.cores", get_global('NUM_CORES'))
  reduced_cores <- floor(currect_cores/3)
  
  if (reduced_cores > 10) {
    reduced_cores = 10
  }
  
  options("mc.cores" = reduced_cores)
}

lapply_fork <- function(X, FUN, cores = NA) {
  # This is the fork approach of parallel lapply:
  # http://dept.stat.lsa.umich.edu/~jerrick/courses/stat701/notes/parallel.html#starting-a-cluster

  if (get_global('PARALLEL')) {
    if (is.na(cores)) {
      cores <- getOption("mc.cores", get_global('NUM_CORES'))
    }

    res <- mclapply(mc.cores = cores, X = X, FUN = FUN)
  } else {
    #run the normal lapply in single core
    res <- lapply(X = X, FUN = FUN)
  }

  #check error
  error_index <- sapply(res,function(x){
    inherits(x,'try-error')
  }) %>% which()

  if(length(error_index) > 0 ){
    error_msg<-''
    for(i in error_index){
      error_msg <- str_c(error_msg,
      str_c("Error from node ",i, ":"),
      res[[i]])
    }

    stop(str_c("\n",error_msg),call. = F)
  }
  
  res
}

lapply_socket <- function(X, FUN, cores = NA, export_objects=c("expressed_genes", "results", "PARALLEL","META_DATA", "OUTPUT_DIR")) {
  # This is the socket approach of parallel lapply:
  # http://dept.stat.lsa.umich.edu/~jerrick/courses/stat701/notes/parallel.html#starting-a-cluster

  if (get("PARALLEL", envir = .GlobalEnv)) {
    # create cluster
    if (is.na(cores)) {
      cores <- getOption("mc.cores", get('NUM_CORES', envir = .GlobalEnv))
    }

    cl <- makeCluster(cores)

    # export variables
    if (length(export_objects) > 0) {
      clusterExport(cl, varlist = export_objects)
      clusterExport(cl, varlist = c('export_objects'), envir = environment())
    }

    clusterEvalQ(cl, {
      source(get('META_DATA',envir = environment()))
    })

    # run the parallel code
    ret <- parSapply(cl = cl, X = X, simplify = FALSE, USE.NAMES = TRUE, FUN=FUN)

    # stop cluster
    stopCluster(cl)

    ret
  } else {
    # run the normal lapply in single core
    sapply(X = X, simplify = FALSE, USE.NAMES = TRUE, FUN=FUN)
  }
}

load_rs_data <- function(file = 'results/Rworkspace/diff_expr.RData') {
  if (!file.exists(file)) {
    stop(file,' does not exist!')
  }
  
  load(file, envir = .GlobalEnv)
}

#### Gene-level D. E. analyses ####

get_total_dds <- function(sample_data, species, filter_low_counts = FALSE, qSVA = FALSE, design_formula = ~1) {
  # Collate count data
  total_count_data <- sample_data %>%
    row.names() %>%
    map(read_counts,species) %>%
    purrr::reduce(inner_join) %>%
    remove_gene_column()
  
  get_deseq2_dataset(
    total_count_data, sample_data,
    filter_low_counts = filter_low_counts, design_formula = ~1, qSVA = qSVA)
}

get_deseq2_dataset <- function(count_data, sample_data, filter_low_counts = TRUE, 
                               design_formula = ~condition, qSVA = FALSE) {
  
  dds <- DESeqDataSetFromMatrix(countData=count_data, colData=sample_data, design=design_formula)
  
  if (filter_low_counts) {
    dds <- dds[rowSums(counts(dds)) > 1, ]
  }

  if(qSVA) {
    dds %<>% get_qsva_dds()
  }

  betaPrior=TRUE
  if(any(dds %>% design() %>% terms.formula %>% attr( "order") >1)) betaPrior=FALSE

  dds %>% DESeq(betaPrior = betaPrior)
}

get_deseq2_results <- function(dds, comparison, condition, condition_base, alpha = 0.05) {

  # if we have interaction term in the formula, we ignore the condition parameter and check the last
  # contrast, which in a 2 by 2 condition case, will be the testing
  # if the condition effect is different in interaction term 1 compared to 2
  if(any(dds %>% design() %>% terms.formula %>% attr( "order")>1)){
    res <- dds %>% results(name=resultsNames(dds) %>% tail(1), alpha=alpha)
  }else{
    res <- dds %>% results(c(comparison, condition, condition_base), alpha=alpha)
  }

  res %>% summary %>% print
  
  res %>% as.data.frame() %>%
    tibble::rownames_to_column(var = "gene") %>%
    dplyr::select(-baseMean, -lfcSE)
}

get_raw_l2fc <- function(dds, sample_data, comparison_samples_filter) {
  sample_data %<>% tibble::rownames_to_column(var = "tmp_sample_name")
  
  all_samples <- sample_data %>% extract2("tmp_sample_name") %>% as.vector

  comparison_sample_data <- filter(sample_data, !!comparison_samples_filter)
  comparison_samples <- comparison_sample_data %>% extract2("tmp_sample_name") %>% as.vector

  base_samples = all_samples %>% setdiff(comparison_samples)
  
  dds %>% 
    get_count_data %>% 
    mutate(comparison_mean = rowMeans(.[comparison_samples]),
           base_mean = rowMeans(.[base_samples]),
           raw_l2fc = log2(comparison_mean / base_mean)) %>%
    dplyr::select(gene, raw_l2fc)
}

get_deseq2_results_name <- function(dds, name, alpha = 0.05) {
  res <- dds %>% results(name = name, alpha = alpha)
  res %>% summary %>% print
  
  res %>% as.data.frame() %>%
    tibble::rownames_to_column(var = "gene") %>%
    dplyr::select(-baseMean, -lfcSE, -stat)
}

# Get differential expression results for a given comparison name
get_res <- function(comparison_name, tpms, species, qSVA=FALSE, 
                    use_tx=FALSE, quant_method='salmon', tx_level=FALSE, alpha=0.05) {
  
  comparison_table <- COMPARISON_TABLE
  x=comparison_table %>% filter(comparison == comparison_name)
  sample_data <- SAMPLE_DATA %>%
    tibble::rownames_to_column(var = "tmp_row_names") %>%
    mutate(!!x$condition_name:= factor(!!parse_expr(x$condition_name))) %>%
    filter(!!parse_expr(x$filter)) %>%
    mutate(sample_name_tmp = tmp_row_names) %>%
    tibble::column_to_rownames(var = "tmp_row_names")
  
  # Ensure that conditions to be used in GSA comparisons are factors with the correct base level set.
  if(any(x$formula %>% as.formula() %>% terms.formula %>% attr( "order") >1)){
    # if there is interaction, we relevel the interaction conditions
    x2 <- INTERACTION_TABLE %>% filter(comparison == comparison_name)
    sample_data[,x2$condition_name_1] %<>% as.factor() %>% relevel(x2$condition_base_1)
    sample_data[,x2$condition_name_2] %<>% as.factor() %>% relevel(x2$condition_base_2)
  }else{
    sample_data[,x$condition_name] %<>% relevel(x$condition_base)
  }
  
  if (use_tx) {
    txi <- get_tximport(sample_data, quant_method, tx_level)
    dds <- DESeqDataSetFromTximport(txi, sample_data, x$formula %>% as.formula())
    dds <- dds[rowSums(counts(dds)) > 1, ]
    dds <- DESeq(dds, betaPrior = TRUE)
  }else{
    dds <- sample_data %>%
      row.names() %>%
      map(read_counts, species) %>%
      purrr::reduce(inner_join) %>%
      remove_gene_column() %>%
      get_deseq2_dataset(sample_data, design_formula = x$formula %>% as.formula(), qSVA=qSVA)
  }
  
  if (use_tx & tx_level) {
    # res contains transcripts, rather than genes
    res <- dds %>%
      get_deseq2_results(x$condition_name, x$condition, x$condition_base, alpha=alpha) %>%
      left_join(dds %>% get_raw_l2fc(sample_data, expr(!!sym(x$condition_name) == !!(x$condition)))) %>%
      dplyr::rename(transcript = gene)
  } else {
    res <- dds %>%
      get_deseq2_results(x$condition_name, x$condition, x$condition_base, alpha = alpha) %>%
      left_join(dds %>% get_raw_l2fc(sample_data, expr(!!sym(x$condition_name) == !!(x$condition))))
  }
  
  # Create a PCA plot and heatmap for this comparison
  if (!use_tx) {
    vst <- dds %>% varianceStabilizingTransformation
    vst %>% plot_pca(intgroup = c(x$condition_name)) %>% print
    vst %>% plot_heat_map(sample_data %>%
                          mutate(sample_info = str_c(!!parse_expr(x$condition_name), 
                                                     sample_name_tmp, sep = ":")) %>%
                          extract2("sample_info"))
    
    start_plot(str_c("pca_", x$comparison))
    vst %>% plot_pca(intgroup=x$condition_name,output_data_table_path=file.path(GRAPHS_DIR,str_c("pca_", x$comparison,'_',SPECIES,'.csv'))) %>% print()
    end_plot()

    start_plot(str_c("heatmap_", x$comparison))
    vst %>% plot_heat_map(sample_data %>%
                            tidyr::unite(col = 'sample_info',c(sample_name_tmp,x$condition_name), 
                                         sep = ":", remove = FALSE) %>%
                            extract2("sample_info"))
    end_plot()
  }
  
  summary_tb_row <- list(
    Comparison = x$comparison, 
    DESeq_model_formula = design(dds) %>% format(),
    Condition_tested = x$condition_name,
    Total_number_of_samples_data = sample_data %>% nrow(),
    Base_level_condition = x$condition_base,
    Number_of_samples_in_base_level_condition = sample_data %>% 
      filter(!!parse_expr(x$condition_name) == x$condition_base) %>% 
      nrow(),
    Sample_names_in_base_level_condition = sample_data %>% 
      filter(!!parse_expr(x$condition_name) == x$condition_base) %>% 
      pull(sample_name_tmp) %>% 
      str_c(collapse = ','),
    Comparison_level_condition = x$condition,
    Number_of_samples_in_comparison_level_condition = sample_data %>% 
      filter(!!parse_expr(x$condition_name) == x$condition) %>% 
      nrow(),
    Sample_names_in_comparison_level_condition = sample_data %>% 
      filter(!!parse_expr(x$condition_name) == x$condition) %>% 
      pull(sample_name_tmp) %>% 
      str_c(collapse = ','),
    p.adj.cutoff = P.ADJ.CUTOFF,
    Up_regulated = res %>% 
      filter(padj < P.ADJ.CUTOFF & log2FoldChange > 0) %>% 
      nrow(),
    Down_regulated = res %>% 
      filter(padj < P.ADJ.CUTOFF & log2FoldChange < 0) %>% 
      nrow(),
    D.E.total = res %>% 
      filter(padj < P.ADJ.CUTOFF) %>% 
      nrow()
  )
  
  list(res=res, dds=dds, summary_tb_row=summary_tb_row)
}

get_count_data <- function(dds, norm = T) {
  dds %>%
    counts(normalized = norm) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "gene")
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
  
  fpkms <- all_counts %>% dplyr::select(gene, dplyr::contains(col_suffix))
  
  # work out average fpkm
  fpkms %>% left_join(get_avg_fpkm(fpkms))
}

.get_avg_fpkm_table <- function(use_tpm=FALSE){
  unit <- ifelse(use_tpm,'_tpm','_fpkm')
  lapply(AVG_FPKM_GROUP,function(g){
    SAMPLE_DATA %>%  tibble::rownames_to_column(var = "tmp_row_names") %>%
      group_by(.dots=g) %>%
      dplyr::summarise(samples=str_c(tmp_row_names, unit, sep = '', collapse = ',')) %>%
      tidyr::unite('avg_name', g, sep='_')
  }) %>% reduce(rbind)
}

get_avg_fpkm <- function(fpkms,use_tpm=FALSE) {
  avg_table <- .get_avg_fpkm_table(use_tpm)
  avg_postfix <- ifelse(use_tpm,'_avg_tpm','_avg_fpkm')
  for (avg in avg_table$avg_name) {
    samples <- avg_table %>%
      filter(avg_name == avg) %>%
      pull(samples) %>%
      str_split(',') %>%
      unlist()
    
    avg_fpkm <- fpkms %>% dplyr::select(one_of(samples)) %>%
      mutate(avg=rowMeans(.)) %>%
      dplyr::pull(avg)
    
    fpkms %<>% mutate(!!str_c(avg,avg_postfix,sep='') := avg_fpkm)
  }
  
  fpkms %>% dplyr::select(gene, contains('avg'))
}

get_avg_tpm <- function(tpms, tx_level) {
  avg_tpm <- get_avg_fpkm(tpms,T)
  avg_tpm %>% dplyr::select(gene, contains('avg'))
}

save_results_by_group <- function(results,use_tx=FALSE) {
  # This function will, for each group defined in the COMPARISON_TABLE, save the comparisons into
  # a CSV file, prefixed with the group name.
  for (g in COMPARISON_TABLE %>% pull(group) %>% unique()) {
    
    comparisons <- COMPARISON_TABLE %>% 
      filter(group == g) %>% 
      pull(comparison)
    
    n_comparisons <- COMPARISON_TABLE %>% 
      filter(group != g) %>% 
      pull(comparison) %>% 
      str_c("^", ., '\\.' ,collapse = '|')
    
    samples_to_include <- SAMPLE_DATA %>%
      filter(!!parse_expr(COMPARISON_TABLE %>% filter(group==g) %>% pull(filter) %>% str_c(collapse = '|'))) %>%
      pull(sample_name) %>% as.vector()
    
    samples_to_exclude <- SAMPLE_DATA %>%
      filter(! sample_name %in% samples_to_include) %>%
      pull(sample_name) %>% as.vector()
    
    samples_to_exclude_pattern <- samples_to_exclude %>% str_c('^',.,sep='',collapse = '|')
    
    ## work out what avg column to be exclude for group
    avg_tb<-.get_avg_fpkm_table(use_tx)
    avg_to_exclude<- lapply(avg_tb$avg_name,function(x){
      samples_in_avg <- avg_tb %>% filter(avg_name==x) %>% pull(samples) %>% strsplit(',') %>% extract2(1) %>% gsub('_fpkm|_tpm','',x=.)
      (samples_in_avg %in% samples_to_include) %>% all
    }) %>% unlist() %>% `!` %>% filter(.data=avg_tb) %>% pull(avg_name)
    
    if(use_tx){
      
      if ('transcript' %in% colnames(results)) {
        columns_included <- c('transcript', 'transcript_length', 'gene', 'number_of_transcript')
        tx_level_str <- "transcript"
      } else {
        columns_included <- c('gene','gene_length', 'max_transcript_length')
        tx_level_str <- "gene"
      }
      
      results %>%
        dplyr::select(
          columns_included, gene_name, chromosome, description, entrez_id, gene_type, everything(),
          -dplyr::contains("_tpm"), -dplyr::ends_with(".stat"),
          -matches(n_comparisons), -(samples_to_exclude)) %>%
        write_csv(file.path(OUTPUT_DIR, "de_gene", str_c(g, "_deseq2_results_count_", SPECIES, "_tx_",tx_level_str,'_',QUANT_METHOD,".csv")), na="")
      
      tpm_output <- results %>%
        dplyr::select(
          columns_included, gene_name, chromosome, description, entrez_id, gene_type,
          dplyr::contains("_tpm"),
          COMPARISON_TABLE %>%
            pull(comparison) %>%
            sapply(FUN = function(x) results %>% colnames() %>% str_which(str_c("^", x, sep =''))) %>% 
            unlist() %>%
            as.vector() %>%
            unique(),
          -dplyr::ends_with(".stat"), -matches(n_comparisons)
        )
      
      if (length(samples_to_exclude) > 0) {
        tpm_output %<>% dplyr::select(-matches(samples_to_exclude_pattern))
      }
      
      if (length(avg_to_exclude) > 0) {
        tpm_output %<>% dplyr::select(-one_of(avg_to_exclude %>% str_c('_avg_tpm')))
      }
      
      tpm_output %>% write_csv(file.path(OUTPUT_DIR, "de_gene", str_c(g ,"_deseq2_results_tpm_", SPECIES, "_tx_",tx_level_str,'_',QUANT_METHOD,".csv")), na="")
      
      SUMMARY_TB %>% filter(Comparison %in% comparisons) %>%
        write_csv(file.path(OUTPUT_DIR, "de_gene", str_c(g ,"_de_summary_", SPECIES, "_tx_",tx_level_str,'_',QUANT_METHOD,".csv")), na="")
    }else{
      # non-tx
      results %>%
        dplyr::select(
          gene, gene_name, chromosome, description, entrez_id, gene_type,
          gene_length, max_transcript_length, everything(),
          -dplyr::contains("_fpkm"), -dplyr::ends_with(".stat"),
          -matches(n_comparisons), -(samples_to_exclude)) %>%
        write_csv(file.path(OUTPUT_DIR, "de_gene", str_c(g, "_deseq2_results_count_", SPECIES, ".csv")), na="")

      fpkm_output <- results %>%
        dplyr::select(
          gene, gene_name, chromosome, description, entrez_id, gene_type,
          gene_length, max_transcript_length,
          dplyr::contains("_fpkm"),
          COMPARISON_TABLE %>%
            pull(comparison) %>%
            sapply(FUN = function(x) results %>% colnames() %>% str_which(str_c("^", x, sep =''))) %>% 
            unlist() %>%
            as.vector() %>%
            unique(),
          -dplyr::ends_with(".stat"), -matches(n_comparisons)
        ) 
      
      if (length(samples_to_exclude) > 0) {
        fpkm_output %<>% dplyr::select(-matches(samples_to_exclude_pattern))
      }
      
      if (length(avg_to_exclude) > 0) {
        fpkm_output %<>% dplyr::select(-one_of(avg_to_exclude %>% str_c('_avg_fpkm')))
      }
      
      fpkm_output %>% write_csv(file.path(OUTPUT_DIR, "de_gene", str_c(g ,"_deseq2_results_fpkm_", SPECIES, ".csv")), na="")
      
      SUMMARY_TB %>% filter(Comparison %in% comparisons) %>%
        write_csv(file.path(OUTPUT_DIR, "de_gene", str_c(g ,"_de_summary_", SPECIES, ".csv")), na="")
    }
  }
}

read_de_results <- function(filename, num_samples, num_conditions, num_comparisons, extra_columns = "") {
  col_types_string <- str_c(
    "ccccii",
    strrep("d", num_samples + num_conditions + num_comparisons * 4),
    extra_columns)
  print(nchar(col_types_string))
  
  read_csv(filename, col_types = col_types_string)
}

#### Quality control checks and plots ####

plot_pca <- function(vst, intgroup=c("condition"), plot_label = TRUE, label_name='name', include_gene = c(),
                     removeBatchEffect = FALSE, batch = NULL, output_data_table_path=NULL){
  
  if (removeBatchEffect) {
    if (is.null(batch)) {
      stop('batch cannot be NULL.')
    }
    
    assay(vst) <- limma::removeBatchEffect(assay(vst), vst %>% extract2(batch))
  }
  
  pca_data <- vst %>% plotPCA2(intgroup = intgroup, returnData = TRUE, include_gene = include_gene)

  # we want to save the data table used for the pca plot for future reference https://github.com/sidbdri/cookiecutter-de_analysis_skeleton/issues/115
  if(!is.null(output_data_table_path)){
    write.csv(pca_data, output_data_table_path)
  }
  
  percent_var <- round(100 * attr(pca_data, "percentVar"))
  
  intgroup.df <- as.data.frame(colData(vst)[, intgroup, drop = FALSE])
  
  if (length(intgroup) > 2) {
    colour_group <- factor(apply(intgroup.df[-ncol(intgroup.df)], 1, paste, collapse = " : "))
    shape_group <- factor(intgroup.df[[ncol(intgroup.df)]])
    p <- pca_data %>% ggplot(aes(PC1, PC2, color=colour_group, shape=shape_group)) + 
      guides(colour=guide_legend(title="group"),
             shape=guide_legend(title=intgroup[length(intgroup)]))    
  }
  else if (length(intgroup) == 2) {
    colour_group <- factor(intgroup.df[[1]])
    shape_group <- factor(intgroup.df[[2]])
    p <- pca_data %>% ggplot(aes(PC1, PC2, color=colour_group, shape=shape_group)) + 
      guides(colour=guide_legend(title=intgroup[1]),
             shape=guide_legend(title=intgroup[2]))
  } else {
    group <- colData(vst)[[intgroup]]
    p <- pca_data %>% ggplot(aes(PC1, PC2, color=group)) +
      guides(color=guide_legend(title=intgroup))
  }
  
  p <- p + geom_point(size=3) +
    xlab(str_c("PC1: ", percent_var[1], "% variance")) +
    ylab(str_c("PC2: ", percent_var[2], "% variance")) +
    theme(legend.position = "right")
  
  if(plot_label)
    p <- p + geom_text(aes(label = !!parse_expr(label_name)), colour="darkgrey", 
                       position=position_nudge(y = 1), size=3)

  p
}

plotPCA2 <- function(object, ...) {
    # This function is a hack of the plotPCA function from DESeq2 package. Instead of using all the genes 
    # for the PCA plot, this function accepts a parameter include_gene=c() which filters the genes in the
    # expression array. This can be use to make a PCA plot for only a subset of gene of interests.
    .local <- function(object, intgroup = "condition", ntop = 500,
                       returnData = FALSE, include_gene=c()) {
      
        count_data <- assay(object)

        if (length(include_gene) > 0) {
            include_gene <- include_gene[which(include_gene %in% (count_data %>% rownames()))]
            count_data <- count_data[include_gene,]
        }

        rv <- rowVars(count_data)
        select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
        pca <- prcomp(t(count_data[select, ]))
        percentVar <- pca$sdev^2/sum(pca$sdev^2)
        
        if (!all(intgroup %in% names(colData(object)))) {
            stop("the argument 'intgroup' should specify columns of colData(dds)")
        }
        
        intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
        group <- if (length(intgroup) > 1) {
            factor(apply(intgroup.df, 1, paste, collapse = ":"))
        }else {
            colData(object)[[intgroup]]
        }
        
        d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group, intgroup.df, name = colnames(object))
        
        if (returnData) {
            attr(d, "percentVar") <- percentVar[1:2]
            return(d)
        }
        
        ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group")) +
            geom_point(size = 3) + xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
            ylab(paste0("PC2: ", round(percentVar[2] *  100), "% variance")) + coord_fixed()
    }
    .local(object,  ...)
}

plot_heat_map <- function(vst, sample_names) {
  distsRL <- vst %>% assay %>% t %>% dist
  
  mat <- distsRL %>% as.matrix()
  rownames(mat) <- colnames(mat) <- sample_names
  
  hc <- distsRL %>% hclust
  hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
  heatmap.2(mat, Rowv = hc %>% as.dendrogram, 
            symm = TRUE, trace = "none",
            col = hmcol %>% rev, margin = c(10, 10))
}

plot_count_distribution <- function(dds, norm = T) {
  counts <- dds %>% 
    get_count_data(norm = norm) %>% 
    melt(id.vars=c("gene"), variable.name = "sample", value.name = "count") 
  
  p <- ggplot(counts, aes(sample, 1 + count)) + 
    geom_violin(aes(fill = sample), scale="width") + 
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
    geom_histogram(binwidth = 0.025) 
  
  print(p)
  p
}

#' Plot the per-sample FPKMs for a gene.
#'
#' @param gene_identifier A string denoting either the gene entrez_id, ensembl_id or gene_symbol.
#' @param result_table A string of the path to the "deseq2_results_fpkm.csv" file, which contains the gene 
#'   FPKM info, or a results table object generated by a diff_expr.R script. If left empty, will read 
#'   "./results/differential_expression/de_gene/deseq2_results_fpkm_{SPECIES}.csv" by default.              
#' @param debug A boolean, whether to print a markdown table of the the gene FPKMs.
#' @param print_graph A boolean, whether to plot the graph before returning the ggplot object.
#' @param feature_group A string vector. The result table should contains FPKM columns which have all the 
#'   features separated by '_', for example "101_WT_Hip_Ctrl_fpkm". This will be split into columns listed in 
#'   feature_group, which in this case can be c('sample_id', 'condition', 'region', 'treatment'),
#' @param filter_string A string, a filter which will be applied when selecting the samples from the results
#'   table. For example: filter="condition=='5xFAD'" or filter="str_detect(sample_meta, 'SC[0-9]')". 
#'   If left NULL, all samples will be used.
#' @param plot_feature A string vector, indicating the usage of aes on the features in feature_group. The 
#'   length of the vector should be the same as the feature_group.
#' @param plot_label A string from feature_group. Which feature to be used as label for the points in the plot.
#' @param plot_x A string from feature_group. Which feature to be used as the x-axis.
#' @return A list containing two items: debug_inf and ggplot object with keys {'info', 'graph'}
#' @examples
#'plot_gene_fpkms('ENSMUSG00000029816', result_table=result_table,debug=FALSE,print_graph=FALSE,
#'              feature_group=c('sample_id','condition','region','treatment'), filter="condition=='5xFAD'",
#'              plot_feature=c('','','','color'),
#'              plot_label="sample_id",plot_x="region")
#'
#' genes<-c('ENSG00000171885','ENSG00000131095')
#' p=''
#' for(g in split(genes, ceiling(seq_along(genes)/4)) ){
#'   for( i in g ){
#'     plot_name<-str_c('p','_',i)
#'     assign(plot_name,plot_gene_fpkms(i, result_table=results, debug=FALSE, print_graph=FALSE))
#'     if(nchar(p)==0){
#'       p<-plot_name
#'     }else{
#'       p<-str_c(p,plot_name,sep = '+')
#'     }
#'   }
#'   if(length(g)<2){
#'     p %>% parse_expr() %>% eval() %>% print()
#'   }else{
#'     str_c(p," plot_layout(ncol = 2)",sep = '+') %>% parse_expr() %>% eval() %>% print()
#'   }
#'   p=''
#' }
plot_gene_fpkms <- function(gene_identifier, result_table = NULL, debug = FALSE, print_graph = FALSE,
                            feature_group = c(), filter_string = '', plot_feature = c(), 
                            plot_label = "", plot_x = "") {
  
  # devtools::install_github("thomasp85/patchwork",force = TRUE)
  # devtools::install_github("slowkow/ggrepel")
  require(patchwork)
  require(ggrepel)
  require(knitr)
  
  if (is.null(result_table)) {
    result_table <- str_c('./results/differential_expression/de_gene/deseq2_results_fpkm_', SPECIES, '.csv')
  }
  
  if (is_string(result_table)) {
    if (!file.exists(result_table)) {
      stop(str_c("result_table [", result_table, "] does not exist."))
    }
    
    result_table <- read.csv(result_table)
  }

  # get the results table, and get the fpkms, and also filter by gene name
  # result should be the fpkms for our selected gene we want to plot
  fpkm_debug <- result_table %>% 
    dplyr::select(gene, gene_name, chromosome, entrez_id,
                  dplyr::ends_with("_fpkm"), 
                  -dplyr::ends_with("avg_fpkm"),
                  -dplyr::ends_with(".stat")) %>%
    dplyr::filter(gene_name==gene_identifier | gene==gene_identifier | entrez_id==gene_identifier)

  # get the gene name from the results table
  gene_name <- fpkm_debug$gene_name %>% as.vector()

  # turn the tibble into long format: each row is a sample and its fpkm etc
  fpkm_debug_long <- fpkm_debug %>% 
    as_tibble() %>%
    melt(id.var=c("gene", "gene_name", "chromosome", "entrez_id"),
         variable.name = 'sample_meta', value.name = 'fpkm')

  # check if featuregroup has been set; if not enter this code
  if (feature_group %>% length() == 0) {
    ## No feature group provided, we are going to plot the FPKM using the sample name and color.
    plot_x <- 'sample_meta'
    plot_label <- 'sample_meta'
    
    if (filter_string!='') {
      fpkm_debug_long %<>% filter(!!parse_expr(filter_string))
    }

    # we force the sample order on x axis
    fpkm_debug_long$sample_meta %<>%factor(levels=unique(fpkm_debug_long$sample_meta))

    p <- fpkm_debug_long %>%
      ggplot(mapping=aes_string(y="fpkm", x = plot_x)) + 
      aes_string(color = plot_x) + 
      geom_point(size = 3) +
      geom_text_repel(aes(label = !!parse_expr(plot_label)),
                      nudge_x = -0.35, direction = "y", hjust = 0.5, 
                      segment.size = 0.1, size = 3) +
      ggtitle(gene_identifier %>% str_c(gene_name, sep = ':')) + 
      theme(legend.position = 'none') +
      theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

  # if featuregroup has been set, instead enter this code
  } else {
    # tell it that we are using sample meta column to label the points on the plot later
    plot_label <- 'sample_meta'

    # remove the trailing fpkm from end of smaple names, resulting in sample names
    # trims only trailing _fpkm in case theres an fpkm in the sample name
    fpkm_debug_long %<>% mutate(sample_meta = sub("_fpkm$", "", sample_meta))

    # join sample data with fpkm table by sample name
    # all sample data variables are now included
    fpkm_debug_long %<>% inner_join(SAMPLE_DATA, by=c("sample_meta" = "sample_name"))

    # We want to plot only these samples
    if (filter_string != '') {
      fpkm_debug_long %<>% filter(!!parse_expr(filter_string))
    }

    # we force the sample order on x axis
    fpkm_debug_long$sample_meta %<>%factor(levels=unique(fpkm_debug_long$sample_meta))

    # set up plot, giving x and y variables
    p <- fpkm_debug_long %>% ggplot(mapping = aes_string(y = "fpkm", x = "sample_meta"))

    # where set, change features of plot to be by our selected feature group
    for (i in which(plot_feature!='')) {
      switch(plot_feature[i],
             'color' = p <- p + aes_string(color=feature_group[i]),
             'shape' = p <- p + aes_string(shape=feature_group[i])
      )
    }

    p <- p + geom_point(size = 3) +
      geom_text_repel(aes(label = !!parse_expr(plot_label)),
                      nudge_x = -0.35, direction = "y", hjust= 0.5, 
                      segment.size = 0.1, size = 3) +
      ggtitle(gene_identifier %>% str_c(gene_name, sep = ':')) + 
      theme(legend.position = "top") +
      theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
  }
  
  if (print_graph) {
    p %>% print()
  }
  
  if(debug) {
    fpkm_debug %>% kable(format = 'markdown', digits = 99) %>% print()
  }
  
  list("info" = fpkm_debug, 'graph' = p)
}

#' This function calls the plot_gene_fpkms function and plots FPKMs for the genes defined in the 
#' MARKER_GENES table.
#' @param result_table Passed directly to the plot_gene_fpkms function; this can either be:
#'   - a string, the path to the "deseq2_results_fpkm.csv", which contains the gene FPKM info, or
#'   - a result table object generated by diff_expr.R script.
#'  If left empty, will by default read: 
#'    "./results/differential_expression/de_gene/deseq2_results_fpkm_{SPECIES}.csv"
check_cell_type <- function(result_table, fpkm_check_cutoff = 5,
                            print_check_log = TRUE, print_fpkm_table = FALSE) {
  
  gene_markers <- GENE_MARKERS %>% pull(SPECIES) %>% split(f = GENE_MARKERS$cell_type)
  fpkm_info <- data.frame()
  message("Cell type check results are saved in [", GRAPHS_DIR,"]")
  
  for (cell_type in gene_markers %>% names()) {
    genes <- gene_markers %>% extract2(cell_type)
    num_genes <- length(genes)
    
    plot_statement <- ''
    
    for (index in seq(num_genes)) {
      plot_name <- str_c(cell_type, '_', genes[index])
      l <- plot_gene_fpkms(genes[index], result_table = result_table, debug = FALSE, print_graph = FALSE, feature_group = CELLTYPE_FEATURE_GROUP, plot_feature = CELLTYPE_PLOT_FEATURE)
      l$graph <- l$graph + 
        geom_hline(yintercept = fpkm_check_cutoff, linetype = "dashed", color = "grey") + 
        geom_hline(yintercept = 0, color = "grey") + 
        ylim(0, NA) + 
        theme_classic() + 
        theme(legend.position = "right",
              axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank())
      
      assign(plot_name, l$graph)
      fpkm_info <- l$info %>% mutate(gene_marker_cell_type = cell_type) %>% rbind(fpkm_info)
      
      plot_statement <- ifelse(index == 1, plot_name, plot_statement %>% str_c(plot_name, sep = '+'))
    }
    
    start_plot(str_c("cell_type_check_", cell_type))
    
    if (num_genes > 1) {
      plot_statement %<>% str_c(" + plot_layout(ncol = ", num_genes, ", guides = \"collect\")") 
    }
    plot_statement %>% parse_expr() %>% eval() %>% print()
    
    end_plot()
  }
  
  ## check fpkms
  if (print_check_log) {
    # gene_marker_cell_type       sample         fpkm    is
    # 1        oligodendrocyte  IN_AS1_fpkm 2.361160e+00 FALSE
    # 2        oligodendrocyte  IN_AS1_fpkm 5.220813e-03 FALSE
    # 3                 neuron  IN_AS1_fpkm 1.467844e-02 FALSE
    # 4                 neuron  IN_AS1_fpkm 1.203420e+00 FALSE
    # 5              microglia  IN_AS1_fpkm 0.000000e+00 FALSE
    # 6              microglia  IN_AS1_fpkm 1.071327e-02 FALSE
    # 7              astrocyte  IN_AS1_fpkm 2.133310e+03  TRUE
    # 8              astrocyte  IN_AS1_fpkm 1.320389e+02  TRUE
    # 9        oligodendrocyte  IN_AS2_fpkm 1.769229e+00 FALSE
    check_result <- fpkm_info %>%
      dplyr::select(contains('fpkm'), gene_marker_cell_type) %>%
      reshape2::melt(id.vars = c('gene_marker_cell_type'), variable.name = 'sample', value.name = 'fpkm') %>%
      mutate(is = fpkm > fpkm_check_cutoff)
    
    # sample       is
    # <fct>        <chr>
    #   1 IN_AS1_fpkm  astrocyte
    # 2 IN_AS2_fpkm  astrocyte
    # 3 IN_AS3_fpkm  astrocyte
    # 4 IN_AG31_fpkm astrocyte
    # 5 IN_AG32_fpkm astrocyte
    # 6 IN_AG33_fpkm astrocyte
    # 7 IN_AE61_fpkm astrocyte
    # 8 IN_AE62_fpkm astrocyte
    # 9 IN_AE63_fpkm astrocyte
    ## For each sample, amount all the cell types, which cell type has the most gene markers passed the cutoff?
    check_result %<>% 
      group_by(sample, gene_marker_cell_type) %>%
      summarise(like = sum(is)/n()) %>%
      # summarise(is=gene_marker_cell_type[which(like == max(like) )] %>% paste(collapse = ' / '))
      summarise(is = gene_marker_cell_type[which(like >=0.5)] %>% paste(collapse = ' / '))
    
    cat("Cell type check result:\n")
    for (i in check_result %>% pull(sample) %>% levels()) {
      str_c(i %>% strsplit('_fpkm') %>% extract2(1),
            ' looks like **',
            check_result %>% filter(sample == i) %>% pull(is) %>% toupper(),
            "**\n") %>% cat()
    }
  }
  
  if (print_fpkm_table) {
    fpkm_info %>% print()
  }
}

# read picard qc matrix for rna degradation
read_median_5prime_to_3prime_bias <- function(samples=c()){
  PICARD_ALIGNMENT_MATRIX=file.path('results/','alignment_metrics',SPECIES)
  MEDIAN_5PRIME_TO_3PRIME_BIAS <- samples %>% lapply(function(sample){
    read_tsv(file.path(PICARD_ALIGNMENT_MATRIX,str_c(sample,'.txt')),comment = '#',n_max = 1,col_type = cols()) %>%
      pull(MEDIAN_5PRIME_TO_3PRIME_BIAS)
  }) %>% unlist() %>% setNames(samples)
}

#### Transcript-level D. E. analyses ####

get_transcripts_to_genes <- function(species = 'human') {
  str_c("data/", species , "_ensembl_{{cookiecutter.ensembl_version}}/tx2gene.tsv") %>%
    read_delim(delim = " ", col_names = c("transcript", "gene"))
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
  
  transcript_tpms %<>% inner_join(get_transcripts_to_genes(SPECIES), by = c("target_id"="transcript"))
  
  transcript_tpms %>% 
    group_by(gene) %>% 
    summarise(tpm = sum(tpm)) %>%
    rename_(.dots = setNames(names(.), c("gene", str_c(sample, "_tpm"))))
}

get_tximport <- function(sample_data, quant_method='salmon', tx_level = TRUE) {
  quant_file <- get_transcript_quant_file(quant_method)
  quant_dirs <- sample_data %>%
    tibble::rownames_to_column(var = "tmp") %>%
    pull("tmp")
  
  quant_files <- str_c("results/",quant_method,"_quant/", SPECIES, "/",  quant_dirs, "/", quant_file)
  names(quant_files) <- quant_dirs
  
  txi <- tximport(quant_files, type=quant_method, txOut = tx_level, 
                  tx2gene=get_transcripts_to_genes(SPECIES), dropInfReps = TRUE)

    txi$Length <- read.csv(quant_files[1],sep = '\t',stringsAsFactors = FALSE) %>%
    dplyr::select(id = 1, length = 2) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(var = 'id')
  
  txi
}

get_total_dds_tximport <- function(sample_data, quant_method = 'salmon', tx_level = TRUE) {
  txi <- get_tximport(sample_data, quant_method, tx_level)
  
  total_dds <- DESeqDataSetFromTximport(txi, sample_data, ~1)
  total_dds <- DESeq(total_dds)
  
  total_dds
}


get_avg_tpm <- function(tpms, tx_level) {
  sample_data <- SAMPLE_DATA
  sample_data %<>% tibble::rownames_to_column(var = "tmp_row_names") %>% 
    group_by(.dots = AVG_FPKM_GROUP) %>%
    summarise(samples = str_c(tmp_row_names, '_tpm', sep = '', collapse = ',')) %>%
    tidyr::unite('avg_name', AVG_FPKM_GROUP, sep = '_')
  
  for (avg in sample_data %>% pull(avg_name)) {
    samples <- sample_data %>%
      filter(avg_name == avg) %>%
      pull(samples) %>%
      str_split(',') %>% 
      unlist()
    
    avg_tpm <- tpms %>% dplyr::select(one_of(samples)) %>%
      mutate(avg_tpm = rowMeans(.)) %>%
      dplyr::pull(avg_tpm)
    
    tpms %<>% mutate(!!str_c(avg, '_avg_tpm', sep='') := avg_tpm)
  }
  
  id_column <- ifelse(tx_level, 'transcript', 'gene')
  
  tpms %>% dplyr::select(!!id_column, contains('avg'))
}

#### Gene ontology enrichment analyses ####

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

  go_results <- go_data %>% GenTable(weight_fisher = result_weight, orderBy = "weight_fisher", topNodes = 150)

  gene_info <- get_gene_info(species)
  go_results$Genes <- sapply(go_results[,c('GO.ID')], function(x) get_significant_genes(x, go_data, gene_info))

  list(go_results = go_results,
       # go_data = go_data,
       # result_classic = result_classic,
       # result_elim = result_elim,
       result_weight = result_weight
  )
}

perform_go_analyses <- function(significant_genes, expressed_genes, comparison_name, file_prefix, species, out_dir="results/differential_expression/go/") {
  if (significant_genes %>% nrow == 0) {
    message("No significant genes supplied.")
    return()
  }

  top_dir <-  file.path(out_dir,species,comparison_name)
  if (!dir.exists(top_dir)) {
    dir.create(top_dir, recursive = TRUE)
  }

  sapply(c("BP", "MF", "CC"), simplify = FALSE, USE.NAMES = TRUE, function(x) {
    ret <- perform_go_analysis(expressed_genes, significant_genes, x, species, top_dir ,comparison_name)
    
    ret %>% 
      extract2('go_results') %>% 
      inner_join(GO_TERMS) %>% 
      dplyr::mutate(Term=FullTerm) %>% 
      dplyr::select(-FullTerm) %>%
      dplyr::rename(annotated_in_background = Annotated, annotated_in_gene_set = Significant,
                    expected_annotated_in_gene_set = Expected, p.value = weight_fisher) %>% # Changing column names in results
      write_csv(file.path(top_dir, str_c(comparison_name, file_prefix, "_go_", x %>% tolower, ".csv")), na="")

    ret
  })
}

#### Reactome pathway analysis ####

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
        write_csv(file.path(top_dir, str_c(comparison_name, file_prefix, "_reactome.csv")), na="")

      ret <- pathways %>% dplyr::select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, geneID)
    }
    
    ret
}

#### Gene set enrichment analysis ####

get_human_vs_species_ortholog_info <- function(species) {
  orthologs <- species %>% 
    str_c("data/", ., "_ensembl_{{cookiecutter.ensembl_version}}/human_orthologs.tsv") %>% 
    read_tsv(col_names = c("species_gene", "human_gene", "type"))
  
  species_genes <- get_gene_info(species)
  human_genes <- get_gene_info("human")
  
  orthologs_entrez <- orthologs %>% 
    left_join(species_genes %>% 
                dplyr::select(gene, entrez_id) %>% 
                dplyr::rename(species_gene = gene, species_entrez_id = entrez_id)) %>%
    left_join(human_genes %>% 
                dplyr::select(gene, entrez_id) %>% 
                dplyr::rename(human_gene = gene, human_entrez_id = entrez_id)) %>%
    dplyr::select(species_entrez_id, human_entrez_id) %>% 
    filter(!is.na(species_entrez_id) & !is.na(human_entrez_id)) %>% 
    distinct()
  
  same_name_entrez <- (species_genes %>% 
                         mutate(gene_name = tolower(gene_name)) %>% 
                         dplyr::select(gene_name, entrez_id) %>% 
                         filter(!is.na(entrez_id)) %>% 
                         dplyr::rename(species_entrez_id = entrez_id)) %>%
    inner_join(human_genes %>% 
                 mutate(gene_name = tolower(gene_name)) %>% 
                 dplyr::select(gene_name, entrez_id) %>% 
                 filter(!is.na(entrez_id)) %>% 
                 dplyr::rename(human_entrez_id = entrez_id)) %>%
    dplyr::select(-gene_name) %>% 
    distinct()
  
  human_species_entrez_mappings <- orthologs_entrez %>% 
    rbind(same_name_entrez) %>% 
    distinct
}


get_gene_sets <- function(species, gene_set_name) {
    str_c("data/",species,"/msigdb/v7.0/", gene_set_name, ".all.v7.0.entrez.gmt.Rdata") %>% load
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

  camera_results %>% tibble::rownames_to_column(var = "GeneSet") %>% filter(FDR < fdr_cutoff) %>%
    write_csv(file.path(top_dir, str_c(comparison_name, "-", gene_set_collection_name, "_sets.csv")), na="")

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
    write_csv(file.path(top_dir, str_c(comparison_name, "-", gene_set_collection_name, "_genes_in_sets.csv")), na="")

  ret[['genes_in_sets']]<- de_results %>% cbind(gene_set_results)

  ret
}

# function to plot heatmap of genes in significant gene sets, per comparison
plot_significant_set_heatmap <- function(set_name, all_sets, comparison, samples_in_comparison, species) {
  # get names of FPKM data columns
  samples_in_comparison %<>% mutate(fpkm_columns = str_c(sample_name, "_fpkm"))
  fpkm_columns <- samples_in_comparison %>% pull(fpkm_columns)
  num_samples <- length(fpkm_columns)

  # construct sample annotation
  rownames(samples_in_comparison) <- NULL
  samples_in_comparison %<>% tibble::column_to_rownames(var = "fpkm_columns")
  annot <- samples_in_comparison %<>% dplyr::select(-sample_name)

  # loop over each gene set category
  for (i in names(all_sets[[set_name]])){
    # get the significant entrez ids
    entrez_ids <- list_of_gene_sets[[set_name]][[i]]

    # get the name of the column of the log 2 fold changes
    comparison_log2fc <- paste(comp, "l2fc", sep = '.')

    # from the global results variable, get rows and columns corresponding to significant entrez IDs
    # and samples from the correct comparison; remove genes with no name
    to_heatmap <- results %>%
      dplyr::select(gene_name, entrez_id, fpkm_columns, comparison_log2fc) %>%
      filter(entrez_id %in% entrez_ids) %>%
      filter(!is.na(gene_name))

    # reorder
    to_heatmap <- to_heatmap %>% arrange(desc(!!sym(comparison_log2fc)))

    # declare the path to the heatmaps, based on the gene set category and comparison
    # create the dir if it doesnt exist
    heatmap_path <- paste("results/differential_expression/gene_set_tests/", species, comparison, set_name, "/", sep = "/")

    ifelse(!dir.exists(file.path(heatmap_path)), dir.create(file.path(heatmap_path), recursive = TRUE), FALSE)

    # TODO: currently, the column to rownames call complains about duplicate row names (i.e. gene names)
    # I have removed duplicates - is this how we want to do this? Is there a better way?
    to_heatmap_unique <- distinct(to_heatmap, gene_name, .keep_all = TRUE)

    # prepare heatmap data so the row names are the gene IDs and remove the entrez column
    heatmap_data<-to_heatmap_unique %>% tibble::column_to_rownames(var="gene_name") %>% dplyr::select(-entrez_id, -comparison_log2fc)

    # divide each row by the mean of that row
    heatmap_data <- t(apply(heatmap_data, 1, function(x) x/mean(x)))
    heatmap_data <- log2(heatmap_data)

    heatmap_data[is.infinite(heatmap_data)] <- NA
    heatmap_data[is.nan(heatmap_data)] <- NA
    heatmap_data <- subset(heatmap_data,rowSums(!is.na(heatmap_data)) == num_samples)

    if(nrow(heatmap_data)!=0){
      max_data <- max(heatmap_data, na.rm=TRUE)
      min_data <- -min(heatmap_data, na.rm=TRUE)
      range <- min(max_data, min_data)

      start_plot(prefix = i, path = heatmap_path)
      pheatmap(heatmap_data,
      breaks = seq(-range, range, length.out = 100),
      cluster_rows = FALSE, cluster_cols = FALSE,
      border_color = NA, show_rownames = (heatmap_data %>% nrow()) < 100,
      annotation_col = annot)
      end_plot()
    }
  }
}

#' This function reads in all the GSA results and tracks the given gene sets
#' @example:
#' track_gene_sets(target_terms=c('GO_RIBOSOME','GO_PROTEASOME_COMPLEX','GO_TRANSLATIONAL_INITIATION'), category='GO',comparison_table=COMPARISON_TABLE)
track_gene_sets <- function(target_terms = c('GO_RIBOSOME', 'GO_PROTEASOME_COMPLEX', 'GO_TRANSLATIONAL_INITIATION'), 
                            category = 'GO',
                            gs_results = get_global('GS_results'),
                            comparison_table = COMPARISON_TABLE,
                            print_table = TRUE, 
                            output_table_file = NA,
                            left_out_comparison = c(),
                            heat_map.fdr.midpoint = 0.05){
  
  ## create a master table for all comparisons and the gene sets results
  gsa_res_tb <- COMPARISON_TABLE %>% 
    pull(comparison) %>% 
    extract(which(!. %in% left_out_comparison)) %>% 
    set_names(.) %>%
    sapply(simplify = FALSE, USE.NAMES = TRUE, function(x) {
      gs_results %>% 
        extract2(x) %>% 
        extract2(category) %>%
        tibble::rownames_to_column('GeneSet') %>% 
        dplyr::mutate(comparison = x)
    }) %>% 
    reduce(rbind)
  
  ## we only keep gene sets of interests
  gsa_res_tb <- gsa_res_tb %>% 
    filter(GeneSet %in% target_terms) %>%
    mutate(log10fdr=log10(FDR)) # for plotting
  
  gsa_res_tb$comparison %<>% factor(levels = COMPARISON_TABLE %>% pull(comparison) %>% rev)
  gsa_res_tb$GeneSet %<>% factor(levels = target_terms)
  
  if (print_table) {
    print(gsa_res_tb %>% arrange(GeneSet, FDR))
  }
  
  # we output the FDR table to csv
  if (!is.na(output_table_file)) {
    gsa_res_tb %>% 
      arrange(GeneSet, FDR) %>% 
      dplyr::select(GeneSet, FDR, comparison) %>%
    reshape(idvar = "comparison", timevar = "GeneSet", direction = "wide") %>%
    write.csv(file = output_table_file)
  }
  
  gsa_res_tb %>% 
    ggplot(aes(GeneSet, comparison)) +
    geom_tile(aes(fill = log10fdr)) + 
    scale_fill_gradient2(low = "red", high = "white", mid = "white",
                         midpoint = log10(heat_map.fdr.midpoint)) +
    geom_text(aes(label = if_else(FDR < heat_map.fdr.midpoint, 
                                  ifelse(Direction=='Up', "UP", "DOWN"),
                                  "")), alpha = 0.75, size=3) +
    labs(fill ="log10(FDR)") +
    labs(title=str_c('FDR cutoff = ',heat_map.fdr.midpoint),
         xlab="Gene Set", ylab="Comparison") + 
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = -90, size=7), 
          panel.border = element_blank(),panel.background = element_blank())
}

track_go <- function(target_terms = c('GO:0051492', 'GO:0010811'),
                            category = 'BP',
                            gs_results = get_global('GO_results'),
                            comparison_table = COMPARISON_TABLE,
                            print_table = TRUE,
                            output_table_file = NA,
                            left_out_comparison = c(),
                            heat_map.p.midpoint = 0.05){

  ## create a master table for all comparisons and the gene sets results
  gsa_res_tb <- COMPARISON_TABLE %>%
    pull(comparison) %>%
    extract(which(!. %in% left_out_comparison)) %>%
    set_names(.) %>%
    sapply(simplify = FALSE, USE.NAMES = TRUE, function(x) {
      gs_results %>%
        extract2(x) %>%  extract2(str_c(x,'.all')) %>%
        extract2(category) %>%  extract2('go_results') %>%
        dplyr::mutate(comparison = x)
    }) %>%
    reduce(rbind)

  ## we only keep gene sets of interests
  gsa_res_tb <- gsa_res_tb %>%
    filter(GO.ID %in% target_terms) %>%
    mutate(log10p=log10(as.numeric(weight_fisher))) # for plotting

  gsa_res_tb$comparison %<>% factor(levels = COMPARISON_TABLE %>% pull(comparison) %>% rev)
  gsa_res_tb$GO.ID %<>% factor(levels = target_terms)

  if (print_table) {
    print(gsa_res_tb %>% arrange(GO.ID, log10p))
  }

  # we output the FDR table to csv
  if (!is.na(output_table_file)) {
    gsa_res_tb %>%
      arrange(GO.ID, log10p) %>%
      dplyr::select(GO.ID, weight_fisher, comparison) %>%
      reshape(idvar = "comparison", timevar = "GO.ID", direction = "wide") %>%
      write.csv(file = output_table_file)
  }

  gsa_res_tb %>%
  ggplot(aes(GO.ID, comparison)) +
    geom_tile(aes(fill = log10p)) +
    scale_fill_gradient2(low = "red", high = "white", mid = "white",
    midpoint = log10(heat_map.p.midpoint)) +
  # geom_text(aes(label = if_else(log10p < log10(heat_map.p.midpoint),
  #                               ifelse(Direction=='Up', "UP", "DOWN"),
  #                               "")), alpha = 0.75, size=3) +
    labs(fill ="log10(p)") +
    labs(title=str_c('p cutoff = ',heat_map.p.midpoint),
    xlab="GO.ID", ylab="Comparison") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = -90, size=7),
    panel.border = element_blank(),panel.background = element_blank())
}

#### Quality surrogate variable analysis ####

qsva <- function(degradationMatrix, mod = matrix(1, ncol = 1, nrow = ncol(degradationMatrix))) {
  degPca = prcomp(t(log2(degradationMatrix + 1)))
  k = num.sv(log2(degradationMatrix + 1), mod, seed={{ range(1, 1000) | random }})
  degPca$x[, seq_len(k)]
}

get_quality_surrogate_variables <- function(dds) {
  design_formula <- dds %>% design()
  sample_data <- dds %>% colData()
  sample_names <- sample_data %>% rownames()
  #@todo adding of the interaction formula might affect this, needs to come back and check later
  quality_surrogate_variables <- sample_names %>%
    str_c("results/read_counts/", ., ".", SPECIES, ".dm.tsv") %>%
    read.degradation.matrix(
      sampleNames = sample_names,
      totalMapped = dds %>% counts %>% colSums,
      readLength = 75,
      type = "region_matrix_single",
      BPPARAM = SerialParam()) %>%
    qsva(mod = design_formula %>% model.matrix(sample_data))

  quality_surrogate_variables %>% print

  quality_surrogate_variables
}

get_qsva_dds <- function(dds) {
  design_formula <- dds %>% design()
  sample_data <- dds %>% colData()
  quality_surrogate_variables <- get_quality_surrogate_variables(dds)

  if (is.vector(quality_surrogate_variables)) {
    quality_surrogate_variables %<>% data.frame(qSVA = .)
  }

  colData(dds) %<>% cbind(quality_surrogate_variables)

  design(dds) <- design_formula %>%
    terms() %>% 
    attr("term.labels") %>%
    c(quality_surrogate_variables %>% colnames, .) %>% 
    reformulate
  
  dds
}

#### Sargasso misassignment percentages ####

get_gene_counts_and_lengths <- function(samples, species, gene_lengths, sum_counts = FALSE) {
  res <- samples %>%  
    map(read_counts, species) %>%
    purrr::reduce(inner_join)
  
  if (sum_counts) {
    # res %<>% mutate(counts = rowSums(.[,-1])) %>% dplyr::select(gene, counts)
    res %<>% mutate(counts = rowSums(.[,-1]))
  }
  
  res %>% inner_join(gene_lengths)
}

calculate_total_reads_for_species <- function(samples, species, sum_counts = FALSE) {
  res <- samples %>%  
    map_dbl(function(sample) {
      species %>% map(~read_counts(sample, .) %>%
                        dplyr::pull(2) %>%  
                        sum()) %>%
        purrr::reduce(sum) 
    })   
  
  if (sum_counts) {
    res %<>% set_names(samples)
    res['counts']=sum(res)
  } else {
    res %<>% set_names(samples)
  }
  res
}

fpkm_formula <- function(read_counts, gene_lengths, total_reads) {
  10^9 * as.double(read_counts) / as.double(gene_lengths) / as.double(total_reads)
}

calculate_fpkm <- function(samples, gene_counts_and_lengths, total_reads){
  samples %>%  
    map(function(sa) {
      read_counts <- gene_counts_and_lengths %>% pull(sa)
      gene_lengths <- gene_counts_and_lengths %>% pull(max_transcript_length)
      gene_counts_and_lengths %>%
        dplyr::select(gene) %>%
        mutate(!!sa:=10^9 * as.double(read_counts) /
                 as.double(gene_lengths) /
                 as.double(total_reads[[sa]]))
    }) %>% reduce(inner_join)
}

read_overall_filtering_summary <-function(){
    read_csv(file.path("results","sargasso","filtered_reads","overall_filtering_summary.txt"))
}

calculate_per_gene_misassignment_percentages <- function(species_of_interest,
                                                         target_samples,target_species,
                                                         reference_samples,reference_species,
                                                         paired=FALSE,
                                                         gene_lengths){

  # In a mix-species sample, for example samples contains human/mouse/rat(HMR) cells,
  # we want to eastimate, for a 'species_of_interest', the amounnt of reads due to misassignment.
  # This is done by using 'reference samepls' does not contains cells from 'species_of_interest'
  # species_of_interest='mouse'  # the species we want to eastimate the misassignment
  # target_samples="01_HMR,02_HMR,03_HMR,04_HMR"
  # target_species="human,mouse,rat"
  # reference_samples="09_Pc,10_Pc,11_Pc,12_Pc,05_EC,06_EC,07_EC,08_EC"
  # reference_species="rat;rat;rat;rat;human;human;human;human"

  # x <- COMPARISON_TABLE %>% filter(comparison == comparison_name)
  # y <- MISASSIGNMENT_SAMPLE_REFERENCE_TABLE %>%
  #   filter(comparison_name == x$comparison) %>%
  #   filter(condition == x$condition_base)
  # target_samples <- SAMPLE_DATA %>% filter(!!parse_expr(x$filter)) %>%
  #   filter(!!parse_expr(x$condition_name) == x$condition) %>%
  #   pull(sample_name) %>%
  #   str_c(collapse = ',')
  # target_samples=y$target_samples
  # species_of_interest=y$species_of_interest
  # target_species=y$target_species
  # reference_samples=y$reference_samples
  # reference_species=y$reference_species

  target_samples %<>% strsplit(.,',') %>% unlist()
  target_species %<>% strsplit(.,',') %>% unlist()
  reference_samples  %<>% strsplit(.,',') %>% unlist()
  reference_species  %<>% strsplit(.,';') %>% unlist()

  # we workout all species involved
  all_species=c(target_species,
                species_of_interest,
                strsplit(reference_species,',') %>% unlist()
  ) %>% unique

  # For each 'species_of_interest' gene, calculate an approximate gene “expression” in 'reference_samples'(rat)
  # containing only mouse and human material, measured in fragments per kilobase
  # per million mapped reads: FPKM_mm_hs(i)= 10^9*c(i) / (l(i)*R)
  # where c(i) is the number of reads assigned to rat gene in the mouse plus human
  # samples, l(i) is the length of the longest transcript of the gene, and R is the
  # total number of reads assigned to all genes of all species in the mouse plus
  # human samples.

  ## step 30 in paper
  # we group the reference sample by reference species
  # $human
  # [1] "05_EC" "06_EC" "07_EC" "08_EC"
  #
  # $mouse
  # [1] "13_As"  "14_As" "15_As_replacement" "16_As"
  refs <- split(reference_samples,reference_species)

  ## If the reference samples are paired with the target samples,
  ## we check if they have the same number of sample.
  if(paired){
    sample_all_paired = lapply(refs,function(x){
      length(x)==length(target_samples)}
    ) %>%
      unlist() %>% all

    if(!sample_all_paired)
      stop('Difference number of reference samples / target samples, cannot match pair.')

    ## We print out the pairing
    lapply(1:length(target_samples),function(i){
      c(setNames(c(target_samples[i]),str_c(target_species,collapse = ',')),sapply(refs,extract2,i)) %>% t() %>% as.data.frame()
    })%>%reduce(rbind) %>% print()
  }

  # for each group of reference samples, we work out the species_of_interest fpkm (misassignment)
  species_of_interest_in_reference_sample_fpkm <- refs %>% lapply(function(rsa){ # rsa: reference_sample
    # gene               `05_EC` `06_EC` `07_EC` `08_EC` counts gene_length max_transcript_length
    # <chr>                <dbl>   <dbl>   <dbl>   <dbl>  <dbl>       <dbl>                 <dbl>
    # 1 ENSRNOG00000049070     315     522     630     663   2130        1025                  1025
    # 2 ENSRNOG00000054139     189     197     210     294    890        2325                  2325
    # 3 ENSRNOG00000054846      75      91     143     114    423        1527                  1527
    ref_gene_counts_and_lengths <- get_gene_counts_and_lengths(
      rsa, species_of_interest, gene_lengths, sum_counts = TRUE)

    # 05_EC     06_EC     07_EC     08_EC    counts
    # 36191427  31824279  41255399  46717964 155989069
    ref_total_reads <- calculate_total_reads_for_species(
      rsa, all_species, sum_counts = TRUE)

    # This is using the total number of reads(counts) to workout the fpkm.
    # The per-sample fpkm is just for debugging purpose
    # gene               `05_EC` `06_EC` `07_EC` `08_EC`  counts
    # <chr>                <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
    # 1 ENSRNOG00000049070  8.49   16.0    14.9    13.8    13.3
    # 2 ENSRNOG00000055837  3.55    3.36    3.71    3.01    3.39
    # 3 ENSRNOG00000054139  2.25    2.66    2.19    2.71    2.45
    # 4 ENSRNOG00000054846  1.36    1.87    2.27    1.60    1.78
    # 5 ENSRNOG00000046657  1.69    1.69    1.52    0.984   1.43
    # 6 ENSRNOG00000011175  0.130   0.148   0.0916  0.222   0.151
    # 7 ENSRNOG00000030807  0.0285  0.162   0.125   0.199   0.132
    c('counts',rsa) %>%
      calculate_fpkm( ref_gene_counts_and_lengths, ref_total_reads) %>%
      dplyr::arrange(gene)
  })

  ## step 31 on paper
  # work out the species_of_interest fpkm in the target samples
  # gene               `01_HMR` `02_HMR` `03_HMR` `04_HMR` counts gene_length max_transcript_length
  # <chr>                 <dbl>    <dbl>    <dbl>    <dbl>  <dbl>       <dbl>                 <dbl>
  # 1 ENSRNOG00000000001       52       41       44       38    175        1416                  1416
  # 2 ENSRNOG00000000007        0        0        1        0      1        3866                  3381
  # 3 ENSRNOG00000000008      593      726      808     1093   3220        1747                  1747
  # 4 ENSRNOG00000000009        1        4        2        0      7        1361                  1361
  # 5 ENSRNOG00000000010        4        0        0        0      4        2444                  2444
  target_gene_counts_and_lengths <- get_gene_counts_and_lengths(
    target_samples, species_of_interest, gene_lengths, sum_counts = TRUE) %>%
    dplyr::arrange(gene)

  # 01_HMR    02_HMR    03_HMR    04_HMR
  # 149055415 122723195  77948579 115554113
  target_total_reads <- calculate_total_reads_for_species(
    target_samples,species_of_interest, sum_counts = FALSE)

  # gene               `01_HMR` `02_HMR` `03_HMR` `04_HMR`
  # <chr>                 <dbl>    <dbl>    <dbl>    <dbl>
  # 1 ENSRNOG00000000001   1.63      1.48   0.948      0.881
  # 2 ENSRNOG00000000007   0         0      0.00902    0
  # 3 ENSRNOG00000000008  15.0      21.2   14.1       20.5
  # 4 ENSRNOG00000000009   0.0325    0.150  0.0448     0
  # 5 ENSRNOG00000000010   0.0725    0      0          0
  # 6 ENSRNOG00000000012   0.178     0      0.215      0
  # 7 ENSRNOG00000000017  16.5      12.1   15.9       10.3
  species_of_interest_in_target_sample_fpkm <- target_samples %>%
    calculate_fpkm(target_gene_counts_and_lengths, target_total_reads) %>%
    arrange(gene)


  ## step 32 on paper
  ## For each target sample, calculate the approximate ratio of
  ## the reference_species ('mouse plus human’) to the species_of_interest (rat) RNA in the sample.
  sargasso_filtering_summary <- read_overall_filtering_summary()


  ## for each reference species, we workout a P(% of misassignment)
  P<-names(species_of_interest_in_reference_sample_fpkm) %>% set_names(.) %>% lapply(function(s){
    # refrence species, this will hopefully handle the mono/co reference samples
    rsps <- strsplit(s,',') %>% unlist()
    rsa <- refs[[s]] # reference samlple

    #FPKM_mm_hsi, depend on 'paired' parameter, avg OR the per-sample fpkm will be used
    # gene                 avg `05_EC` `06_EC` `07_EC` `08_EC`
    # <chr>              <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
    # 1 ENSRNOG00000000001     0       0       0       0       0
    # 2 ENSRNOG00000000007     0       0       0       0       0
    # 3 ENSRNOG00000000008     0       0       0       0       0
    f1 <- species_of_interest_in_reference_sample_fpkm[[s]] %>%
      dplyr::rename(avg=counts)

    # d
    # 04_HMR    02_HMR    03_HMR    01_HMR
    # 0.4278663 0.2358971 0.8291291 0.2187826
    d <- sargasso_filtering_summary %>%
      dplyr::select(str_c('Assigned-Reads-',species_of_interest),str_c('Assigned-Reads-',rsps)) %>%
      mutate(toTargetRatio = (rowSums(.)-.[[1]])/.[[1]]) %>%
      dplyr::select(toTargetRatio) %>%
      mutate(Sample=sargasso_filtering_summary$Sample) %>%
      dplyr::select(Sample, toTargetRatio) %>%
      filter(Sample %in% target_samples) %>%
      tibble::deframe() %>% extract(target_samples)


    # FPKM_rni
    # gene               `01_HMR` `02_HMR` `03_HMR` `04_HMR`
    # <chr>                 <dbl>    <dbl>    <dbl>    <dbl>
    # 1 ENSRNOG00000000001   1.63      1.48   0.948      0.881
    # 2 ENSRNOG00000000007   0         0      0.00902    0
    # 3 ENSRNOG00000000008  15.0      21.2   14.1       20.5
    # 4 ENSRNOG00000000009   0.0325    0.150  0.0448     0
    # 5 ENSRNOG00000000010   0.0725    0      0          0
    f2 <- species_of_interest_in_target_sample_fpkm

    if (!all(f1$gene == f2$gene)) {
      stop('gene in reference are not the same as gene in target.')
    }

    if(paired){
      ## for each refernce, we have 1 fpkm
      #                       05_EC 06_EC 07_EC 08_EC
      # ENSRNOG00000000001     0     0     0     0
      # ENSRNOG00000000007     0     0     0     0
      # ENSRNOG00000000008     0     0     0     0
      # ENSRNOG00000000009     0     0     0     0
      # ENSRNOG00000000010     0     0     0     0
      # ENSRNOG00000000012     0     0     0     0
      fpkm1 = f1 %>% dplyr::select(-avg) %>%
        tibble::column_to_rownames(var = 'gene') %>%
        as.matrix()
      numerator <- (( fpkm1 %>% as.matrix() %>% t()) * d )  %>% t()
    }else{
      ## we have one fpkm for all ref samples
      #                     avg
      # ENSRNOG00000000001   0
      # ENSRNOG00000000007   0
      # ENSRNOG00000000008   0
      # ENSRNOG00000000009   0
      # ENSRNOG00000000010   0
      # ENSRNOG00000000012   0
      fpkm1 = f1 %>% dplyr::select(gene,avg) %>%
        tibble::column_to_rownames(var = 'gene')
      numerator <- (fpkm1 %>% as.matrix()) %*% (d %>% as.matrix() %>% t())
    }


    #                         01_HMR     02_HMR       03_HMR     04_HMR
    # ENSRNOG00000000001  1.62657851  1.4761767  0.947681184  0.8813083
    # ENSRNOG00000000007  0.00000000  0.0000000  0.009020439  0.0000000
    # ENSRNOG00000000008 15.03476793 21.1866088 14.105591117 20.5463534
    # ENSRNOG00000000009  0.03254444  0.1498372  0.044817198  0.0000000
    # ENSRNOG00000000010  0.07249261  0.0000000  0.000000000  0.0000000
    # ENSRNOG00000000012  0.17842088  0.0000000  0.214991667  0.0000000
    denominator <- f2 %>%
      tibble::column_to_rownames(var = 'gene') %>%
      dplyr::select(names(d)) %>%
      as.matrix()

    # gene               `01_HMR` `02_HMR` `03_HMR` `04_HMR`     p
    # <chr>                 <dbl>    <dbl>    <dbl>    <dbl> <dbl>
    # 1 ENSRNOG00000000001        0        0        0        0     0
    # 2 ENSRNOG00000000007      NaN      NaN        0      NaN     0
    # 3 ENSRNOG00000000008        0        0        0        0     0
    # 4 ENSRNOG00000000009        0        0        0      NaN     0
    # 5 ENSRNOG00000000010        0      NaN      NaN      NaN     0
    per_sample_p <- (100 * numerator /  denominator) %>%
      as.data.frame() %>% setNames(target_samples) %>%
      tibble::rownames_to_column('gene') %>%
      mutate_at(.vars = target_samples, list(~replace(., is.infinite(.), NaN))) %>%
      mutate(p = rowMeans(.[,target_samples], na.rm=TRUE)) %>%
      as_tibble()

    list(p=per_sample_p,d=d)
  })

  # $human
  # $human$p
  # # A tibble: 32,883 x 6
  # gene               `01_HMR` `02_HMR` `03_HMR` `04_HMR`     p
  # <chr>                 <dbl>    <dbl>    <dbl>    <dbl> <dbl>
  #   1 ENSRNOG00000000001        0        0        0        0     0
  # 2 ENSRNOG00000000007      NaN      NaN        0      NaN     0
  # 3 ENSRNOG00000000008        0        0        0        0     0
  # 4 ENSRNOG00000000009        0        0        0      NaN     0
  # 5 ENSRNOG00000000010        0      NaN      NaN      NaN     0
  # # … with 32,873 more rows
  #
  # $human$d
  # 01_HMR     02_HMR     03_HMR     04_HMR
  # 0.45306509 0.45660661 0.06097529 0.13513955
  #
  # $mouse
  # $mouse$p
  # # A tibble: 32,883 x 6
  # gene               `01_HMR` `02_HMR` `03_HMR` `04_HMR`     p
  # <chr>                 <dbl>    <dbl>    <dbl>    <dbl> <dbl>
  #   1 ENSRNOG00000000001        0        0        0        0     0
  # 2 ENSRNOG00000000007      NaN      NaN        0      NaN     0
  # 3 ENSRNOG00000000008        0        0        0        0     0
  # 4 ENSRNOG00000000009        0        0        0      NaN     0
  # 5 ENSRNOG00000000010        0      NaN      NaN      NaN     0
  # # … with 32,873 more rows
  #
  # $mouse$d
  # 01_HMR   02_HMR   03_HMR   04_HMR
  # 4.570748 4.239137 1.206085 2.337179

  # gene               p.human p.mouse     p
  # <chr>                <dbl>   <dbl> <dbl>
  # 1 ENSRNOG00000000001       0       0     0
  # 2 ENSRNOG00000000007       0       0     0
  # 3 ENSRNOG00000000008       0       0     0
  # 4 ENSRNOG00000000009       0       0     0
  # 5 ENSRNOG00000000010       0       0     0
  P_df<-P %>% lapply(extract2,'p') %>%
    purrr::reduce(left_join,by='gene',suffix=c(str_c('.',names(.)))) %>%
    dplyr::select(gene,starts_with('p')) %>%
    mutate(p=rowSums(.[2:ncol(.)]))

  # gene               human_counts `05_EC` `06_EC` `07_EC` `08_EC` mouse_counts `13_As` `14_As` `15_As_replacement` `16_As` `01_HMR` `02_HMR` `03_HMR` `04_HMR`
  # <chr>                     <dbl>   <dbl>   <dbl>   <dbl>   <dbl>        <dbl>   <dbl>   <dbl>               <dbl>   <dbl>    <dbl>    <dbl>    <dbl>    <dbl>
  # 1 ENSRNOG00000000001            0       0       0       0       0            0       0       0                   0       0   1.63      1.48   0.948      0.881
  # 2 ENSRNOG00000000007            0       0       0       0       0            0       0       0                   0       0   0         0      0.00902    0
  # 3 ENSRNOG00000000008            0       0       0       0       0            0       0       0                   0       0  15.0      21.2   14.1       20.5
  debug_df<-species_of_interest_in_reference_sample_fpkm %>% map2(names(.),.,function(s,df){
    df %>% dplyr::rename_at(vars(counts),funs(str_c(s,'_',.)))
  }) %>% purrr::reduce(left_join) %>% left_join(species_of_interest_in_target_sample_fpkm)

  # $human
  # 04_HMR     02_HMR     03_HMR     01_HMR
  # 0.13513955 0.45660661 0.06097529 0.45306509
  #
  # $mouse
  # 04_HMR   02_HMR   03_HMR   01_HMR
  # 2.337179 4.239137 1.206085 4.570748
  D_df <- P %>% lapply(extract2,'d')

  ## we work out the composition of the reference sample
  # refs_composition <- lapply(refs,function(ss){
  refs_composition <- sargasso_filtering_summary %>%
    # dplyr::filter(Sample %in% ss) %>%
    dplyr::select(Sample,starts_with('Assigned-Reads')) %>%
    tidyr::pivot_longer(-Sample,names_to='type') %>%
    group_by(Sample) %>%
    mutate(prec=round(100*value/sum(value),digits=1)) %>%
    tidyr::pivot_wider(Sample,names_from = type,values_from = prec ) %>%
    dplyr::rename_at(vars(-Sample),gsub,pattern='Assigned-Reads-',replacement='',fixed=TRUE) %>%
    ungroup()
  # })


  list(
    P=P_df %>% left_join(debug_df),
    D=D_df,
    refs_composition=refs_composition
  )
}

get_misassignment_percentages <- function(comparison_name, gene_lengths) {

  x <- COMPARISON_TABLE %>% filter(comparison == comparison_name)

  ret <- list()

  #for condition
  y <- MISASSIGNMENT_SAMPLE_REFERENCE_TABLE %>%
    filter(comparison_name == x$comparison) %>%
    filter(condition == x$condition)

  if (nrow(y) > 0) {
    # target_samples <- SAMPLE_DATA %>% filter(!!parse_expr(x$filter)) %>%
    #   filter(!!parse_expr(x$condition_name) == x$condition) %>%
    #   pull(sample_name) %>%
    #   str_c(collapse = ',')

    P_condition <- calculate_per_gene_misassignment_percentages(y$species_of_interest,
                                                                y$target_samples, y$target_species,
                                                                y$reference_samples, y$reference_species,
                                                                y$paired,
                                                                gene_lengths)
    ret[['P_condition']] = P_condition$P %>% dplyr::select(gene,p)
    ret[['P_condition_debug']]=P_condition
    ret[['condition_reference_samples']] = y$reference_samples

  } else {
    ret[['P_condition']] = NA
    ret[['P_condition_debug']] = NA
    ret[['condition_reference_samples']] = NA
  }



  #for condition_base
  y <- MISASSIGNMENT_SAMPLE_REFERENCE_TABLE %>%
    filter(comparison_name == x$comparison) %>%
    filter(condition == x$condition_base)

  if (nrow(y) > 0) {
    # target_samples <- SAMPLE_DATA %>% filter(!!parse_expr(x$filter)) %>%
    #   filter(!!parse_expr(x$condition_name) == x$condition_base) %>%
    #   pull(sample_name) %>%
    #   str_c(collapse = ',')
    P_condition_base <- calculate_per_gene_misassignment_percentages(y$species_of_interest,
                                                                     y$target_samples, y$target_species,
                                                                     y$reference_samples, y$reference_species,
                                                                     y$paired,
                                                                     gene_lengths)
    ret[['P_condition_base']] = P_condition_base$P %>% dplyr::select(gene,p)
    ret[['P_condition_base_debug']] = P_condition_base
    ret[['condition_base_reference_samples']] = y$reference_samples
  } else {
    ret[['P_condition_base']] = NA
    ret[['P_condition_base_debug']] =NA
    ret[['condition_base_reference_samples']] = NA
  }
  ret
}


check_sample_bam <- function(samples,species=SPECIES,chr='2',start=75505857,end=75510000){
  library('Rsamtools')
  library(GenomicAlignments)
  library(parallel)
  library(ggplot2)
  library(dplyr)
  bamRanges=GRanges(chr, IRanges(start,end))

  bamFiles <- list.files('results/final_bams/',pattern = '.bam$',full.names = T) %>%
  grep(pattern = paste(samples,collapse = '|'),value = T)
  bamIndexFile <- list.files('results/final_bams/',pattern = '.bam.bai',full.names = T) %>%
  grep(pattern = paste(samples,collapse = '|'),value = T)
  bamExperiment <-list(description="",created=date())
  bv <- BamViews(bamFiles, bamIndicies=bamIndexFile,bamRanges=bamRanges, bamExperiment=bamExperiment)
  reads<-readGAlignments(bv)

  olap1 <- endoapply(reads, subsetByOverlaps, bamRanges)
  olap1 <- lapply(olap1, "seqlevels<-", value=as.character(seqnames(bamRanges)))
  cvg <- endoapply(olap1, coverage,
  shift=-start(ranges(bamRanges[1])),
  width=width(ranges(bamRanges[1])))

  coverage_tb <- lapply(names(cvg), function(sample){
    a=cvg[[sample]][[1]]
    sapply(c(1:length(runLength(a))),function(x){
      rep(runValue(a)[x],runLength(a)[x])
    }) %>% unlist() %>% set_names(c(start(ranges(bamRanges))[1]:end(ranges(bamRanges))[1])) %>%
      as.data.frame() %>% setNames(c('count')) %>%
      tibble::rownames_to_column('position') %>%
      mutate(sample=sample,position=as.numeric(position))
  }) %>% purrr::reduce(rbind)
  ggplot2::ggplot(data=coverage_tb, aes(x=position, y=count)) +
  # geom_line(linetype="dashed", size=0.2) +
    geom_area() +
  # geom_bar(stat="identity")
  # geom_point() +
    facet_wrap(~sample)
}


# This function plots the read distribution of a chromosome region of bam files
# It can be used as a pre-'igv' check
# one can test the code with the following:
# >setwd('/srv/data/results/nrf2_ich_jamie_loan/aaaaaaaaarggggghhhhh-959851d2469757607aa3e1f8b0ec1e1cc533cf42/20210805')
# >check_sample_bam(samples=c('14_KO_Ma_ICH','14_KO_Mg_ICH'),species='mouse',chr='2',start=75505857,end=75505957) %>% plot
check_sample_bam <- function(samples=c('14_KO_Ma_ICH','14_KO_Mg_ICH'),species='mouse',chr='2',start=75505857,end=75510000){
  library('Rsamtools')
  library(GenomicAlignments)
  library(parallel)
  library(ggplot2)
  library(dplyr)
  bamRanges=GRanges(chr, IRanges(start,end))

  bamFiles <- list.files('results/final_bams/',pattern = '.bam$',full.names = T) %>%
  grep(pattern = paste(samples,collapse = '|'),value = T)
  bamIndexFile <- list.files('results/final_bams/',pattern = '.bam.bai',full.names = T) %>%
  grep(pattern = paste(samples,collapse = '|'),value = T)
  bamExperiment <-list(description="",created=date())
  bv <- BamViews(bamFiles, bamIndicies=bamIndexFile,bamRanges=bamRanges, bamExperiment=bamExperiment)
  reads<-readGAlignments(bv)

  olap1 <- endoapply(reads, subsetByOverlaps, bamRanges)
  olap1 <- lapply(olap1, "seqlevels<-", value=as.character(seqnames(bamRanges)))
  cvg <- endoapply(olap1, coverage,
  shift=-start(ranges(bamRanges[1])),
  width=width(ranges(bamRanges[1])))

  coverage_tb <- lapply(names(cvg), function(sample){
    a=cvg[[sample]][[1]]
    sapply(c(1:length(runLength(a))),function(x){
      rep(runValue(a)[x],runLength(a)[x])
    }) %>% unlist() %>% set_names(c(start(ranges(bamRanges))[1]:end(ranges(bamRanges))[1])) %>%
      as.data.frame() %>% setNames(c('count')) %>%
      tibble::rownames_to_column('position') %>%
      mutate(sample=sample,position=as.numeric(position))
  }) %>% purrr::reduce(rbind)
  ggplot2::ggplot(data=coverage_tb, aes(x=position, y=count)) +
  # geom_line(linetype="dashed", size=0.2) +
    geom_area() +
  # geom_bar(stat="identity")
  # geom_point() +
    facet_wrap(~sample)
}

# plot avg fpkm for each comparison
# https://github.com/sidbdri/cookiecutter-de_analysis_skeleton/issues/163
plot_scatter_fpkm <- function(results){
  COMPARISON_TABLE %>% pull(comparison) %>% set_names(.) %>% lapply(function(comparison_name){
    x=COMPARISON_TABLE %>% filter(comparison == comparison_name)
    same_in_base<-  SAMPLE_DATA %>%
      filter(!!parse_expr(x$condition_name) == x$condition_base) %>%
      pull(sample_name) %>% str_c('_fpkm',sep = '')
    same_in_condition <-  SAMPLE_DATA %>%
      filter(!!parse_expr(x$condition_name) == x$condition) %>%
      pull(sample_name) %>% str_c('_fpkm',sep = '')

    result_for_plot<-results %>%
    dplyr::select(one_of(c(same_in_base,same_in_condition)),
    padj=str_c(comparison_name, '.padj'),
    l2fc=str_c(comparison_name, '.l2fc'))

    result_for_plot$avg_fpkm_base <- result_for_plot %>% dplyr::select(one_of(same_in_base)) %>%
      mutate(avg_1=rowMeans(.)) %>% pull(avg_1)
    result_for_plot$avg_fpkm_condition <- result_for_plot %>% dplyr::select(one_of(same_in_condition)) %>%
      mutate(avg_1=rowMeans(.)) %>% pull(avg_1)


    start_plot(str_c("scatter_fpkm_", x$comparison))
    p <- result_for_plot %>%
    ggplot(aes(x=avg_fpkm_condition, y=avg_fpkm_base)) +
      geom_point(data = result_for_plot %>% dplyr::filter(padj >= P.ADJ.CUTOFF), shape = 4, colour="black", alpha=0.25) +
      geom_point(data = result_for_plot %>% dplyr::filter(padj < P.ADJ.CUTOFF & l2fc > 0), shape = 4, colour="red") +
      geom_point(data = result_for_plot %>% dplyr::filter(padj < P.ADJ.CUTOFF & l2fc < 0), shape = 4, colour="blue") +
      scale_x_log10() +
      scale_y_log10() +
      xlab(x$condition) + ylab(x$condition_base) +
      theme_classic()
    plot(p)
    end_plot()
    'success'
  })
}

plot_genes_fpkm <- function(result_table,genes,print_fpkm_table = FALSE) {
  num_genes <- length(genes)
  if(num_genes>4){
    warning('more than 4 genes are selected to plot gene fpkm across all samples! label might not be visible in the plot.')
  }

  fpkm_info <- data.frame()
  old.opt<-getOption("ggrepel.max.overlaps")
  options(ggrepel.max.overlaps = 20)

  for (index in seq(num_genes)) {
    plot_name <- str_c('fpkm_plot_', result_table%>%filter(gene==genes[index]) %>% pull(gene_name),sep = '')
    plot_name <- str_replace(plot_name, "-", "_")
    l <- plot_gene_fpkms(gene_identifier = genes[index],result_table = results,debug = FALSE, print_graph = FALSE,
    feature_group=TOPDEGENE_FEATURE_GROUP,plot_feature=TOPDEGENE_PLOT_FEATURE)
    l$graph <- l$graph +
      geom_hline(yintercept = 0, color = "grey") +
      ylim(0, NA) +
      theme_classic() +
      theme(legend.position = "right",
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank())

    assign(plot_name, l$graph)
    fpkm_info <- l$info %>% rbind(fpkm_info)
    plot_statement <- ifelse(index == 1, plot_name, plot_statement %>% str_c(plot_name, sep = '+'))
  }

  options(ggrepel.max.overlaps = old.opt)
  if (num_genes > 1) {
    plot_statement %<>% str_c(" + plot_layout(ncol = ", round(num_genes/2), ", guides = \"collect\")")
  }

  if (print_fpkm_table) {
    fpkm_info %>% print()
  }

  plot_statement %>% parse_expr() %>% eval() %>% print()
}


plot_gene_percentage <- function(count_matrix,gene_set_list){
  total_count_per_sample <- count_matrix %>% as.data.frame() %>% mutate_all(as.numeric) %>%
    summarise(across(everything(), ~ sum(., is.na(.), 0)))

  tb <- total_count_per_sample %>% tidyr::pivot_longer(cols = everything(),names_to='sample',values_to = 'total')
  for(gs in names(gene_set_list)){
    total_goi_count_per_sample <- count_matrix %>% as.data.frame() %>% mutate_all(as.numeric) %>%
      dplyr::filter(rownames(.) %in% gene_set_list[[gs]]) %>%
      summarise(across(everything(), ~ sum(., is.na(.), 0)))
    tb %<>% left_join(total_goi_count_per_sample %>% tidyr::pivot_longer(cols = everything(),names_to='sample',values_to = gs))
  }

  p <- tb %>%
    mutate(others=total-rowSums(across(head(names(gene_set_list),1):tail(names(gene_set_list),1)), na.rm = T)) %>%
    dplyr::select(-total) %>%
    tidyr::pivot_longer(cols = -sample,names_to='type',values_to='count') %>%
    ggplot( aes(fill=type, y=count, x=sample)) +
    geom_bar(position="stack", stat="identity") +
    scale_y_continuous(expand = expansion(mult = c(0, .05))) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  p
}

ontology_find_all_children_terms <- function(term,parent2children){
  if(!is.na(term)){
    children <- parent2children %>% extract2(term) %>% unname()
    return(c(children,lapply(children,ontology_find_all_children_terms,parent2children) %>%
      unlist() %>% discard(is.na) %>% unique)
    )
  }
}
