SALMON <- "salmon"

run_deseq <- function(comparison_name, results_tbl, fpkms, species, qSVA, use_tx=FALSE, tx_level=FALSE, ...){

  res <- get_res(comparison_name, fpkms, species, qSVA = qSVA, use_tx=use_tx, tx_level=tx_level)

  results_tb <- results_tbl %>%
    left_join(res[[1]]) %>%
    dplyr::rename(!!str_c(comparison_name, '.l2fc') := log2FoldChange,
                  !!str_c(comparison_name, '.raw_l2fc') := raw_l2fc,
                  !!str_c(comparison_name, '.stat') := stat,
                  !!str_c(comparison_name, '.pval') := pvalue,
                  !!str_c(comparison_name, '.padj') := padj)

  # for interaction we remove raw_l2fc column
  if (comparison_name %in% INTERACTION_TABLE$comparison) {
    results_tb %<>% dplyr::select(-str_c(comparison_name, '.raw_l2fc'))
  }

  P <- NULL
  if (MISASSIGNMENT_PERCENTAGE) {
    P <- get_misassignment_percentages(comparison_name, gene_lengths)

    if (!is.na(P$condition_reference_samples)) {
      results_tb %<>% left_join(
        P$P_condition %>%
          dplyr::select(gene, !!str_c(comparison_name, '.perc.', COMPARISON_TABLE %>%
            filter(comparison == comparison_name) %>%
            pull(condition)) := p))
    }

    if (!is.na(P$condition_base_reference_samples)) {
      results_tb %<>% left_join(
        P$P_condition_base %>%
          dplyr::select(gene,!!str_c(comparison_name, '.perc.', COMPARISON_TABLE %>%
            filter(comparison == comparison_name) %>%
            pull(condition_base)) := p))
    }

    res$summary_tb_row %<>% as.data.frame() %>%
      mutate(Misassignment_samples_in_comparison_level_condition = P$condition_reference_samples) %>%
      mutate(Misassignment_samples_in_base_level_condition = P$condition_base_reference_samples)

  }

  p_plot <- plot_pvalue_distribution(results_tb, str_c(comparison_name,'.pval'))

  ## return the results and merge them later
  list(comparison_name = comparison_name,
       res = res$res,
       dds = res$dds,
       vst = res$vst,
       results_tb = results_tb,
       summary_tb = res$summary_tb_row,
       p_plot = p_plot,
       misassignment_percentage = P)

}

# Get differential expression results for a given comparison name
get_res <- function(comparison_name, tpms, species, qSVA=FALSE,
                    use_tx=FALSE, tx_level=FALSE, alpha=0.05) {

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
    txi <- get_tximport(sample_data, species, tx_level)
    dds <- DESeqDataSetFromTximport(txi, sample_data, x$formula %>% as.formula())
    dds <- dds[rowSums(counts(dds)) > 1, ]
    
    betaPrior <- TRUE
    if (any(dds %>% design() %>% terms.formula %>% attr( "order") > 1)) betaPrior <- FALSE
    
    dds <- DESeq(dds, betaPrior = betaPrior)
  }else{
    dds <- sample_data %>%
      row.names() %>%
      map(read_counts, species) %>%
      purrr::reduce(inner_join) %>%
      remove_gene_column() %>%
      get_deseq2_dataset(sample_data, design_formula = x$formula %>% as.formula(), qSVA=qSVA)
  }

  res <- dds %>%
    get_deseq2_results(x$condition_name, x$condition, x$condition_base, alpha = alpha) %>%
    left_join(dds %>% get_raw_l2fc(sample_data, expr(!!sym(x$condition_name) == !!(x$condition))))

  if (use_tx & tx_level) {
    # res contains transcripts, rather than genes
    res %<>% dplyr::rename(transcript = gene)
  }

  # Create a PCA plot and heatmap for this comparison
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

  # plot pca for all feature in the sample_data table
  # https://github.com/sidbdri/cookiecutter-de_analysis_skeleton/issues/198
  plot_name=str_c("pca_all_features_",comparison_name,sep = '')
  start_plot(plot_name,num_plots = num_features)

  if (global_exists(plot_name)) {
    rm_global(plot_name)
  }

  sample_data %>% dplyr::select(-species,-sample_name) %>% colnames() %>%
    walk(function(feature) {
      # only including features which aren't the same for all samples in that comparison
      if(sample_data %>% pull(feature) %>% unique() %>% length() > 1){
        vst %>%
          plot_pca(intgroup = c(feature), FALSE) %>%
          add_to_patchwork(plot_var_name = plot_name)
      }
    })

  get(plot_name) %>% print
  end_plot()

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

  list(res=res, dds=dds, vst=vst, summary_tb_row=summary_tb_row)
}

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

.get_avg_fpkm_table <- function(){
  lapply(AVG_FPKM_GROUP,function(g){
    SAMPLE_DATA %>%  tibble::rownames_to_column(var = "tmp_row_names") %>%
      group_by(.dots=g) %>%
      dplyr::summarise(samples=str_c(tmp_row_names, '_fpkm', sep = '', collapse = ','),.groups = 'drop') %>%
      tidyr::unite('avg_name', g, sep='_')
  }) %>% reduce(rbind)
}

get_avg_fpkm <- function(fpkms) {
  avg_table <- .get_avg_fpkm_table()
  for (avg in avg_table$avg_name) {
    samples <- avg_table %>%
      filter(avg_name == avg) %>%
      pull(samples) %>%
      str_split(',') %>%
      unlist()

    avg_fpkm <- fpkms %>% dplyr::select(one_of(samples)) %>%
      mutate(avg=rowMeans(.)) %>%
      dplyr::pull(avg)

    fpkms %<>% mutate(!!str_c(avg,'_avg_fpkm',sep='') := avg_fpkm)
  }

  fpkms %>% dplyr::select(gene, contains('avg'))
}

get_avg_tpm <- function(tpms) {
  get_avg_fpkm(tpms)
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
    avg_tb<-.get_avg_fpkm_table()
    avg_to_exclude<- lapply(avg_tb$avg_name,function(x){
      samples_in_avg <- avg_tb %>% filter(avg_name==x) %>% pull(samples) %>% strsplit(',') %>% extract2(1) %>% gsub('_fpkm|_tpm','',x=.)
      (samples_in_avg %in% samples_to_include) %>% all
    }) %>% unlist() %>% `!` %>% filter(.data=avg_tb) %>% pull(avg_name)

    if(use_tx){
      results %>%
        dplyr::select(
          any_of(c('transcript', 'transcript_length', 'gene', 'number_of_transcript','gene_length','max_transcript_length')),
          gene_name, chromosome, description, entrez_id, gene_type, everything(),
          -dplyr::contains("_tpm"), -dplyr::ends_with(".stat"),
          -matches(n_comparisons), -all_of((samples_to_exclude))) %>%
        write_csv(file.path(DE_OUT_DIR, str_c(g, "_deseq2_results_count_", SPECIES, ".csv")), na="")
      
      tpm_output <- results %>%
        dplyr::select(
          any_of(c('transcript', 'transcript_length', 'gene', 'number_of_transcript','gene_length','max_transcript_length')),
          gene_name, chromosome, description, entrez_id, gene_type,
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
      
      tpm_output %>% write_csv(file.path(DE_OUT_DIR, str_c(g ,"_deseq2_results_tpm_", SPECIES, ".csv")), na="")
      
      SUMMARY_TB %>% filter(Comparison %in% comparisons) %>%
        write_csv(file.path(DE_OUT_DIR, str_c(g ,"_de_summary_", SPECIES, ".csv")), na="")
    }else{
      # non-tx
      results %>%
        dplyr::select(
          gene, gene_name, chromosome, description, entrez_id, gene_type,
          gene_length, max_transcript_length, everything(),
          -dplyr::contains("_fpkm"), -dplyr::ends_with(".stat"),
          -matches(n_comparisons), -(samples_to_exclude)) %>%
        write_csv(file.path(DE_OUT_DIR, str_c(g, "_deseq2_results_count_", SPECIES, ".csv")), na="")

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

      fpkm_output %>% write_csv(file.path(DE_OUT_DIR, str_c(g ,"_deseq2_results_fpkm_", SPECIES, ".csv")), na="")

      SUMMARY_TB %>% filter(Comparison %in% comparisons) %>%
        write_csv(file.path(DE_OUT_DIR, str_c(g ,"_de_summary_", SPECIES, ".csv")), na="")
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

#### Transcript-level D. E. analyses ####

get_transcripts_to_genes <- function(species = 'human') {
  str_c("data/", species , "_ensembl_{{cookiecutter.ensembl_version}}/tx2gene.tsv") %>%
    read_delim(delim = " ", col_names = c("transcript", "gene"), col_types = cols())
}

get_tximport <- function(sample_data, species, tx_level = TRUE) {
  quant_dirs <- sample_data %>%
    tibble::rownames_to_column(var = "tmp") %>%
    pull("tmp")

  quant_files <- str_c("results/salmon_quant/", species, "/",  quant_dirs, "/quant.sf")
  names(quant_files) <- quant_dirs

  txi <- tximport(quant_files, type = "salmon", txOut = tx_level,
                  tx2gene = get_transcripts_to_genes(species), dropInfReps = TRUE)

  txi$Length <- read.csv(quant_files[1],sep = '\t',stringsAsFactors = FALSE) %>%
    dplyr::select(id = 1, length = 2) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(var = 'id')

  txi
}

get_total_dds_tximport <- function(sample_data, species, tx_level = TRUE, design_formula=~1) {
  txi <- get_tximport(sample_data, species, tx_level)

  total_dds <- DESeqDataSetFromTximport(txi, sample_data, design_formula)
  total_dds <- DESeq(total_dds)

  total_dds
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

# function to generate plot to check the interaction effect
# like the one https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#interactions
# plot_interaction(dds,gene='ENSMUSG00000035486',group=c("condition","genotype"))
plot_interaction <- function(dds, gene='ENSMUSG00000035486',group=c("condition","genotype")){
  d <- plotCounts(dds, gene=gene, intgroup=group, returnData=TRUE)
  ggplot(d, aes_string(x=group[2], y="count",group=group[1])) +
    facet_wrap(as.formula(paste("~", group[1]))) +
    geom_point(position=position_jitter(w=0.1,h=0)) + ggtitle('ENSMUSG00000035486') +
    xlab(group[2]) + ylab("log2(counts+1)") +
    stat_summary(fun.y=mean, geom="line", colour="red", size=0.8)
}
