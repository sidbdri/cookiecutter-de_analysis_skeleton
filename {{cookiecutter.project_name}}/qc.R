# This function plots the read distribution of a chromosome region of bam files
# It can be used as a pre-'igv' check
# one can test the code with the following:
# >setwd('/srv/data/results/nrf2_ich_jamie_loan/aaaaaaaaarggggghhhhh-959851d2469757607aa3e1f8b0ec1e1cc533cf42/20210805')
# >check_sample_bam(samples=c('14_KO_Ma_ICH','14_KO_Mg_ICH'),species='mouse',chr='2',start=75505857,end=75505957) %>% plot
check_sample_bam <- function(samples = c('14_KO_Ma_ICH','14_KO_Mg_ICH'), 
                             species = 'mouse', 
                             chr = '2', start = 75505857, end = 75510000, 
                             bin_width=1) {
  
  bamRanges <- GRanges(chr, IRanges(start,end))
  
  bamFiles <- file.path('results/final_bams/', str_c(samples, species, 'bam', sep = '.'))
  bamIndexFile <- file.path('results/final_bams/', str_c(samples, species, 'bam.bai', sep = '.'))
  bamExperiment <- list(description = "", created = date())
  bv <- BamViews(bamFiles, bamIndicies = bamIndexFile, bamRanges = bamRanges, bamExperiment = bamExperiment)
  reads <- readGAlignments(bv)
  
  olap1 <- endoapply(reads, subsetByOverlaps, bamRanges)
  olap1 <- lapply(olap1, "seqlevels<-", value = as.character(seqnames(bamRanges)))
  cvg <- endoapply(olap1, coverage,
                   shift = -start(ranges(bamRanges[1])),
                   width = width(ranges(bamRanges[1])))
  
  coverage_tb <- lapply(names(cvg), function(sample){
    a = cvg[[sample]][[1]]
    sapply(c(1:length(runLength(a))),function(x){
      rep(runValue(a)[x],runLength(a)[x])
    }) %>% unlist() %>% as.vector() %>% set_names(c(start(ranges(bamRanges))[1]:end(ranges(bamRanges))[1])) %>%
      as.data.frame() %>% setNames(c('count')) %>%
      tibble::rownames_to_column('position') %>%
      mutate(sample = sample, position = as.numeric(position))
  }) %>% purrr::reduce(rbind)
  
  coverage_tb %>% mutate(position = round(position/bin_width,digits = 0)) %>%
    group_by(sample, position) %>%
    summarise(count = sum(count)) %>% ungroup()  %>%
    # we order the output plots by the samples
    mutate(sample = factor(sample, ordered = T, labels = sub('\\.bam$', '',  basename(bamFiles)), levels = basename(bamFiles))) %>%
    ggplot2::ggplot(aes(x = position, y = count)) +
    # geom_line(linetype="dashed", size=0.2) +
    geom_area() +
    # geom_bar(stat="identity")
    # geom_point() +
    facet_wrap(~sample)
}

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
  if (!is.null(output_data_table_path)) {
    write.csv(pca_data, output_data_table_path)
  }
  
  percent_var <- round(100 * attr(pca_data, "percentVar"))
  
  intgroup.df <- as.data.frame(colData(vst)[, intgroup, drop = FALSE])
  
  if (length(intgroup) > 2) {
    colour_group <- factor(apply(intgroup.df[-ncol(intgroup.df)], 1, paste, collapse = " : "))
    shape_group <- factor(intgroup.df[[ncol(intgroup.df)]])
    p <- pca_data %>% ggplot(aes(PC1, PC2, color = colour_group, shape = shape_group)) + 
      guides(colour = guide_legend(title = "group"),
             shape = guide_legend(title = intgroup[length(intgroup)]))    
  } else if (length(intgroup) == 2) {
    colour_group <- factor(intgroup.df[[1]])
    shape_group <- factor(intgroup.df[[2]])
    p <- pca_data %>% ggplot(aes(PC1, PC2, color = colour_group, shape = shape_group)) + 
      guides(colour = guide_legend(title = intgroup[1]),
             shape = guide_legend(title = intgroup[2]))
  } else {
    group <- colData(vst)[[intgroup]]
    p <- pca_data %>% ggplot(aes(PC1, PC2, color = group)) +
      guides(color = guide_legend(title = intgroup))
  }
  
  p <- p + geom_point(size = 3) +
    xlab(str_c("PC1: ", percent_var[1], "% variance")) +
    ylab(str_c("PC2: ", percent_var[2], "% variance")) +
    theme(legend.position = "right")
  
  if (plot_label)
    p <- p + geom_text(aes(label = !!parse_expr(label_name)), colour = "darkgrey", 
                       position = position_nudge(y = 1), size = 3)
  
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
    melt(id.vars = c("gene"), variable.name = "sample", value.name = "count") 
  
  p <- ggplot(counts, aes(sample, 1 + count)) + 
    geom_violin(aes(fill = sample), scale = "width") + 
    geom_boxplot(width = .1, outlier.shape = NA) + 
    coord_trans(y = "log10") + 
    scale_y_continuous(breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000)) +
    guides(fill = FALSE) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p)
}

plot_genes_fpkm <- function(result_table,genes,print_fpkm_table = FALSE) {
  num_genes <- length(genes)
  if (num_genes > 4) {
    warning('more than 4 genes are selected to plot gene fpkm across all samples! label might not be visible in the plot.')
  }
  
  fpkm_info <- data.frame()
  old.opt <- getOption("ggrepel.max.overlaps")
  options(ggrepel.max.overlaps = 20)
  
  for (index in seq(num_genes)) {
    plot_name <- str_c('fpkm_plot_', result_table %>% filter(gene == genes[index]) %>% pull(gene), sep = '')
    
    l <- plot_gene_fpkms(gene_identifier = genes[index],result_table = result_table,debug = FALSE, print_graph = FALSE,
                         feature_group = TOPDEGENE_FEATURE_GROUP, plot_feature = TOPDEGENE_PLOT_FEATURE)
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
      summarise(is = gene_marker_cell_type[which(like >= 0.5)] %>% paste(collapse = ' / '))
    
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
    dplyr::filter(gene_name == gene_identifier | gene == gene_identifier | entrez_id == gene_identifier)
  
  # get the gene name from the results table
  gene_name <- fpkm_debug$gene_name %>% as.vector()
  
  # turn the tibble into long format: each row is a sample and its fpkm etc
  fpkm_debug_long <- fpkm_debug %>% 
    as_tibble() %>%
    melt(id.var = c("gene", "gene_name", "chromosome", "entrez_id"),
         variable.name = 'sample_meta', value.name = 'fpkm')
  
  # check if featuregroup has been set; if not enter this code
  if (feature_group %>% length() == 0) {
    ## No feature group provided, we are going to plot the FPKM using the sample name and color.
    plot_x <- 'sample_meta'
    plot_label <- 'sample_meta'
    
    if (filter_string != '') {
      fpkm_debug_long %<>% filter(!!parse_expr(filter_string))
    }
    
    # we force the sample order on x axis
    fpkm_debug_long$sample_meta %<>% factor(levels = unique(fpkm_debug_long$sample_meta))
    
    p <- fpkm_debug_long %>%
      ggplot(mapping = aes_string(y = "fpkm", x = plot_x)) + 
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
    fpkm_debug_long %<>% inner_join(SAMPLE_DATA, by = c("sample_meta" = "sample_name"))
    
    # We want to plot only these samples
    if (filter_string != '') {
      fpkm_debug_long %<>% filter(!!parse_expr(filter_string))
    }
    
    # we force the sample order on x axis
    fpkm_debug_long$sample_meta %<>% factor(levels = unique(fpkm_debug_long$sample_meta))
    
    # set up plot, giving x and y variables
    p <- fpkm_debug_long %>% ggplot(mapping = aes_string(y = "fpkm", x = "sample_meta"))
    
    # where set, change features of plot to be by our selected feature group
    for (i in which(plot_feature != '')) {
      switch(plot_feature[i],
             'color' = p <- p + aes_string(color = feature_group[i]),
             'shape' = p <- p + aes_string(shape = feature_group[i])
      )
    }
    
    p <- p + geom_point(size = 3) +
      geom_text_repel(aes(label = !!parse_expr(plot_label)),
                      nudge_x = -0.35, direction = "y", hjust = 0.5, 
                      segment.size = 0.1, size = 3) +
      ggtitle(gene_identifier %>% str_c(gene_name, sep = ':')) + 
      theme(legend.position = "top") +
      theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
  }
  
  if (print_graph) {
    p %>% print()
  }
  
  if (debug) {
    fpkm_debug %>% kable(format = 'markdown', digits = 99) %>% print()
  }
  
  list("info" = fpkm_debug, 'graph' = p)
}

# plot avg fpkm for each comparison
# https://github.com/sidbdri/cookiecutter-de_analysis_skeleton/issues/163
plot_scatter_fpkm <- function(results){
  COMPARISON_TABLE %>% pull(comparison) %>% set_names(.) %>% lapply(function(comparison_name){
    x <- COMPARISON_TABLE %>% filter(comparison == comparison_name)
    same_in_base <- SAMPLE_DATA %>%
      filter(!!parse_expr(x$filter)) %>%
      filter(!!parse_expr(x$condition_name) == x$condition_base) %>%
      pull(sample_name) %>% str_c('_fpkm',sep = '')
    same_in_condition <- SAMPLE_DATA %>%
      filter(!!parse_expr(x$filter)) %>%
      filter(!!parse_expr(x$condition_name) == x$condition) %>%
      pull(sample_name) %>% str_c('_fpkm',sep = '')
    
    result_for_plot <- results %>%
      dplyr::select(one_of(c(same_in_base, same_in_condition)),
                    padj = str_c(comparison_name, '.padj'),
                    l2fc = str_c(comparison_name, '.l2fc'))
    
    result_for_plot$avg_fpkm_base <- result_for_plot %>% dplyr::select(one_of(same_in_base)) %>%
      mutate(avg_1 = rowMeans(.)) %>% pull(avg_1)
    result_for_plot$avg_fpkm_condition <- result_for_plot %>% dplyr::select(one_of(same_in_condition)) %>%
      mutate(avg_1 = rowMeans(.)) %>% pull(avg_1)
    
    result_for_plot %<>% filter(avg_fpkm_condition > 0 & avg_fpkm_base > 0)
    
    start_plot(str_c("scatter_fpkm_", x$comparison))
    p <- result_for_plot %>%
      ggplot(aes(x = avg_fpkm_condition, y = avg_fpkm_base)) +
      geom_point(data = result_for_plot %>% dplyr::filter(is.na(padj)), shape = 4, colour = "grey", alpha = 0.5) +
      geom_point(data = result_for_plot %>% dplyr::filter(padj >= P.ADJ.CUTOFF), shape = 4, colour = "black", alpha = 0.25) +
      geom_point(data = result_for_plot %>% dplyr::filter(padj < P.ADJ.CUTOFF & l2fc > 0), shape = 4, colour = "red") +
      geom_point(data = result_for_plot %>% dplyr::filter(padj < P.ADJ.CUTOFF & l2fc < 0), shape = 4, colour = "blue") +
      scale_x_log10() +
      scale_y_log10() +
      xlab(x$condition) + ylab(x$condition_base) +
      theme_classic()
    plot(p)
    end_plot()
    'success'
  })
}

plot_gene_percentage <- function(count_matrix,gene_set_list,use_percentage=TRUE){
  total_count_per_sample <- count_matrix %>% as.data.frame() %>% mutate_all(as.numeric) %>%
    summarise(across(everything(), ~ sum(., is.na(.), 0)))
  
  tb <- total_count_per_sample %>% tidyr::pivot_longer(cols = everything(), names_to = 'sample', values_to = 'total')
  for (gs in names(gene_set_list)) {
    total_goi_count_per_sample <- count_matrix %>% as.data.frame() %>% mutate_all(as.numeric) %>%
      dplyr::filter(rownames(.) %in% gene_set_list[[gs]]) %>%
      summarise(across(everything(), ~ sum(., is.na(.), 0)))
    tb %<>% left_join(total_goi_count_per_sample %>% tidyr::pivot_longer(cols = everything(), names_to = 'sample', values_to = gs))
  }
  bar_position <- ifelse(use_percentage,'fill','stack')
  
  p <- tb %>%
    mutate(others = total - rowSums(across(head(names(gene_set_list), 1):tail(names(gene_set_list), 1)), na.rm = T)) %>%
    dplyr::select(-total) %>%
    tidyr::pivot_longer(cols = -sample, names_to = 'type',values_to = 'count') %>%
    ggplot(aes(fill = type, y = count, x = sample)) +
    geom_bar(position = bar_position, stat = "identity") +
    scale_y_continuous(expand = expansion(mult = c(0, .05))) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  p
}

get_nuclear_encoded_mitochondrial_genes <- function() {
  mitochondrial_terms <- c('GO:0005739', 
                           ontology_find_all_children_terms('GO:0005739', as.list(GO.db::GOCCCHILDREN))) %>% 
    unique()
  
  topGO::annFUN.org("CC", 
                      mapping = switch(SPECIES, 
                                       mouse = "org.Mm.eg.db", 
                                       rat = "org.Rn.eg.db", 
                                       human = "org.Hs.eg.db"), 
                      ID = "ensembl") %>%
    extract(mitochondrial_terms) %>% 
    unlist() %>% 
    unique()
}

ontology_find_all_children_terms <- function(term, parent2children){
  if (!is.na(term)) {
    children <- parent2children %>% extract2(term) %>% unname()
    return(c(children,lapply(children,ontology_find_all_children_terms,parent2children) %>%
               unlist() %>% discard(is.na) %>% unique)
    )
  }
}

plot_pvalue_distribution <- function(results, pvalue_column) {
  pvals <- results %>%
    filter_at(c(pvalue_column),~ !is.na(.)) %>%
    dplyr::select(pvalue_column)
  
  p <- ggplot(pvals, aes_string(pvalue_column)) + 
    geom_histogram(binwidth = 0.025) 
  
  print(p)
  p
}

plot_volcano <- function(results_table, comparison_name) {
  results_table %<>% 
    mutate(l2fc = !!sym(str_c(comparison_name, ".l2fc")),
           minus_log10_pval = -log10(!!sym(str_c(comparison_name, ".padj"))),
           sig = case_when(
             l2fc < 0 & minus_log10_pval > -log10(P.ADJ.CUTOFF) ~ "down",
             l2fc > 0 & minus_log10_pval > -log10(P.ADJ.CUTOFF) ~ "up",
             TRUE ~ "notsig")) 
  
  top_5_down <- results_table %>% filter(l2fc < 0 & gene_name != "") %>% arrange(desc(minus_log10_pval)) %>% slice_head(n = 5)
  top_5_up <- results_table %>% filter(l2fc > 0 & gene_name != "") %>% arrange(desc(minus_log10_pval)) %>% slice_head(n = 5)
  
  p <- results_table %>% 
    ggplot(aes(x = l2fc, y = minus_log10_pval, color = sig)) +
    geom_point(alpha = 0.5) + 
    geom_hline(yintercept = -log10(P.ADJ.CUTOFF), color = "grey", linetype = "dashed") + 
    geom_vline(xintercept = 0, color = "grey", linetype = "dashed") + 
    scale_color_manual(values = c("blue", "black", "red")) + 
    geom_text_repel(data = top_5_down, aes(label = gene_name), color = "black", force_pull = 0.5) +
    geom_text_repel(data = top_5_up, aes(label = gene_name), color = "black", force_pull = 0.5) +
    xlab("Log2 fold change") + ylab("Adjusted p-value") +
    theme_classic() + 
    theme(legend.position = "none")
  
  print(p)
}