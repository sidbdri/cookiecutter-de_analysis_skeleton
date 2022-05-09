source("load_packages.R")
source("enrichment_tests.R")
source("utility_functions.R")

SPECIES <- "unknown_species"
source(str_c('meta_data_', SPECIES, '.R'))

library("WGCNA")
conflict_prefer("cor", "WGCNA")

RESULTS_DIR <- "results/wgcna/"
dir.create(RESULTS_DIR, recursive = TRUE)

########################################
#              FUNCTIONS               #
########################################

get_dds <- function() {
  count_data <- SAMPLE_NAMES %>%
    map(read_counts, SPECIES) %>%
    purrr::reduce(inner_join) %>%
    tibble::column_to_rownames(var="gene")
  
  dds <- DESeqDataSetFromMatrix(countData = count_data, colData = SAMPLE_DATA, design = ~1)
  dds <- dds[rowSums(counts(dds)) > 1, ]
  dds %>% DESeq
}

plot_soft_threshold_graphs <- function(expression_data) {
  powers = c(c(1:10), seq(from = 12, to = 30, by = 2))
  
  sft <- expression_data %>% t %>% 
    pickSoftThreshold(powerVector = powers, networkType = "signed", blockSize = 100000, verbose = 5)
  
  signed_r2 <- -sign(sft$fitIndices[,3])*sft$fitIndices[,2]
  mean_connectivity <- sft$fitIndices[,5]
  
  plot(powers, signed_r2,
       xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit,signed R^2", type = "n",
       main = paste("Scale independence"));
  text(powers, signed_r2, labels = powers, col = "red")
  abline(h = 0.90, col = "red")
  
  plot(powers, mean_connectivity,
       xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
       main = paste("Mean connectivity"))
  text(powers, mean_connectivity, labels = powers, col = "red")  
}

perform_wgcna <- function(expression_data, power=16) {
  MAX_BLOCK_SIZE <- 20000
  
  net = expression_data %>% t %>% 
    blockwiseModules(power = power, networkType = "signed", 
                     TOMType = "signed", minModuleSize = 30,
                     reassignThreshold = 0, mergeCutHeight = 0.25,
                     numericLabels = TRUE, pamRespectsDendro = FALSE,
                     maxBlockSize = MAX_BLOCK_SIZE, verbose = 5)
  
  print(table(net$colors))
  
  mergedColors = labels2colors(net$colors)
  plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                      "Module colors", dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  
  list(net$colors, orderMEs(net$MEs))
}

get_experimental_variables <- function() {
  vars <- SAMPLE_DATA %>% dplyr::select(all_of(WGCNA_EXPERIMENTAL_VARIABLES))
  
  for (var in WGCNA_EXPERIMENTAL_VARIABLES) {
    vars %<>% mutate(!!var := !!parse_expr(var) %>% as.factor() %>% as.integer())
  }
  
  vars
}

get_eigengene_variable_correlations <- function(module_eigengenes) {
  exp_vars <- get_experimental_variables()
  eigengene_exp_vars_correlation <- module_eigengenes %>% WGCNA::cor(exp_vars, use = 'p')
  correlation_p_vals <- eigengene_exp_vars_correlation %>% corPvalueStudent(SAMPLE_DATA %>% nrow)
  
  list(eigengene_exp_vars_correlation, correlation_p_vals)
}

display_eigengene_variable_correlations <- function(
  eigengene_exp_vars_correlation, correlation_p_vals,
  write_to_file=FALSE) {
  
  textMatrix <- eigengene_exp_vars_correlation %>% 
    signif(2) %>% 
    paste(" (", correlation_p_vals %>% signif(1), ")", sep = "");
  
  dim(textMatrix) <- dim(eigengene_exp_vars_correlation)
  
  if (write_to_file) {
    png(filename = str_c(RESULTS_DIR, "eigengene_correlation.png"),
        width = 800, height = 800)
  }
  
  labeledHeatmap(Matrix = eigengene_exp_vars_correlation,
                 xLabels = get_experimental_variables() %>% names,
                 yLabels = names(module_eigengenes),
                 ySymbols = names(module_eigengenes),
                 colorLabels = FALSE, colors = blueWhiteRed(50),
                 textMatrix = textMatrix, 
                 setStdMargins = FALSE, zlim = c(-1,1),
                 main = paste("Eigengene-experimental variable relationships"))
  
  if (write_to_file) {
    dev.off()
  }
}

get_gene_eigengene_correlations <- function(expression_data, module_eigengenes) {
  gene_eigengene_correlations = expression_data %>% t %>% 
    cor(module_eigengenes, use = "p")
  
  correlation_p_vals = gene_eigengene_correlations %>% 
    corPvalueStudent(SAMPLE_DATA %>% nrow)
  
  gene_eigengene_correlations %<>% as.data.frame
  correlation_p_vals %<>% as.data.frame
  
  modNames <- module_eigengenes %>% names %>% substring(3)
  names(gene_eigengene_correlations) <- "MM" %>% paste(modNames, sep = "")
  names(correlation_p_vals) <- "p.MM" %>% paste(modNames, sep = "")
  
  list(gene_eigengene_correlations, correlation_p_vals)
}

get_gene_variable_correlations <- function(expression_data, variable, var_name) {
  variable %<>% as.data.frame
  names(variable) <- var_name
  
  correlations <- expression_data %>% t %>% cor(variable, use = "p")
  correlation_p_vals = corPvalueStudent(correlations, SAMPLE_DATA %>% nrow)
  
  correlations %<>% as.data.frame 
  correlation_p_vals %<>% as.data.frame
  
  names(correlations) <- paste("GC.", variable %>% names, sep = "")
  names(correlation_p_vals) <- paste("p.GC.", variable %>% names, sep = "")
  
  list(correlations, correlation_p_vals)
}

plot_gene_module_variable_correlations <- function(
  module, modules_to_genes, module_eigengenes, 
  gene_eigengene_correlations, gene_variable_correlations,
  var_name, write_to_file=FALSE) {
  
  mod_names <- module_eigengenes %>% names %>% substring(3)
  
  module_column <- module %>% match(mod_names)
  genes_in_module <- modules_to_genes == module
  
  if (write_to_file) {
    png(filename = str_c(RESULTS_DIR, "M", module, "_", var_name, "_correlation.png"),
        width = 800, height = 800)
  }
  
  verboseScatterplot(gene_eigengene_correlations[genes_in_module, module_column],
                     gene_variable_correlations[genes_in_module, 1],
                     xlab = paste("Module membership in module", module),
                     ylab = paste("Gene correlation for", var_name),
                     main = paste("Gene correlation vs. module membership\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module + 1)
  
  if (write_to_file) {
    dev.off() 
  }
}

plot_module_eigengene_values <- function(module, module_eigengenes, var_name, write_to_file=FALSE) {
  plot_info <- module_eigengenes %>% cbind(SAMPLE_DATA %>% tibble::rownames_to_column(var="name"))
  
  if (write_to_file) {
    png(filename = str_c(RESULTS_DIR, "M", module, "_", var_name, "_eigengene_values.png"),
        width = 800, height = 800)
  }
  
  p <- ggplot(data = plot_info, aes_string(x = str_c(var_name), y = str_c("ME", module))) +
    geom_point() +
    theme_classic()
  
  print(p)
  
  if (write_to_file) {
    dev.off() 
  }
}

plot_module_eigengene_values_per_sample <- function(module, module_eigengenes, write_to_file=FALSE) {
  plot_info <- module_eigengenes %>% cbind(SAMPLE_DATA %>% tibble::rownames_to_column(var="name"))
  
  if (write_to_file) {
    png(filename = str_c(RESULTS_DIR, "M", module, "_eigengene_per_sample.png"),
        width = 800, height = 800)
  }
  
  p <- ggplot(data=plot_info, aes_string(x = "sample_name", y = str_c("ME", module))) +
    geom_point() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = -90, hjust = 0))
  
  print(p)
  
  if (write_to_file) {
    dev.off() 
  }
}

get_genes_for_module <- function(module, modules_to_genes, gene_eigengene_correlations, gene_info) {
  column <- str_c("MM", module)
  
  gene_eigengene_correlations %>% 
    as.data.frame %>% 
    tibble::rownames_to_column(var = "gene") %>% 
    extract(modules_to_genes == module,) %>%
    inner_join(gene_info) %>%
    dplyr::select_("gene", "gene_name", "description", "chromosome", column) %>%
    arrange_(paste0("desc(", column, ")"))
}

perform_go_analyses <- function(significant_genes, expressed_genes, file_prefix) {
  c("BP", "MF", "CC") %>% walk(
    function(x) {
      perform_go_analysis(expressed_genes, significant_genes, x, SPECIES) %>% 
        extract2(1) %>% 
        write_csv(str_c(RESULTS_DIR, file_prefix, "_go_", x %>% tolower, ".csv"),na = "")
    }
  )
}

#####

SAMPLE_NAMES <- rownames(SAMPLE_DATA)

# Use DESeq2 to calculate normalised counts
dds <- get_dds()

# Set a threshold to exclude lowly expressed genes
expressed_genes <- rowSums(counts(dds, norm = T)) > 100

# Perform variance stabilizing transformation on counts, so that they can be input to WGCNA
vst <- dds %>% varianceStabilizingTransformation

# Now exclude lowly-expressed genes
gene_expression <- (vst %>% assay)[expressed_genes, ]

# Plot graphs to help pick a soft-threshold power
gene_expression %>% plot_soft_threshold_graphs

# Perform WGCNA analysis
wgcna_network <- gene_expression %>% perform_wgcna(power = 20)
modules_to_genes <- wgcna_network[[1]]
module_eigengenes <- wgcna_network[[2]]

# Calculate correlations between module eigengenes and experimental variables
eigengene_variable_correlation_data <- get_eigengene_variable_correlations(module_eigengenes)
eigengene_variable_correlations <- eigengene_variable_correlation_data[[1]]
eigengene_variable_correlation_p_vals <- eigengene_variable_correlation_data[[2]]

display_eigengene_variable_correlations(
  eigengene_variable_correlations, eigengene_variable_correlation_p_vals,
  write_to_file = F)

# Calculate correlations between gene expression and module eigengenes
gene_eigengene_correlation_data <- get_gene_eigengene_correlations(
  gene_expression, module_eigengenes)

gene_eigengene_correlations <- gene_eigengene_correlation_data[[1]]

# Calculate correlations between gene expression and experimental variables
gene_variable_correlations <- WGCNA_EXPERIMENTAL_VARIABLES %>% 
  map(function(var) {
    gene_expression %>% 
      get_gene_variable_correlations(SAMPLE_DATA %>% pull(!!parse_expr(var)) %>% as.factor() %>% as.integer(), var) %>% 
      extract2(1)
  }
)

# Construct table of per-gene info for output
mod_numbers <- module_eigengenes %>% names %>% substring(3)

output <- NULL

for (module in seq(0, module_eigengenes %>% colnames %>% length - 1)) {
  genes_in_module <- modules_to_genes == module
  module_column <- str_c("MM", module)
  
  gecs_for_module <- gene_eigengene_correlations[genes_in_module, ] %>%
    dplyr::select_(module_column) %>%
    rename_("eigengene_cor" = module_column)
  
  if (is.null(output)) {
    output <- gecs_for_module 
  } else {
    output %<>% rbind(gecs_for_module) 
  }
}

output %<>% tibble::rownames_to_column(var = "gene") %>%
  inner_join(gene_expression %>% 
               as.data.frame %>% 
               tibble::rownames_to_column(var = "gene") %>%
               cbind(modules_to_genes) %>% 
               mutate(module = str_c(modules_to_genes)) %>% 
               dplyr::select(gene, module))

for (gv_cor in gene_variable_correlations) {
  output %<>% inner_join(gv_cor %>% tibble::rownames_to_column(var = "gene"))
}

gene_info <- get_gene_info(SPECIES)
output %<>% inner_join(gene_info) 

output %<>% 
  dplyr::select(gene, gene_name, description, chromosome, module, 
                eigengene_cor, starts_with("GC")) %>% 
  write_csv(str_c(RESULTS_DIR, "genes_to_modules.csv"),na = "")

# For each module:
# (i)   Scatter plots of gene-variable correlations against gene-eigengene correlations
# (ii)  Plot of per-condition eigengene values 
# (iii) Perform GO analyses for genes in module

gene_universe <- expressed_genes %>% as.data.frame %>% tibble::rownames_to_column(var="gene")
colnames(gene_universe) <- c("gene", "expressed")
gene_universe %<>% filter(expressed == TRUE)

for (module in seq(0, module_eigengenes %>% colnames %>% length - 1)) {
  for (i in seq_along(WGCNA_EXPERIMENTAL_VARIABLES)) {
    plot_gene_module_variable_correlations(
      module, modules_to_genes, module_eigengenes,
      gene_eigengene_correlations, gene_variable_correlations[[i]],
      WGCNA_EXPERIMENTAL_VARIABLES[i], write_to_file = TRUE)
    
    plot_module_eigengene_values(module, module_eigengenes, WGCNA_EXPERIMENTAL_VARIABLES[i], write_to_file = TRUE)
  }

  plot_module_eigengene_values_per_sample(module, module_eigengenes, write_to_file = TRUE)
  
  get_genes_for_module(module, modules_to_genes, gene_eigengene_correlations, gene_info) %>%
    perform_go_analyses(gene_universe, str_c("M", module))
}
