set_global <- function(value, variable) {
  assign(variable, envir = .GlobalEnv, value)
}

get_global <- function(variable) {
  get(variable, envir = .GlobalEnv)  
}

global_exists <- function(variable) {
  exists(variable, where = .GlobalEnv)
}

rm_global <- function(variable) {
  rm(list = c(variable), envir = .GlobalEnv)
}

read_counts <- function(sample, species) {
  counts_file_name <- str_c("results/read_counts/", sample, ".", species,".counts")
  counts_file_name %>% read_tsv(col_names = c("gene", str_c(sample)))
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
    filter(row_number() == 1) %>% 
    ungroup
}

get_gene_lengths <- function(species) {
  species %>% 
    str_c("data/", ., "_ensembl_{{cookiecutter.ensembl_version}}/gene_lengths.csv") %>%
    read_csv
}

start_plot <- function(prefix, width = 12, height = 12, path = GRAPHS_DIR, num_plots = 1) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }

  .adjust_pdf_size <- function(num_plots) {
    num_features <- num_plots
    num_row <- sqrt(num_features) %>% ceiling()
    num_column <- num_features/num_row %>% ceiling()
    c('width' = max(num_row/2, 1), 'height' = max(num_column/2, 1))
  }

  sf <- .adjust_pdf_size(num_plots)

  if (PLOT_TO_FILE) {
    file.path(path, str_c(prefix, "_", SPECIES, ".pdf")) %>%
      pdf(width = width*sf['width'], height = height*sf['height'])
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

load_rs_data <- function(file = 'results/Rworkspace/diff_expr.RData') {
  if (!file.exists(file)) {
    stop(file,' does not exist!')
  }
  
  load(file, envir = .GlobalEnv)
}