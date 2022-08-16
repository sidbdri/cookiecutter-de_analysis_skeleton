check_formulas <- function() {
  for (r in COMPARISON_TABLE %>% rownames()) {
    row <- COMPARISON_TABLE[r,]
    f <- row$formula %>% as.formula() %>% terms()
    condition <- row$condition_name

    # we check factor in the design has more than 1 levels
    sample_data <- SAMPLE_DATA %>% filter(!!parse_expr(row$filter))
    f_levels_chec <- sapply(labels(f)[which(attr(f,"order")==1)], function(x) {sample_data %>% pull(x) %>% unique()}, simplify = F)
    if (any(sapply(f_levels_chec, length) == 1)) {
      print(f_levels_chec)
      stop('comparison ', row$comparison,': contrasts can be applied only to factors with 2 or more levels.')
    }

    # we check if design matrix if full rank
    if (!model.matrix(row$formula %>% as.formula(), data = sample_data) %>% is.fullrank) {
      print(sample_data)
      print(model.matrix(row$formula %>% as.formula(), data = sample_data))
      stop('comparison ', row$comparison, ' design matrix is not full rank..., formula: ', row$formula %>% as.formula())
    }
    
    # if there is an interaction, we check if we have the interaction detail in the INTERACTION_TABLE,
    # otherwise, we check if we have the deciding condition correctly setup
    if (any(f %>% attr( "order") > 1)) {
      if (!(row$comparison %in% INTERACTION_TABLE$comparison)) {
        print(row)
        stop("interaction formula was detected but interaction details cannot be found in the INTERACTION_TABLE.")
      }
    } else {
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
    x <- COMPARISON_TABLE %>% filter(comparison == comparison_name)
    sample_data <- SAMPLE_DATA %>% filter(!!parse_expr(x$filter))
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
  invalid <- COMPARISON_TABLE %>% pull(comparison) %>% grepl('^[0-9]', x = ., perl = T)
  if (any(invalid)) {
    message('Comparison cannot start with a number. The following comparison names are invalid: ')
    cat(COMPARISON_TABLE %>% pull(comparison) %>% extract(which(invalid)) %>% paste(collapse = '\n'))
  }

  invalid <- COMPARISON_TABLE %>% pull(comparison) %>% grepl('+-', x = ., perl = T)
  if (any(invalid)) {
    message('It is not recommanded to have +/- symbol in comparison name. The following comparison names are not recommanded: ')
    cat(COMPARISON_TABLE %>% pull(comparison) %>% extract(which(invalid)) %>% paste(collapse = '\n'))
  }
}

# read picard qc matrix for rna degradation
read_median_5prime_to_3prime_bias <- function(samples = c()) {
  PICARD_ALIGNMENT_MATRIX <- file.path('results/', 'alignment_metrics', SPECIES)
  MEDIAN_5PRIME_TO_3PRIME_BIAS <- samples %>% lapply(function(sample) {
    read_tsv(file.path(PICARD_ALIGNMENT_MATRIX,str_c(sample,'.txt')),
             comment = '#',
             n_max = 1,
             col_type = cols()) %>%
      pull(MEDIAN_5PRIME_TO_3PRIME_BIAS)
  }) %>% unlist() %>% setNames(samples)
}