source("meta_data_functions.R")

#### Sample definition ####

SAMPLE_NAMES <- c(condition1, condition2, etc) %>%
  outer(c(rep1, rep2, etc), str_c, sep = "todo") %>%
  t %>%
  as.vector

# Remember to order the samples sensibly for output in results spreadsheets (i.e. appropriately
# grouped by some combination of conditions).
SAMPLE_DATA <- data.frame(
  condition = ...,
  species = SPECIES,
  row.names = SAMPLE_NAMES
) %>% 
  tibble::rownames_to_column(var = "tmp_row_names") %>%
  filter(species == !!SPECIES) %>%
  arrange(TODO) %>% 
  tibble::column_to_rownames(var = "tmp_row_names")

# https://github.com/sidbdri/cookiecutter-de_analysis_skeleton/issues/103
# read in picard alignment matrix and extract the MEDIAN_5PRIME_TO_3PRIME_BIAS column for PCA plot later
SAMPLE_DATA$read_bias <- read_median_5prime_to_3prime_bias(SAMPLE_DATA$sample_name)

#### Average FPKM groups ####

# Specify what average FPKMs are needed in the results table. These usually correspond to column names 
# in the SAMPLE_DATA table. There is no limit to the number of groups. For example, for a sample table
# as follows:
# 
# | type      |  condition_2 |
# |--------------------------|
# | cortical  |  correction  |
# | motor     |  correction  |
# | cortical  |  mutant      |
# | motor     |  mutant      |
#
# AVG_FPKM_GROUP = list(c('type','condition_2'),c('type'))
#
# This will generate an average FPKM for all combinations of each the variables in each vector, resulting 
# in 6 columns, i.e.:
# cortical_correction_avg_fpkm, cortical_mutant_avg_fpkm, motor_correction_avg_fpkm, 
# motor_mutant_avg_fpkm, cortical_avg_fpkm, motor_avg_fpkm
AVG_FPKM_GROUP <- list(c(), c())

## we reorder the sample_data table base on the first PCA feature(s) https://github.com/sidbdri/cookiecutter-de_analysis_skeleton/issues/110
SAMPLE_DATA %<>% dplyr::arrange_at(vars(AVG_FPKM_GROUP[[1]]))

#### Comparisons table ####

# An example can be found here:
# https://github.com/sidbdri/cookiecutter-sargasso-de_analysis_skeleton
# If the group column contains more than one group, the results will be saved into different CSVs by group.
# ** N.B. No group name should be a substring of another group name. **
COMPARISON_TABLE <- tribble(
  ~comparison, ~formula, ~condition_name, ~condition, ~condition_base, ~filter, ~group,
  #"P10_Ctx_KO_vs_WT", "~genotype", "genotype", "KO", "WT", "age=='P10' & region=='Ctx'",group_1
)

# This table is for the comparisons that contains a interaction term in their formula
# http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#interactions
INTERACTION_TABLE <- tribble(
  ~comparison, ~condition_name_1, ~condition_1, ~condition_base_1,~condition_name_2, ~condition_2, ~condition_base_2,
  #"LBC_THp_BA17_Good_vs_Poor_interaction","sex","M","F","cognitive","Good","Poor"
)

# tell rMATS which comparison contains paired samples, and the column name in the SAMPLE_DATA which indicates the pairing
# only need to include those comparisons which have sample level pairing
RMATS_SAMPLE_PAIR_TABLE <- tribble(
  ~comparison, ~pair_column,
  # "P10_Ctx_KO_vs_WT","pair"
)

# to make sure comparison name does not start with a number https://github.com/sidbdri/cookiecutter-de_analysis_skeleton/issues/142
check_comparison_name()

# This is to make sure the DESeq2 formulas are ordered so that the 'deciding' 
# condition is the last item in the formula
check_formulas()

# This is to check the samples in each comparison to make sure the 'filter'
# is construct correctly in the COMPARISON_TABLE
check_samples()

#### Sargasso misassignment table ####

# This table specifies the samples used to estimate the misassignment percentage for each condition. Leave the 
# table blank (rows=0) to skip the estimation. Otherwise, two columns will be generated in the results CSV:
#
# <comparison_name>.percentage_misassignment_condition, and 
# <comparison_name>.percentage_misassignment_condition_base
#
# The following explanation takes the HMR example from the paper,
# where we want to estimate the rat reads misassigned from human/mouse
# the reference samples are the co-cultures(HM), the target samples are the triple co-culture(HMR).
#
# Columns: 
# ~condition: matches the condition/condition_base column in the COMPARISON_TABLE.
# ~species_of_interest: "rat"
# ~target_samples: comma separated triple co-culture sample names. "01_HMR,02_HMR"
# ~target_species: comma separated species of the target_samples. "human,mouse,rat"
# ~reference_samples: comma separated co-cultures sample names. "01_HM,02_HM"
# ~reference_species: comma separated species of the reference_samples, one for each refrence sample, separated by ';'.  "human,mouse; human,mouse"
# ~paired: boolean value indicate whether the reference_samples are paired with the target_samples.
#          In the case of paried
#               - the order of the reference sample reflects the pairing.
#               - In the case of using mono-culture samples, the order of each mono species reflects the pairing.
#               - per-reference_sample fpkm(fpkm_mm_hs) will be used.
#          In the case of not paried
#               - the order of the reference sample does not matter.
#               - the reference sample fpkm(fpkm_mm_hs) is calculated by mergeing all reference samples.
MISASSIGNMENT_SAMPLE_REFERENCE_TABLE <- tribble(
  ~comparison_name, ~condition, ~species_of_interest, ~target_samples, ~target_species, ~reference_samples, ~reference_species, ~paired,
  #"tri_vs_mono_pericytes", "triple",SPECIES,"01_HMR,02_HMR","human,mouse,rat","01_M,02_M,01_H,02_H","mouse;mouse;human;human", T,
  #"tri_vs_mono_pericytes", "triple",SPECIES,"01_HMR,02_HMR","human,mouse,rat","01_HM,02_HM","human,mouse;human,mouse", T
)

## This helps you to create part of the MISASSIGNMENT_SAMPLE_REFERENCE_TABLE
# COMPARISON_TABLE %>% pull(comparison) %>% set_names(.) %>%  lapply(
#   function(comparison_name) {
#     x=COMPARISON_TABLE %>% filter(comparison==comparison_name)
#     str_c(x$comparison,
#           x$condition,
#           "SPECIES",
#           SAMPLE_DATA %>% filter(!!parse_expr(x$filter)) %>%
#             filter(!!parse_expr(x$condition_name) == x$condition) %>%
#             pull(sample_name) %>%
#             str_c(collapse = ','),
#           sep = '\",\"') %>% cat('\n',sep = '')
#     str_c(x$comparison,
#           x$condition_base,
#           "SPECIES",
#           SAMPLE_DATA %>% filter(!!parse_expr(x$filter)) %>%
#             filter(!!parse_expr(x$condition_name) == x$condition_base) %>%
#             pull(sample_name) %>%
#             str_c(collapse = ','),
#           sep = '\",\"') %>% cat('\n',sep = '')
#   })

#### PCA and heatmap features ####

# Set up PCA and heatmap plots
FEATURES_FOR_ALL_SAMPLES_PCA <- c('condition')
FEATURES_FOR_ALL_SAMPLES_HEATMAP <- c('condition', 'treatment')

# Set up celltype graph features,  change to match one or more column in SAMPLE_DATA
# The following two vectors must be the same length
# CELLTYPE_PLOT_FEATURE options: "color", "shape"
CELLTYPE_FEATURE_GROUP <- c('condition')
CELLTYPE_PLOT_FEATURE <- c('color')

# Set up top DE gene fpkm graph features, change to match one or more column in SAMPLE_DATA
# The following two vectors must be the same length
# TOPDEGENE_PLOT_FEATURE options: "color", "shape"
TOPDEGENE_FEATURE_GROUP <- c('condition')
TOPDEGENE_PLOT_FEATURE <- c('color')

SUMMARY_TB <- setNames(data.frame(matrix(ncol = 14, nrow = 0)),
                     c("Comparison", 
                       "DESeq_model_formula", 
                       "Condition_tested", 
                       "Total_number_of_samples_data",
                       "Base_level_condition", 
                       "Number_of_samples_in_base_level_condition",
                       "Sample_names_in_base_level_condition",
                       "Comparison_level_condition", 
                       "Number_of_samples_in_comparison_level_condition",
                       "Sample_names_in_comparison_level_condition",
                       "p.adj.cutoff", 
                       "Up_regulated", 
                       "Down_regulated", 
                       "D.E.total"))

# Add columns for reference samples used in calculating misassignment percentage for conditions in comparisons
if (MISASSIGNMENT_SAMPLE_REFERENCE_TABLE %>% nrow() > 0) {
  SUMMARY_TB$Misassignment_samples_in_base_level_condition <- character(0)
  SUMMARY_TB$Misassignment_samples_in_comparison_level_condition <- character(0)
}

#### Cell-type gene markers ####

# This table is used to perform a QC check on cell-type specific genes. The FPKM of each gene in the table 
# will be plotted for each sample. This table is currently generated manually but once we update to R 3.5
# we will be able to use, for example, biomaRt to query orthology:
# https://support.bioconductor.org/p/46475/
GENE_MARKERS <- tribble(
  ~human, ~mouse, ~rat, ~gene_name,~cell_type,
  "ENSG00000066336", "ENSMUSG00000002111", "ENSRNOG00000012172", "Sfpi1", "microglia",
  "ENSG00000204472", "ENSMUSG00000024397", "ENSRNOG00000000853", "Aif1", "microglia",
  "ENSG00000100146", "ENSMUSG00000033006", "ENSRNOG00000011305", "Sox10", "oligodendrocyte",
  "ENSG00000197971", "ENSMUSG00000041607", "ENSRNOG00000016516", "Mbp", "oligodendrocyte",
  "ENSG00000102003", "ENSMUSG00000031144", "ENSRNOG00000059720", "Syp", "neuron",
  "ENSG00000167281", "ENSMUSG00000025576", "ENSRNOG00000003386", "Rbfox3", "neuron",
  "ENSG00000171885", "ENSMUSG00000024411", "ENSRNOG00000016043", "Aqp4", "astrocyte",
  "ENSG00000131095", "ENSMUSG00000020932", "ENSRNOG00000002919", "Gfap", "astrocyte",
  "ENSG00000120156", "ENSMUSG00000006386", "ENSRNOG00000008587", "Tek", "endothelial",
  "ENSG00000184113", "ENSMUSG00000041378", "ENSRNOG00000045811", "Cldn5", "endothelial",
  "ENSG00000067048", "ENSMUSG00000069045", "ENSRNOG00000002501", "Ddx3y", "sex",
  "ENSG00000229807", "ENSMUSG00000086503", "", "Xist", "sex"
)

#### WGCNA ####

# After performing WGCNA, look for correlations between module eigengenes etc. and these variables from
# the SAMPLE_DATA table
WGCNA_EXPERIMENTAL_VARIABLES = c("sex", "treatment")
