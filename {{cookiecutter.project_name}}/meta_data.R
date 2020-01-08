source("load_packages.R")
source("common_functions.R")

SPECIES <- "unknown_species"
P.ADJ.CUTOFF <- 0.05
NUM_CORES <- 30

SAMPLE_NAMES <- c(condition1, condition2, etc) %>%
  outer(c(rep1, rep2, etc), str_c, sep="todo") %>%
  t %>%
  as.vector

SAMPLE_DATA <- data.frame(
  condition=...,
  species=SPECIES,
  row.names=SAMPLE_NAMES
) %>% 
  tibble::rownames_to_column(var = "tmp_row_names") %>%
  filter(species==!!SPECIES) %>%
  tibble::column_to_rownames(var = "tmp_row_names")



# https://github.com/sidbdri/cookiecutter-de_analysis_skeleton/issues/103
# read in picard alignment matrix and extract the MEDIAN_5PRIME_TO_3PRIME_BIAS column for PCA plot later
SAMPLE_DATA$read_bias=read_median_5prime_to_3prime_bias(SAMPLE_DATA$sample_name)

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

# An example can be found here:
# https://github.com/sidbdri/cookiecutter-sargasso-de_analysis_skeleton
# If the group column contains more than one group, the results will be saved into different CSVs by group.
COMPARISON_TABLE <- tribble(
  ~comparison, ~formula, ~condition_name, ~condition, ~condition_base, ~filter, ~species, ~group,
  #"P10_Ctx_KO_vs_WT", "~genotype", "genotype", "KO", "WT", "age=='P10' & region=='Ctx'",group_1
)

# This is to make sure the DESeq2 formulas are ordered so that the 'deciding' 
# condition is the last item in the formula
check_formulas()

# This is to check the samples in each comparison to make sure the 'filter'
# is construct correctly in the COMPARISON_TABLE
check_samples()

# This table specifies the samples used to estimate the misassignment percentage for each condition. Leave the 
# table blank (rows=0) to skip the estimation. Otherwise, two columns will be generated in the results CSV:
#
# <comparison_name>.percentage_misassignment_condition, and 
# <comparison_name>.percentage_misassignment_condition_base
#
# The ~condition column matches the condition/condition_base column in the COMPARISON_TABLE.
# The ~misassignment_samples_filter will be used in dplyr::filter to find reference samples for the estimation.
# ~reference_species are the species in the reference samples.
#
# For example,
# "AMCon", "cells=='NA' & treatment=='Con'", "human,mouse",
#
# This row specifies that, for the rat samples in condition 'AMCon', we will use the human/mouse samples 
# specified to estimate the misassignment percentage.
MISASSIGNMENT_SAMPLE_REFERENCE_TABLE <- tribble(
  ~condition, ~misassignment_samples_filter, ~reference_species,
 #"AMCon", "cells=='NA' & treatment=='Con'", "human,mouse",
)

# Set up PCA and heatmap plots
PCA_FEATURE <- c('condition')
HEAT_MAP_FEATURE <- c('condition', 'treatment')

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
  "ENSG00000133636", "ENSMUSG00000019890", "ENSRNOG00000004179", "Nts", "Dorsal",
  "ENSG00000003137", "ENSMUSG00000063415", "ENSRNOG00000015076", "Cyp26b1", "Dorsal",
  "ENSG00000185551", "ENSMUSG00000030551", "ENSRNOG00000010308", "Nr2f2", "ventral_hippocampus",
  "ENSG00000140848", "ENSMUSG00000034361", "ENSRNOG00000043286", "Cpne2", "ventral_hippocampus",
)
