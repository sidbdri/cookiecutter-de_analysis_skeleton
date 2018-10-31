source("load_packages.R")
source("common_functions.R")

SPECIES="{{cookiecutter.species}}"

SAMPLE_NAMES <- c(condition1, condition2, etc) %>%
  outer(c(rep1, rep2, etc), str_c, sep="todo") %>%
  t %>%
  as.vector

SAMPLE_DATA <- data.frame(
  condition=...,
  species=SPECIES,
  row.names=SAMPLE_NAMES %>% str_c(species, sep = '.')
) %>% 
  tibble::rownames_to_column(var = "tmp_row_names") %>%
  filter(species==!!SPECIES) %>%
  tibble::column_to_rownames(var = "tmp_row_names")

# Specify what average FPKMs are needed in the results table. These usually 
# correspond to column names in the SAMPLE_DATA table. There is no limit of the 
# number of groups. For example, for a sample table as follows:
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
# This will generate an average fpkm for all combinations of each the variables 
# in each vector resulting in 6 columns, i.e.:
# cortical_correction_avg_fpkm, cortical_mutant_avg_fpkm, motor_correction_avg_fpkm, 
# motor_mutant_avg_fpkm, cortical_avg_fpkm, motor_avg_fpkm
AVG_FPKM_GROUP = list(c(),c())

# An example can be found here:
# https://github.com/sidbdri/cookiecutter-sargasso-de_analysis_skeleton
COMPARISON_TABLE<-tribble(
~comparison, ~formula, ~condition_name, ~condition, ~condition_base, ~filter,
#"P10_Ctx_KO_vs_WT", "~genotype", "genotype", "KO", "WT", "age=='P10' & region=='Ctx'",
)

# This is to make sure the DESeq2 formulas are ordered so that the 'deciding' 
# condition is the last item in the formula
check_formulas()


# This table specificed the samples used to estimate the misassign percentage for each condition.
# Leave the table blank (row=0) to skip the estimation. Otherwise, two columns will be generated
# in the result csv,
# <comparison_name>.percentage_misassignment_condition and <comparison_name>.percentage_misassignment_condition_base
# The ~condition column matchs the condition/condition_base column in the COMPARISON_TABLE
# The ~misassignment_samples_filter will be used in dplyr::filter to find out reference samples for the estimation.
# The ~misassignment_sample_species is the species in the reference samples.
# For example,
# "AMCon", "cells=='NA' & treatment=='Con'", "human,mouse",
# This row specify that, for the rat sample with condition 'AMCon',
# we will use the human/mouse samples to estimate the misassign percentage.
MISASSIGNMENT_SAMPLE_REFERENCE_TABLE <- tribble(
  ~condition, ~misassignment_samples_filter, ~misassignment_sample_species,
 #"AMCon", "cells=='NA' & treatment=='Con'", "human,mouse",
)

# Set up PCA and heatmap plots
PCA_FEATURE<-c('condition')
HEAT_MAP_FEATURE<-c('condition','treatment')

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

# add column for reference samples used in calculating misassignment percentage for conditions in comparisons
if(MISASSIGNMENT_SAMPLE_REFERENCE_TABLE %>% nrow() >0){
  SUMMARY_TB$Misassignment_samples_in_base_level_condition  <- character(0)
  SUMMARY_TB$Misassignment_samples_in_comparison_level_condition  <- character(0)
}