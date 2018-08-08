source("load_packages.R")
source("common_functions.R")

species={{cookiecutter.species}}

SAMPLE_NAMES <- c(condition1, condition2, etc) %>%
  outer(c(rep1, rep2, etc), str_c, sep="todo") %>%
  t %>%
  as.vector

# n.b. Ensure that conditions to be used in GSA comparisons are factors with
# the correct base level set.
SAMPLE_DATA <- data.frame(
  condition=...,
  species=species,
  row.names=SAMPLE_NAMES%>%str_c(species,sep = '.')
)

SAMPLE_DATA %<>% tibble::rownames_to_column(var = "tmp_row_names") %>%
                 filter(species==!!species) %>%
                 tibble::column_to_rownames(var = "tmp_row_names")


# To specify what avg fpkms are needed in the result table
# This usually corespond to column names in the SAMPLE_DATA table
# There is no limit of the number of groups
# Example:
# For a table as follow:
# | type      |  condition_2 |
# |--------------------------|
# | cortical  |  correction  |
# | motor     |  correction  |
# | cortical  |  mutant      |
# | motor     |  mutant      |
# AVG_FPKM_GROUP = c('type','condition_2')
# This will generate avg fpkm for all combinations of the two variable in the 
# SAMPLE_DATA table, resulting 4 columns 
# cortical_correction_avg_fpkm, cortical_mutant_avg_fpkm, motor_correction_avg_fpkm, motor_mutant_avg_fpkm
AVG_FPKM_GROUP = c()

#example can be found https://github.com/sidbdri/cookiecutter-sargasso-de_analysis_skeleton
COMPARISON_TABLE<-tribble(
~comparison, ~formula, ~condition_name, ~condition, ~condition_base, ~filter,
#"P10_Ctx_KO_vs_WT", "~genotype", "genotype", "KO", "WT", "age=='P10' & region=='Ctx'",
)

SUMMARY_TB<-setNames(data.frame(matrix(ncol = 14, nrow = 0)),
                     c("Comparison", "DESeq_model_formula", "Condition_tested", "Total_number_of_samples_data",
                       "Base_level_condition", "Number_of_samples_in_base_level_condition",
                       "Sample_names_in_base_level_condition",
                       "Comparison_level_condition", "Number_of_samples_in_comparison_level_condition",
                       "Sample_names_in_comparison_level_condition",
                       "p.adj.cutoff", "Up_regulated", "Down_regulated", "D.E.total"))


#####