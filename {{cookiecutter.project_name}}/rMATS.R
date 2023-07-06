source("load_packages.R")
source("meta_data_functions.R")

SPECIES <- "unknown_species"
META_DATA <- stringr::str_c('meta_data_',SPECIES,'.R')
source(META_DATA)

# we are using 4.1.0 version
RMATS_PATH<- '/opt/rmats-turbo/'

# paired-end read or single-end read
{% if cookiecutter.paired_end_read == "yes" %}
rMAT_PARA_t = 'paired'
{% else %}
rMAT_PARA_t = 'single'
{% endif %}

# FDR cutoff in the final result table
P.ADJ.CUTOFF <- 0.05

# read length, this can be detected from the sample fastq file. In most cases, the read length are the same across samples,
# but it could be different.
rMAT_PARA_readLength = system("find -L ./data -name '*.fastq.gz' | xargs -I % bash -c 'zcat % | head -2 | tail -1 | tr -d \"\n\"  | wc -c'",
                              ignore.stderr=T,intern=T) %>% unique()
rMAT_PARA_veriable_readLength=(length(rMAT_PARA_readLength) > 1)

# When read length across samples are not the same, an error will stop the script.
# uncomment the below line and specific a read length manually.
# rMAT_PARA_readLength <- '75'
if(length(rMAT_PARA_readLength) > 1) stop('sample has different read length, plese specify manually!')

# detect the strandness setting from the strand.txt file in the project root folder
libType=switch(system('bash -c "cat strand.txt"',intern=T),
               '0' = "fr-unstranded", '1' = 'fr-firststrand', '2' = 'fr-secondstrand', '-1' = 'unknown')
if(libType=='unknown') stop('rMATS libType is unknown!')

# This is the number of threads used in each rMATS run
# rMAT_PARA_nthread={{cookiecutter.number_total_threads}}
# we hard code this to 1 because it crashes when use more than 1 core
rMAT_PARA_nthread=1

# event with avg count less than avg_count_cutoff will be filter out from the result table
avg_count_cutoff=5

# add some gene level statistic to summary table
SUMMARY_TB %<>% mutate(Up_regulated_gene=integer(),
                      Down_regulated_gene=integer(),
                      D.E.total_gene=integer())

# the generated bash script that contains the command to call rMATS
tmp_script = str_c("rMATS_",SPECIES,".sh",sep = '')

# final result folder
output_folder = file.path('results','differential_expression','de_rmats')
if (!dir.exists(output_folder)) dir.create(output_folder,recursive=TRUE)

# this function prepare the rMATS required files and generated the command to call rMATS
generate_rmats_cmd <- function(sample_data,species,cmp_name){
  x=COMPARISON_TABLE %>% filter(comparison==cmp_name)

  #this filder contains the rMATS output
  top_dir<-file.path("results","rMATS",species,cmp_name)
  tmp_dir<-file.path(top_dir,'tmp')

  # we need to remove this folder because the 'segmentation fault' error
  # when running rmats with more than 1 core in a folder which container pervious rmats runs
  if (dir.exists(top_dir)){
    message('Found pervious rMATS folder for comparison ',cmp_name,' at ', top_dir,'. Removing it to avoid segmentation fault error!')
    unlink(top_dir,recursive = TRUE,force = TRUE)
  }
  dir.create(tmp_dir,recursive=TRUE)

  # create the rMATS sample files
  file.create(file.path(top_dir,"b1.txt"))
  file.create(file.path(top_dir,"b2.txt"))

  # so we know the absolute path
  b1<-file.path(top_dir,"b1.txt") %>% normalizePath()
  b2<-file.path(top_dir,"b2.txt") %>% normalizePath()


  reps<-sample_data %>% tibble::rownames_to_column(var = "tmp_row_names") %>%
    group_by(!!parse_expr(x$condition_name)) %>%
    mutate(bam_file=str_c("results/final_bams/",tmp_row_names,".",species,".bam",sep = '') %>% normalizePath)

  # check if this comparison has sample level pairing
  sample_pairing <- RMATS_SAMPLE_PAIR_TABLE %>% filter(comparison==cmp_name) %>% pull(pair_column)
  if(length(sample_pairing)>0){
    #this will order the samples by the 'pair' column so the bam file are in order
    reps %<>% arrange_at(sample_pairing,.by_group = TRUE)
  }

  reps %<>% summarise(smaple_bams=str_c(bam_file,collapse = ','))

  reps %>%
    filter(!!parse_expr(x$condition_name) == x$condition_base) %>%
    pull(smaple_bams) %>%
    write(b2)

  reps %>%
    filter(!!parse_expr(x$condition_name) == x$condition) %>%
    pull(smaple_bams) %>%
    write(b1)

  # we hard code python3.8 here as rMATS was build in python3.8 on server and in docker,
  # so when working in Rstudio, the system(cmd) will call python3.8 on sidb instead of the sidb default python2
  cmd <- str_c( "python3.8",
                file.path(RMATS_PATH,'rmats.py'),
                "--b1", b1,
                "--b2", b2,
                "--gtf", str_c("data/",dir(path = "data/", pattern = str_c(species,"_ensembl_*"))) %>%
                  list.files(pattern = '.gtf',full.names = T) %>% normalizePath,
                "--od",top_dir %>% normalizePath(),
                "-t", rMAT_PARA_t,
                "--nthread",rMAT_PARA_nthread,
                "--readLength", rMAT_PARA_readLength,
                "--tmp",tmp_dir %>% normalizePath(),
                "--libType",libType,
                ifelse(length(sample_pairing)>0,'--paired-stats',''),
                ifelse(rMAT_PARA_veriable_readLength,'--variable-read-length',''),
                sep = " "
  )
  message("rMATS command generated for comparison ", cmp_name,":\n",cmd)
  cmd
}


# we generate the bash script for rMATS commands
write('
#!/bin/bash
#trap "exit" INT TERM
#trap "kill 0" EXIT
export WORKON_HOME=${HOME}/{{cookiecutter.virtualenv_home}}
export PROJECT_HOME=${HOME}/{{cookiecutter.projects_base}}
source /usr/local/bin/virtualenvwrapper.sh
workon {{cookiecutter.project_name}}
',tmp_script)


# Write command to bash script
COMPARISON_TABLE %>% pull(comparison) %>% walk ( function(x){
  x=COMPARISON_TABLE %>% filter(comparison==x)
  
  sample_data <- SAMPLE_DATA %>%
    tibble::rownames_to_column(var = "tmp_row_names") %>%
    mutate(!!x$condition_name:= factor(!!parse_expr(x$condition_name))) %>%
    filter(!!parse_expr(x$filter)) %>%
    tibble::column_to_rownames(var = "tmp_row_names")

  cmd = sample_data %>% generate_rmats_cmd(SPECIES, x$comparison)
  write(cmd,tmp_script,append = T)
})

cmd = str_c('bash ',tmp_script,sep = '')
# note that this system call is going to use the python outside the VE. Python3.8 is hard coded in the bash command
# so hopefully the dependency on SIDB and docker OUTSIDE the VE has been fulfilled.
system(cmd)

# When finished, we have the rMATS output, we parse the result
source("utility_functions.R")
source("deseq.R")

total_dds_data <- get_total_dds(SAMPLE_DATA,species=SPECIES)
gene_info <- get_gene_info(SPECIES)
gene_lengths <- read_csv(str_c("data/",SPECIES,"_ensembl_",{{cookiecutter.ensembl_version}},"/gene_lengths.csv", sep=""))
results <- total_dds_data %>% get_count_data()
fpkms <- results %>% get_fpkms(gene_lengths, colnames(results) %>% tail(-1), "_fpkm")

results %<>%
  inner_join(fpkms) %>%
  inner_join(gene_info) %>%
  inner_join(gene_lengths)

# We do not need counts for samples for rMATS output
results %<>% dplyr::select(-one_of(SAMPLE_DATA %>% rownames()))

# Load rMATS results and join with result table
COMPARISON_TABLE %>% pull(comparison) %>% walk ( function(x){
  x=COMPARISON_TABLE %>% filter(comparison==x)

  sample_data <- SAMPLE_DATA %>%
    tibble::rownames_to_column(var = "tmp_row_names") %>%
    mutate(!!x$condition_name:= factor(!!parse_expr(x$condition_name))) %>%
    filter(!!parse_expr(x$filter)) %>%
    tibble::column_to_rownames(var = "tmp_row_names")

  top_dir = str_c("results/rMATS/",SPECIES,"/",x$comparison)
  out_dir = str_c(output_folder,SPECIES,x$comparison,sep = '/')

  if (!dir.exists(out_dir)) dir.create(out_dir,recursive=TRUE)

  for (event in c("SE","MXE","A3SS","A5SS","RI")) {
    for(type in c('JC','JCEC')){

      as.event.file=paste(event,'.MATS.',type,'.txt',sep = '')
      # This table contains two 'ID' column, thus removing the redundant one
      # Note that read_tsv might behave differently, early version don't have the name_repir parameter and
      # will rename the duplicated column in a different way. This might cause error!!
      # The name_repair='minimal' will not rename duplicated column, thus allowing read_tsv to read the file
      # without error.
      result.table <- read_tsv(file.path(top_dir,as.event.file),col_types = cols(),name_repair='minimal') %>%
                        dplyr::select(-1) %>%  # we remove the first ID column, as the 2nd ID column position can vary
                        dplyr::select(ID,everything()) %>%
                        suppressWarnings()

      if(nrow(result.table) > 0 ){
              message("For comparison", x$comparison, ", processing file: ",as.event.file,"...")
              # replace ',' with ';' to avoid excel problem in the result
              result.table %<>%
                mutate_at(.vars = vars("IJC_SAMPLE_1","SJC_SAMPLE_1","IJC_SAMPLE_2","SJC_SAMPLE_2","IncLevel1","IncLevel2"),
                          .funs = (~gsub(",",";",.)))

              #workout avg count
              result.table %<>%  dplyr::select(ID,IJC_SAMPLE_1,SJC_SAMPLE_1,IJC_SAMPLE_2,SJC_SAMPLE_2) %>%
                tidyr::unite(tmp,c("IJC_SAMPLE_1","SJC_SAMPLE_1","IJC_SAMPLE_2","SJC_SAMPLE_2"),sep=';') %>%
                mutate(avg_count=strsplit(tmp,';') %>% lapply(as.numeric) %>% lapply(mean) %>% unlist()) %>%
                dplyr::select(-tmp) %>%
                left_join(result.table,by='ID')

              #workout pre comparison avg fpkm
              avg<- fpkms %>% dplyr::select(one_of(str_c(sample_data %>% rownames(),'_fpkm'))) %>%
                mutate(avg_fpkm=rowMeans(.)) %>%
                dplyr::pull(avg_fpkm)

              result.table <- get("results",envir=.GlobalEnv) %>%
                dplyr::select(gene, gene_name, chromosome, description, entrez_id,
                              gene_type, gene_length, max_transcript_length,
                              one_of(str_c(sample_data %>% rownames(),'_fpkm'))) %>%
                mutate(avg_fpkm=avg) %>%
                inner_join(result.table,by=c('gene' = 'GeneID')) %>%
                dplyr::select(everything(),-geneSymbol,-chr)

              #filter low counts
              result.table %<>% filter(avg_count > avg_count_cutoff) %>% dplyr::select(-avg_count)

              #fill summary table
              SUMMARY_TB <- get("SUMMARY_TB", envir = .GlobalEnv) %>% rbind(
                list(
                  Comparison=x$comparison%>% str_c(str_c(event,type,sep = '_'),sep = '_'),
                  DESeq_model_formula = design(get("total_dds_data", envir = .GlobalEnv)) %>% format(),
                  Condition_tested = x$condition_name,
                  Total_number_of_samples_data=sample_data %>% nrow(),
                  Base_level_condition=x$condition_base,
                  Number_of_samples_in_base_level_condition=sample_data %>% filter(!!parse_expr(x$condition_name)==x$condition_base)%>% nrow(),
                  Sample_names_in_base_level_condition=sample_data %>% tibble::rownames_to_column(var = "tmp_row_names") %>% filter(!!parse_expr(x$condition_name)==x$condition_base)%>% pull(tmp_row_names) %>% str_c(collapse = ','),
                  Comparison_level_condition=x$condition,
                  Number_of_samples_in_comparison_level_condition=sample_data %>% filter(!!parse_expr(x$condition_name)==x$condition)%>% nrow(),
                  Sample_names_in_comparison_level_condition=sample_data %>% tibble::rownames_to_column(var = "tmp_row_names") %>% filter(!!parse_expr(x$condition_name)==x$condition)%>% pull(tmp_row_names) %>% str_c(collapse = ','),
                  p.adj.cutoff=P.ADJ.CUTOFF,
                  Up_regulated= result.table %>% filter(FDR < P.ADJ.CUTOFF & IncLevelDifference < 0 ) %>% nrow(),
                  Down_regulated=result.table %>% filter(FDR < P.ADJ.CUTOFF & IncLevelDifference > 0) %>% nrow(),
                  D.E.total=result.table %>% filter(FDR < P.ADJ.CUTOFF) %>% nrow(),
                  Up_regulated_gene = result.table %>% filter(FDR < P.ADJ.CUTOFF & IncLevelDifference < 0)  %>% pull(gene) %>% unique() %>% length(),
                  Down_regulated_gene = result.table %>% filter(FDR < P.ADJ.CUTOFF & IncLevelDifference > 0) %>% pull(gene) %>% unique() %>% length(),
                  D.E.total_gene = result.table %>% filter(FDR < P.ADJ.CUTOFF) %>% pull(gene) %>% unique() %>% length()
                ) %>% as.data.frame(stringsAsFactors = FALSE)
              )
              assign("SUMMARY_TB", SUMMARY_TB,envir = .GlobalEnv)

              # save results
              result.table %>%
                dplyr::select(-PValue,-FDR, everything(),
                              -one_of(str_c(sample_data %>% rownames(),'_fpkm')),
                              -IncFormLen,-SkipFormLen) %>%
                write_csv(file.path(out_dir,paste('AS_results_fpkm_',paste(type,event,SPECIES,sep='_'),".csv",sep='')),na = "")
            }else{
              message("For comparison", x$comparison, ", processing file: ",as.event.file,"... file is empry, skip!")
            } # if(nrow(result.table) < 0)
          } # for (event in c("SE","MXE","A3SS","A5SS","RI"))
        } # for(type in c('JC','JCEC')){
}) # COMPARISON_TABLE %>% pull(comparison) %>% walk ( function(x){

# reformat/save summary table
SUMMARY_TB %>%
  dplyr::select(-DESeq_model_formula) %>%
  write_csv(file.path(output_folder,str_c("AS_summary_",SPECIES,".csv")),na = "")


# # maser package is useful to inspect rMATS results
# library(maser)
# library(rtracklayer)
#
# # function to load rMATS results for a comparison
# load_rmats_results<-function(comparison_name,rmats_res_ftype='JCEC'){
#   x=COMPARISON_TABLE %>% dplyr::filter(comparison==comparison_name)
#   rmats_res_path=file.path('results','rMATS',SPECIES,x$comparison)
#   rmats_res_conditions=c(x$condition, x$condition_base)
#   rmats_res<-maser(rmats_res_path, rmats_res_conditions, ftype = rmats_res_ftype)
#   return(rmats_res)
# }
#
# # load gtf file
# ens_gtf <- str_c("data/",dir(path = "data/", pattern = str_c(SPECIES,"_ensembl_*"))) %>%
#   list.files(pattern = '.gtf',full.names = T) %>% normalizePath %>%
#   rtracklayer::import.gff()
#
# # load rMATS results
# rmats_res<-load_rmats_results('KO_vs_WT')
#
# # inspect rMATS results
# # head(summary(rmats_res, type = "SE")[, 1:8])
#
# # filter out rMATS results with low coverage
# rmats_res_filt <- filterByCoverage(rmats_res, avg_reads = 5)
# maser::volcano(rmats_res_filt, fdr = 0.05, deltaPSI = 0.1, type = "SE")
#
# # get top events
# rmats_res_top <- topEvents(rmats_res_filt, fdr = 0.05, deltaPSI = 0.1)
# maser::dotplot(rmats_res_top, type = "SE")
#
# # get events for a specific gene
# rmats_res_UNC13A <- geneEvents(rmats_res_filt, geneS = "UNC13A", fdr = 0.05, deltaPSI = 0.1)
# maser::display(rmats_res_UNC13A, "SE")
# plotGenePSI(rmats_res_UNC13A, type = "SE", show_replicates = TRUE)
#
# # find event id in the list of events
# rmats_res_UNC13A@SE_events
# # plot specific event
# maser::plotTranscripts(rmats_res_UNC13A, type = "SE", event_id = 86072,gtf = ens_gtf, zoom = FALSE, show_PSI = TRUE)
