source('meta_data.R')

species="{{cookiecutter.species}}"

RMATS_402 = '/opt/rMATS.4.0.2'
RMATS_325 = '/home/xinhe/Projects/rMATS3/MATS/rMATS_Paired.sh'

rMAT_PARA_t = 'paired'
rMAT_PARA_nthread=8
rMAT_PARA_tstat=8
rMAT_PARA_readLength=75
rMAT_PARA_cstat = 0.05

generate_rmats_count_cmd <- function(sample_data,species,cmp_name){

  x=COMPARISON_TABLE %>% filter(comparison==cmp_name)

  top_dir<-str_c("results/rMATS/",species,"/",cmp_name)
  if (!dir.exists(top_dir)) {
    dir.create(top_dir,recursive=TRUE)
  }
  
  file.create(str_c(top_dir,"/b1.txt"))
  file.create(str_c(top_dir,"/b2.txt"))
  
  b1<-str_c(top_dir,"/b1.txt") %>% normalizePath()
  b2<-str_c(top_dir,"/b2.txt") %>% normalizePath()
  
  reps<-sample_data %>% 
    group_by(!!parse_expr(x$condition_name)) %>% 
    mutate(bam_file=str_c("results/final_bams/",sample_name,".",species,".bam",sep = '') %>% normalizePath) %>%
    summarise(replicates=str_c(bam_file,collapse = ',')) 
  
  reps %>% 
    filter(!!parse_expr(x$condition_name) == x$condition_base) %>% 
    pull(replicates) %>%
    write(b1)
  
  reps %>% 
    filter(!!parse_expr(x$condition_name) == x$condition) %>% 
    pull(replicates) %>%
    write(b2)
  
  cmd_count <- str_c( str_c("cd", RMATS_402,sep=' '), " && python rmats.py",
               "--b1", b1, 
               "--b2", b2,
               "--gtf", str_c("data/",dir(path = "data/", pattern = str_c(species,"_ensembl_*"))) %>% 
                 list.files(pattern = '.gtf',full.names = T) %>% normalizePath,
               "--od",top_dir %>% normalizePath(),
               "-t", rMAT_PARA_t,
               "--nthread",rMAT_PARA_nthread,
               "--tstat",rMAT_PARA_tstat,
               "--readLength", rMAT_PARA_readLength,
               "--statoff",
               "--cstat", rMAT_PARA_cstat,
               "--libType fr-unstranded",sep = " "
  )
  
  cmd_stat <- generate_rmats_stat_cmd(species,cmp_name)

  cmd = str_c("(", cmd_count," && ",cmd_stat," ) &",sep ="")

  cmd
}

generate_rmats_stat_cmd <- function(species,comparison){
    top_dir<-str_c("results/rMATS/",species,"/",comparison) %>% normalizePath()
    countFiles<-top_dir %>% str_c(c("JCEC.raw.input.A3SS.txt","JC.raw.input.A3SS.txt",
                                    "JCEC.raw.input.MXE.txt","JC.raw.input.MXE.txt",
                                    "JCEC.raw.input.SE.txt","JC.raw.input.SE.txt",
                                    "JC.raw.input.A5SS.txt","JCEC.raw.input.A5SS.txt",
                                    "JC.raw.input.RI.txt","JCEC.raw.input.RI.txt"),sep = '/')
    cmd=c()
    for (cf in countFiles) {
        #cf = '/home/xinhe/Projects/als_test/results/rMATS/human/mutant_vs_correction_cortical_salmon/JC.raw.input.A3SS.txt'
        outdir = dirname(cf) %>% str_c(gsub("raw.input.|.txt", "",basename(cf)), sep = '/')

        if (!dir.exists(outdir)) {
            dir.create(outdir,recursive=TRUE)
        }
        cmd = c(cmd,str_c(cf,outdir,0.05,1,sep = ' ') )
    }
    str_c("echo ", cmd %>% str_c(collapse = ' '), "| xargs -d ' ' -n 4 -P 10 -t ", RMATS_325)
}


perform_rmats_stat <- function(species,comparison){
  cmd_stat <- generate_rmats_stat_cmd(species,comparison)
  cmd_stat %>% system()
}



#####start script
write('#!/bin/bash
trap "exit" INT TERM
trap "kill 0" EXIT',tmp_script)


#write cmd to bash script to run rMATS 4.0.2 and 3.2.5
COMPARISON_TABLE %>% pull(comparison) %>% walk ( function(x){
  x=COMPARISON_TABLE %>% filter(comparison==x)
  
  sample_data <- SAMPLE_DATA %>%
    tibble::rownames_to_column(var = "tmp_row_names") %>%
    mutate(!!x$condition_name:= factor(!!parse_expr(x$condition_name))) %>%
    filter(!!parse_expr(x$filter)) %>%
    tibble::column_to_rownames(var = "tmp_row_names")

  cmd = sample_data %>% generate_rmats_count_cmd(species, x$comparison)
  write(cmd,tmp_script,append = T)
})

write("wait",tmp_script,append = T)

cmd = str_c('bash ',tmp_script,sep = '')
system(cmd)

#when finish, parse the result