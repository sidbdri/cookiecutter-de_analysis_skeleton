source('meta_data.R')

SPECIES="{{cookiecutter.species}}"

RMATS_402 = '/opt/rMATS.4.0.2'
RMATS_325 = '/opt/rMATS.3.2.5/MATS/rMATS_Paired.sh'

rMAT_PARA_t = 'paired'
#This is for rMATS4 statistic. we are not using this atm
rMAT_PARA_tstat=1
rMAT_PARA_readLength=75
rMAT_PARA_cstat = 0.1


#This is the number of thread used in each rMATS run
rMAT_PARA_nthread={{cookiecutter.number_total_threads}}
#This is the number of rMATS running in parallal 
XARG_PARA_nthread=1
rMAT3_PARA_splicing_difference = 0.1

#averge count must above a threadhold
avg_count_cutoff=5

sSUMMARY_TB %<>% mutate(Up_regulated_gene=integer(),
                      Down_regulated_gene=integer(),
                      D.E.total_gene=integer())

tmp_script = str_c("rMATS_",SPECIES,".sh",sep = '')

output_folder = 'results/differential_expression/de_rmats'
if (!dir.exists(output_folder)) dir.create(output_folder,recursive=TRUE)

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
  
  reps<-sample_data %>% tibble::rownames_to_column(var = "tmp_row_names") %>%
    group_by(!!parse_expr(x$condition_name)) %>% 
    mutate(bam_file=str_c("results/final_bams/",tmp_row_names,".",species,".bam",sep = '') %>% normalizePath) %>%
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

  cmd = str_c( cmd_count," && ",cmd_stat, sep ="")

  cmd
}

generate_rmats_stat_cmd <- function(species,comparison){
    top_dir<-str_c("results/rMATS/",species,"/",comparison) %>% normalizePath()
    countFiles<-top_dir %>% str_c(c("JC.raw.input.A3SS.txt","JC.raw.input.MXE.txt",
                                    "JC.raw.input.SE.txt","JC.raw.input.A5SS.txt",
                                    "JC.raw.input.RI.txt"), sep = '/')
    cmd=c()
    for (cf in countFiles) {
        #cf = '/home/xinhe/Projects/als_test/results/rMATS/human/mutant_vs_correction_cortical_salmon/JC.raw.input.A3SS.txt'
        outdir = dirname(cf) %>% str_c(gsub("raw.input.|.txt", "",basename(cf)), sep = '/')

        if (!dir.exists(outdir)) {
            dir.create(outdir,recursive=TRUE)
        }
        cmd = c(cmd,str_c(cf,outdir,rMAT3_PARA_splicing_difference,rMAT_PARA_nthread,sep = ' ') )
    }
    str_c("echo ", cmd %>% str_c(collapse = ' '), "| xargs -d ' ' -n 4 -P " , XARG_PARA_nthread, " -t ", RMATS_325)
}


perform_rmats_stat <- function(species,comparison){
  cmd_stat <- generate_rmats_stat_cmd(species,comparison)
  cmd_stat %>% system()
}



#####start script
write('#!/bin/bash
#trap "exit" INT TERM
#trap "kill 0" EXIT',tmp_script)


#write cmd to bash script to run rMATS 4.0.2 and 3.2.5
COMPARISON_TABLE %>% pull(comparison) %>% walk ( function(x){
  x=COMPARISON_TABLE %>% filter(comparison==x)
  
  sample_data <- SAMPLE_DATA %>%
    tibble::rownames_to_column(var = "tmp_row_names") %>%
    mutate(!!x$condition_name:= factor(!!parse_expr(x$condition_name))) %>%
    filter(!!parse_expr(x$filter)) %>%
    tibble::column_to_rownames(var = "tmp_row_names")

  cmd = sample_data %>% generate_rmats_count_cmd(SPECIES, x$comparison)
  write(cmd,tmp_script,append = T)
})


cmd = str_c('bash ',tmp_script,sep = '')
system(cmd)

#when finish, parse the result

####################################
#when finish, parse the result

total_dds_data <- get_total_dds(SAMPLE_DATA,species=SPECIES)

gene_info <- get_gene_info(SPECIES)
gene_lengths <- read_csv(str_c("data/",SPECIES,"_ensembl_",{{cookiecutter.ensembl_version}},"/gene_lengths.csv", sep=""))

results <- total_dds_data %>% get_count_data()

# fpkms %<>% mutate(
#   P10_Ctx_KO_fpkm_avg = !!get_avg_fpkm(filter="age=='P10' & genotype=='KO' & region=='Ctx'"),
#   P10_Piri_KO_fpkm_avg = !!get_avg_fpkm(filter="age=='P10' & genotype=='KO' & region=='Piri'")
# # )

fpkms <- results %>%
get_fpkms(gene_lengths, colnames(results) %>% tail(-1), "_fpkm")


results %<>%
inner_join(fpkms) %>%
    inner_join(gene_info) %>%
    inner_join(gene_lengths)

# we do not need counts for samples
results %<>% dplyr::select(-one_of(SAMPLE_DATA %>% rownames()))


#load rMATS results and join with result table
COMPARISON_TABLE %>% pull(comparison) %>% walk ( function(x){
    x=COMPARISON_TABLE %>% filter(comparison==x)

    sample_data <- SAMPLE_DATA %>%
        tibble::rownames_to_column(var = "tmp_row_names") %>%
        mutate(!!x$condition_name:= factor(!!parse_expr(x$condition_name))) %>%
        filter(!!parse_expr(x$filter)) %>%
        tibble::column_to_rownames(var = "tmp_row_names")

    top_dir = str_c("results/rMATS/",SPECIES,"/",x$comparison)
    out_dir = str_c(output_folder,SPECIES,x$comparison,sep = '/')

    if (!dir.exists(out_dir)) {
        dir.create(out_dir,recursive=TRUE)
    }

    for (f in top_dir %>% list.dirs() %>% tail(-1) ) {
        # f='results/rMATS/human/mutant_vs_correction_motor/JC.SE'
        as.event = basename(f) %>% strsplit('.',fixed = T) %>% unlist() %>% extract(2)

        result.table <- read_tsv(str_c(f,'rMATS_Result.txt',sep = '/'))
      
        #replace ',' with ';' to avoid excel problem in the result
        result.table %<>% 
          mutate_at(.vars = vars("IJC_SAMPLE_1","SJC_SAMPLE_1","IJC_SAMPLE_2","SJC_SAMPLE_2","IncLevel1","IncLevel2"),
                    .funs = funs(gsub(",",";",.)))

        #workout avg count
        result.table %<>%  dplyr::select(ID,IJC_SAMPLE_1,SJC_SAMPLE_1,IJC_SAMPLE_2,SJC_SAMPLE_2) %>%
            tidyr::unite(tmp,c("IJC_SAMPLE_1","SJC_SAMPLE_1","IJC_SAMPLE_2","SJC_SAMPLE_2"),sep=';') %>%
            mutate(avg_count=strsplit(tmp,';')%>% lapply(as.numeric) %>%lapply(mean) %>% unlist()) %>%
            dplyr::select(-tmp) %>% left_join(result.table)


        #join gene info
        result.table <- read_tsv(str_c(f,'/..','/fromGTF.',as.event,'.txt')) %>%
          inner_join(result.table)
        
        #workout pre comparison avg fpkm
        avg<- fpkms %>% dplyr::select(one_of(str_c(sample_data %>% rownames(),'_fpkm'))) %>% 
          mutate(avg_fpkm=rowMeans(.)) %>%  
          dplyr::pull(avg_fpkm)

        result.table <-  get("results",envir=.GlobalEnv) %>% 
          dplyr::select(gene, gene_name, chromosome, description, entrez_id, 
                        gene_type, gene_length, max_transcript_length,
                        one_of(str_c(sample_data %>% rownames(),'_fpkm'))) %>%
          mutate(avg_fpkm=avg) %>% 
          inner_join(result.table,by=c('gene' = 'GeneID')) %>%
          dplyr::select(everything(),-geneSymbol,-chr)

        #filter low counts
        result.table %<>% filter(avg_count>avg_count_cutoff) %>% dplyr::select(-avg_count)

        #fill summary table
        SUMMARY_TB <- get("SUMMARY_TB", envir = .GlobalEnv) %>% 
          add_row(Comparison = x$comparison%>% str_c(basename(f),sep = '_'), DESeq_model_formula = design(get("total_dds_data", envir = .GlobalEnv)) %>% format(),
                  Condition_tested = x$condition_name ,
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
          )
        assign("SUMMARY_TB", SUMMARY_TB,envir = .GlobalEnv)

        #save results
        result.table %>% 
          dplyr::select(-PValue,-FDR, everything(), 
                        -one_of(str_c(sample_data %>% rownames(),'_fpkm')),
                        -IncFormLen,-SkipFormLen) %>% 
          write_csv(str_c(out_dir,'/AS_results_fpkm_',basename(f),'_',SPECIES,".csv"))
    }
})

# modify/reformat summary table
SUMMARY_TB %>% dplyr::select(-DESeq_model_formula) %>% 
  write_csv(str_c(output_folder,"/AS_summary_",SPECIES,".csv"))

