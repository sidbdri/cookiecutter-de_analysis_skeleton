#!/usr/bin/env bash

function cleanup {
   echo "Cleaning tmp..."
   rm -rf ${MAIN_DIR}/*.tmp
   echo "Killing all sub-processes..."
   kill -- -$$
}

trap exit INT
trap cleanup EXIT

#set -o nounset
#set -o errexit
#set -o xtrace

## filling in the relevant settings
SARGASSO_TEST_DIR=~/tmp/sargasso_test
RNASEQ_DIR='/srv/data/ghardingham/snakemake_test'
READ1_IDENTIFIER='01_1'
READ2_IDENTIFIER='01_2'
FASTQ_SUFFIX='fastq.gz'
PAIRED_END_READ='no'
SAMPLES='1467_A 1467_V 1468_P'
SAMPLES_ORIGIN='mouse mouse rat'
SPECIES_PARA=()
SPECIES_PARA+=("human /srv/data/genome/human/ensembl-99/STAR_indices/primary_assembly")
SPECIES_PARA+=("mouse /srv/data/genome/mouse/ensembl-99/STAR_indices/primary_assembly")
SPECIES_PARA+=("rat /srv/data/genome/rat/ensembl-99/STAR_indices/toplevel")
MISMATCH_SETTING='0 2 4'
MINMATCH_SETTING='0 2 5'
MUTLIMAP_SETTING='1'

MAPPER_EXECUTABLE=STAR2.7.0f
NUM_TOTAL_THREADS=16
##


function listFilesNoNewLine {
    local DELIMITER=$1
    shift
    local FILES=$@

    local OUTPUT=$(ls -1 $FILES | tr '\n' "${DELIMITER}")
    echo -n ${OUTPUT%$DELIMITER}
}



function addSample2tsv {
    local SAMPLE_TSV=$1 && shift
    local BASE_DIR=$1 && shift
    local READ1_IDENTIFIER=$1 && shift
    local READ2_IDENTIFIER=$1 && shift
    local FASTQ_SUFFIX=$1 && shift
    local PAIRED_READ=$1 && shift
    if [ -z "$READ2_IDENTIFIER" ];then
        local PAIRED_READ=0
    fi
    echo -ne '' > ${SAMPLE_TSV}
    local SAMPLE=$@
    for sample in ${SAMPLE}; do
        echo -ne ${sample}" " >> ${SAMPLE_TSV}
        echo -n $(listFilesNoNewLine "," ${BASE_DIR}/${sample}/*${READ1_IDENTIFIER}.${FASTQ_SUFFIX}) >> ${SAMPLE_TSV}
        if [ "$PAIRED_READ" == "yes" ]; then
            echo -n " "  >> ${SAMPLE_TSV}
            echo -n $(listFilesNoNewLine "," ${BASE_DIR}/${sample}/*${READ2_IDENTIFIER}.${FASTQ_SUFFIX}) >> ${SAMPLE_TSV}
        fi
        echo "" >> ${SAMPLE_TSV}
    done
}




TMP_DIR=${SARGASSO_TEST_DIR}/tmp
mkdir -p  ${TMP_DIR}

SAMPLE_TSV=${SARGASSO_TEST_DIR}/sample.tsv
addSample2tsv ${SAMPLE_TSV} ${RNASEQ_DIR} ${READ1_IDENTIFIER} ${READ2_IDENTIFIER} ${FASTQ_SUFFIX} ${PAIRED_END_READ} ${SAMPLES}



####################################################################################
#### MAPPING
echo "Running Sargasso ...."
for mismatch in ${MISMATCH_SETTING}; do
    for minmatch in ${MINMATCH_SETTING}; do
        for multimap in ${MUTLIMAP_SETTING}; do
            out_dir=${SARGASSO_TEST_DIR}/${mismatch}_${minmatch}_${multimap}
            species_separator rnaseq --mapper-executable ${MAPPER_EXECUTABLE}  --sambamba-sort-tmp-dir=${TMP_DIR} \
                    --mismatch-threshold ${mismatch} --minmatch-threshold ${minmatch} --multimap-threshold ${multimap} --reject-multimaps \
                    --num-threads ${NUM_TOTAL_THREADS} \
                    ${SAMPLE_TSV} ${out_dir} ${SPECIES_PARA[@]}
            cd ${out_dir} && make >sargasso.log 2>&1 &
            wait $(jobs -p)
        done
    done
done




echo "
library(dplyr)
library(ggplot2)
library(magrittr)
library(purrr)
library(readr)
library(reshape2)
library(stringr)

samples=c(\""`echo ${SAMPLES} | sed 's/ /","/g'`"\")
origin=c(\""`echo ${SAMPLES_ORIGIN} | sed 's/ /","/g'`"\") %>% set_names(samples)

tb <- lapply(c("`echo ${MISMATCH_SETTING} | sed 's/ /,/g'`"),function(number_mismatch){
          lapply(c("`echo ${MINMATCH_SETTING} | sed 's/ /,/g'`"),function(min_match){
            lapply(c("`echo ${MUTLIMAP_SETTING} | sed 's/ /,/g'`"),function(number_multimap){
                file=file.path(str_c(number_mismatch,'_',min_match,'_',number_multimap),'filtered_reads','overall_filtering_summary.txt')
                read_csv(file) %>% dplyr::select(Sample,contains('Reads')) %>%
                  mutate(Par=str_c(number_mismatch,min_match,number_multimap,sep = '_')) %>%
                  tidyr::pivot_longer(cols=contains('Reads'),names_to='type',values_to = 'count') %>%
                  mutate(origin=origin[Sample])
            }) %>% purrr::reduce(rbind)
          }) %>% purrr::reduce(rbind)
        }) %>% purrr::reduce(rbind)


count_table <- lapply(c(\""`echo ${SAMPLES} | sed 's/ /","/g'`"\")%>%set_names(.),function(sample){
  file<-list.files(path = '.', pattern = str_c(sample,'."${SAMPLES_ORIGIN%% *}".log.out',sep=''), all.files = FALSE,full.names = TRUE, recursive = TRUE,ignore.case = TRUE)[1]
  readLines(file,n= grep('Number of input reads',readLines(file))) %>% tail(1) %>% strsplit('\t') %>% extract2(1) %>% extract(2)
})

tb %<>% mutate(total_count=count_table[Sample] %>% unlist() %>% as.numeric()) %>%
  mutate(prec=count/total_count) %>%
  mutate_at(vars(one_of(c('Par','type'))),as.factor) %>%
  mutate(Par=factor(Par,level=unique(Par)))

tb %>%
  ggplot(aes(x=Par, y=count)) +
  geom_point(aes(color=Sample)) +
  facet_wrap(~origin + type,scales='free')+ coord_flip()
ggsave('counts.png')

tb %>%
  ggplot(aes(x=Par, y=prec)) +
  geom_point(aes(color=Sample)) +
  facet_wrap(~origin + type,scales='free')+ coord_flip() +
    scale_y_continuous(labels = scales::percent_format(accuracy = 0.01))
ggsave('prec.png')
" > ${SARGASSO_TEST_DIR}/plot.R

Rscript ${SARGASSO_TEST_DIR}/plot.R


