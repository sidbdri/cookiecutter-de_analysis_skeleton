#!/bin/bash

set -o nounset
set -o errexit
set -o xtrace

source functions.sh

MAIN_DIR={{cookiecutter.projects_base}}/{{cookiecutter.project_name}}
DATA_DIR=${MAIN_DIR}/data
RNASEQ_DIR=${DATA_DIR}/rnaseq
ENSEMBL_DIR=${DATA_DIR}/{{cookiecutter.species}}_ensembl_{{cookiecutter.ensembl_version}}
STAR_INDEX=${ENSEMBL_DIR}/{{cookiecutter.assembly_name}}
GTF_FILE=${ENSEMBL_DIR}/{{cookiecutter.gtf_file}}
RESULTS_DIR=${MAIN_DIR}/results

QC_DIR=${RESULTS_DIR}/fastqc
MAPPING_DIR=${RESULTS_DIR}/mapped_reads

NUM_THREADS=16

SAMPLES="{{cookiecutter.rnaseq_samples}}"

# Perform QC on raw reads
mkdir -p ${QC_DIR}

for sample in ${SAMPLES}; do
    output_dir=${QC_DIR}/${sample}
    mkdir -p ${output_dir}

    (zcat ${RNASEQ_DIR}/${sample}/*.sanfastq.gz | fastqc --outdir=${output_dir} stdin) &
done

wait

# Map reads
mkdir -p ${MAPPING_DIR}

for sample in ${SAMPLES}; do
    map_reads ${sample} ${STAR_INDEX} ${NUM_THREADS} $(listFiles , ${RNASEQ_DIR}/${sample}/*_1.sanfastq.gz) $(listFiles , ${RNASEQ_DIR}/${sample}/*_2.sanfastq.gz) ${MAPPING_DIR}
done

# Gather all QC data
multiqc -d -f -m featureCounts -m star -m fastqc results
