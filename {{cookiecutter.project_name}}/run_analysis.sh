#!/bin/bash

set -o nounset
set -o errexit
set -o xtrace

source functions.sh

MAIN_DIR={{cookiecutter.projects_base}}/{{cookiecutter.project_name}}
DATA_DIR=${MAIN_DIR}/data
RNASEQ_DIR=${DATA_DIR}/rnaseq
RESULTS_DIR=${MAIN_DIR}/results

QC_DIR=${RESULTS_DIR}/fastqc

NUM_THREADS=16

SAMPLES="{{cookiecutter.rnaseq_samples}}"

# Perform QC on raw reads
mkdir -p ${QC_DIR}

for sample in ${SAMPLES}; do
    output_dir=${QC_DIR}/${sample}
    mkdir -p ${output_dir}

    (zcat ${RNASEQ_DIR}/${sample}/*.sanfastq.gz | fastqc --outdir=${output_dir} stdin &)
done

# Gather all QC data
multiqc -d -f -m featureCounts -m star -m fastqc results
