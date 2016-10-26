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
SALMON_INDEX=${ENSEMBL_DIR}/{{cookiecutter.salmon_index}}
KALLISTO_INDEX=${ENSEMBL_DIR}/{{cookiecutter.kallisto_index}}
RESULTS_DIR=${MAIN_DIR}/results

QC_DIR=${RESULTS_DIR}/fastqc
MAPPING_DIR=${RESULTS_DIR}/mapped_reads
COUNTS_DIR=${RESULTS_DIR}/read_counts
SALMON_QUANT_DIR=${RESULTS_DIR}/salmon_quant
KALLISTO_QUANT_DIR=${RESULTS_DIR}/kallisto_quant
DIFF_EXPR_DIR=${RESULTS_DIR}/differential_expression

NUM_THREADS=16

SAMPLES="{{cookiecutter.rnaseq_samples}}"

##### Perform QC on raw reads
mkdir -p ${QC_DIR}

for sample in ${SAMPLES}; do
    output_dir=${QC_DIR}/${sample}
    mkdir -p ${output_dir}

    (zcat ${RNASEQ_DIR}/${sample}/*.{{cookiecutter.fastq_suffix}} | fastqc --outdir=${output_dir} stdin) &
done

wait

##### Map reads
mkdir -p ${MAPPING_DIR}

for sample in ${SAMPLES}; do
    map_reads ${sample} ${STAR_INDEX} ${NUM_THREADS} $(listFiles , ${RNASEQ_DIR}/${sample}/*_1.{{cookiecutter.fastq_suffix}}) $(listFiles , ${RNASEQ_DIR}/${sample}/*_2.{{cookiecutter.fastq_suffix}}) ${MAPPING_DIR}
done

##### Count mapped reads
mkdir -p ${COUNTS_DIR}

for sample in ${SAMPLES}; do
    count_reads_for_features ${NUM_THREADS} ${GTF_FILE} ${MAPPING_DIR}/${sample}.bam ${COUNTS_DIR}/${sample}.counts
done

##### Quantify transcript expression with Salmon
mkdir -p ${SALMON_QUANT_DIR}

for sample in ${SAMPLES}; do
    salmon quant -i ${SALMON_INDEX} -l A -1 $(listFiles , ${RNASEQ_DIR}/${sample}/*_1.{{cookiecutter.fastq_suffix}}) -2 $(listFiles , ${RNASEQ_DIR}/${sample}/*_2.{{cookiecutter.fastq_suffix}}) -o ${SALMON_QUANT_DIR}/${sample} --seqBias --gcBias -p ${NUM_THREADS} -g ${GTF_FILE}
done

##### Quantify transcript expression with Kallisto
mkdir -p ${KALLISTO_QUANT_DIR}

for sample in ${SAMPLES}; do
    kallisto quant -i ${KALLISTO_INDEX} -o ${KALLISTO_QUANT_DIR}/${sample} --bias --rf-stranded -t ${NUM_THREADS} $(ls -1 ${RNASEQ_DIR}/${sample}/*_1.{{cookiecutter.fastq_suffix}} | sed -r 's/(.*)/\1 \1/' | sed -r 's/(.* .*)_1.{{cookiecutter.fastq_suffix}}/\1_2.{{cookiecutter.fastq_suffix}}/' | tr '\n' ' ')
done

##### Gather all QC data
multiqc -d -f -m featureCounts -m star -m fastqc -m salmon -m kallisto results

##### Perform differential expression
mkdir -p ${DIFF_EXPR_DIR}

get_gene_lengths ${GTF_FILE} > ${RESULTS_DIR}/gene_lengths.csv

# Construct transcript->gene mapping file for tximport
awk '$3=="transcript" {print $14, $10}' ${GTF_FILE} | sed 's/"//g;s/;//g' > ${RESULTS_DIR}/tx2gene.tsv

Rscript diff_expr.R

clean_de_results ${DIFF_EXPR_DIR}/deseq2_results.csv
clean_de_results ${DIFF_EXPR_DIR}/deseq2_results_fpkm.csv
