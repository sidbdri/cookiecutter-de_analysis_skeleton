#!/bin/bash

set -o nounset
set -o errexit
set -o xtrace

source functions.sh

MAIN_DIR=${HOME}/{{cookiecutter.projects_base}}/{{cookiecutter.project_name}}
DATA_DIR=${MAIN_DIR}/data
PICARD_DATA=${DATA_DIR}/picard
RNASEQ_DIR=${DATA_DIR}/rnaseq
ENSEMBL_DIR=${DATA_DIR}/{{cookiecutter.species}}_ensembl_{{cookiecutter.ensembl_version}}
STAR_INDEX=${ENSEMBL_DIR}/{{cookiecutter.assembly_names[cookiecutter.species]}}
GTF_FILE=${ENSEMBL_DIR}/{{cookiecutter.gtf_files[cookiecutter.species]}}
REF_FLAT=${PICARD_DATA}/{{cookiecutter.rff_files[cookiecutter.species]}}
REFERENCE=${ENSEMBL_DIR}/{{cookiecutter.species}}_{{cookiecutter.assembly_names[cookiecutter.species]}}.fa
SALMON_INDEX=${ENSEMBL_DIR}/{{cookiecutter.salmon_index}}
KALLISTO_INDEX=${ENSEMBL_DIR}/{{cookiecutter.kallisto_index}}
RESULTS_DIR=${MAIN_DIR}/results

QC_DIR=${RESULTS_DIR}/fastqc
MAPPING_DIR=${RESULTS_DIR}/mapped_reads
PICARD_DIR=${RESULTS_DIR}/alignment_metrics
COUNTS_DIR=${RESULTS_DIR}/read_counts
SALMON_QUANT_DIR=${RESULTS_DIR}/salmon_quant
KALLISTO_QUANT_DIR=${RESULTS_DIR}/kallisto_quant
DIFF_EXPR_DIR=${RESULTS_DIR}/differential_expression

NUM_THREADS=16

SAMPLES="{{cookiecutter.rnaseq_samples}}"
PAIRED_END_READ="{{cookiecutter.paired_end_read}}"

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
    if [ $PAIRED_END_READ = "yes" ]; then
        map_reads ${sample} ${STAR_INDEX} ${NUM_THREADS} $(listFiles , ${RNASEQ_DIR}/${sample}/*{{cookiecutter.read1_identifier}}.{{cookiecutter.fastq_suffix}}) $(listFiles , ${RNASEQ_DIR}/${sample}/*{{cookiecutter.read2_identifier}}.{{cookiecutter.fastq_suffix}}) ${MAPPING_DIR}
    else
        map_reads ${sample} ${STAR_INDEX} ${NUM_THREADS} $(listFiles , ${RNASEQ_DIR}/${sample}/*.{{cookiecutter.fastq_suffix}}) "" ${MAPPING_DIR}
    fi
done

##### Run Picard alignment metrics summary
mkdir -p ${PICARD_DIR}

grep rRNA ${GTF_FILE} | cut -s -f 1,4,5,7,9 > ${PICARD_DATA}/intervalListBody.txt

for sample in ${SAMPLES}; do
    sambamba view -H ${MAPPING_DIR}/${sample}.sorted.bam > ${PICARD_DATA}/${sample}_header.txt
    cat ${PICARD_DATA}/${sample}_header.txt ${PICARD_DATA}/intervalListBody.txt > ${PICARD_DATA}/${sample}.txt

    picard_rnaseq_metrics ${sample} ${MAPPING_DIR} ${PICARD_DIR} ${REF_FLAT} ${PICARD_DATA}
done

##### Count mapped reads
mkdir -p ${COUNTS_DIR}

first_sample="TRUE"

for sample in ${SAMPLES}; do
    [[ "${first_sample}" == "FALSE" ]] || {
        first_sample="FALSE"
        count_reads_for_features_strand_test ${NUM_THREADS} ${GTF_FILE} ${MAPPING_DIR}/${sample}.bam ${COUNTS_DIR}/strand_test.${sample}.counts
    }

    count_reads_for_features ${NUM_THREADS} ${GTF_FILE} ${MAPPING_DIR}/${sample}.bam ${COUNTS_DIR}/${sample}.counts
done

##### Quantify transcript expression with Salmon
mkdir -p ${SALMON_QUANT_DIR}


for sample in ${SAMPLES}; do
    if [ $PAIRED_END_READ = "yes" ]; then
        salmon quant -i ${SALMON_INDEX} -l A -1 $(listFiles ' ' ${RNASEQ_DIR}/${sample}/*{{cookiecutter.read1_identifier}}.{{cookiecutter.fastq_suffix}}) -2 $(listFiles ' ' ${RNASEQ_DIR}/${sample}/*{{cookiecutter.read2_identifier}}.{{cookiecutter.fastq_suffix}}) -o ${SALMON_QUANT_DIR}/${sample} --seqBias --gcBias -p ${NUM_THREADS} -g ${GTF_FILE}
    else
        salmon quant -i ${SALMON_INDEX} -l A -r $(listFiles ' ' ${RNASEQ_DIR}/${sample}/*.{{cookiecutter.fastq_suffix}}) -o ${SALMON_QUANT_DIR}/${sample} --seqBias --gcBias -p ${NUM_THREADS} -g ${GTF_FILE}
   fi
done

##### Quantify transcript expression with Kallisto
mkdir -p ${KALLISTO_QUANT_DIR}

for sample in ${SAMPLES}; do
    if [ $PAIRED_END_READ = "yes" ]; then
        kallisto quant -i ${KALLISTO_INDEX} -o ${KALLISTO_QUANT_DIR}/${sample} --bias --rf-stranded -t ${NUM_THREADS} $(ls -1 ${RNASEQ_DIR}/${sample}/*{{cookiecutter.read1_identifier}}.{{cookiecutter.fastq_suffix}} | sed -r 's/(.*)/\1 \1/' | sed -r 's/(.* .*){{cookiecutter.read1_identifier}}.{{cookiecutter.fastq_suffix}}/\1{{cookiecutter.read2_identifier}}.{{cookiecutter.fastq_suffix}}/' | tr '\n' ' ')
    else
        kallisto quant -i ${KALLISTO_INDEX} -o ${KALLISTO_QUANT_DIR}/${sample} --bias --single -t ${NUM_THREADS} $(ls -1 ${RNASEQ_DIR}/${sample}/*.{{cookiecutter.fastq_suffix}} )
    fi
done

##### Gather all QC data
multiqc -d -f -m featureCounts -m star -m fastqc -m salmon -m kallisto -m picard -x differential_expression results

##### Perform differential expression
mkdir -p ${DIFF_EXPR_DIR}/go
mkdir -p ${DIFF_EXPR_DIR}/graphs
mkdir -p ${DIFF_EXPR_DIR}/gsa

get_gene_lengths <(tail -n +6 ${GTF_FILE}) > ${RESULTS_DIR}/gene_lengths.csv

# Construct transcript->gene mapping file for tximport
awk '$3=="transcript" {print $14, $10}' ${GTF_FILE} | sed 's/"//g;s/;//g' > ${RESULTS_DIR}/tx2gene.tsv

Rscript diff_expr.R

for de_results in ${DIFF_EXPR_DIR}/*.csv; do
    clean_de_results ${de_results}
done
