#!/bin/bash

NUM_TOTAL_THREADS={{cookiecutter.number_total_threads}}

MAIN_DIR=${HOME}/{{cookiecutter.projects_base}}/{{cookiecutter.project_name}}
DATA_DIR=${MAIN_DIR}/data
RESULTS_DIR=${MAIN_DIR}/results
LOG_DIR=${RESULTS_DIR}/logs

SPECIES=()
ENSEMBL_DIR=()
GTF_FILE=()

{% for s in cookiecutter.species.split(' ') %}
SPECIES+=({{ s }})
ENSEMBL_DIR+=(${DATA_DIR}/{{ s }}_ensembl_{{cookiecutter.ensembl_version}})
GTF_FILE+=(${DATA_DIR}/{{ s }}_ensembl_{{cookiecutter.ensembl_version}}/{{cookiecutter.gtf_files[s]}})
{% endfor %}

# make sure we are in project ve https://github.com/sidbdri/cookiecutter-de_analysis_skeleton/issues/108
workon {{cookiecutter.project_name}}

echo "Running get_gene_lengths for species ...."
for species in ${!SPECIES[@]}; do
    if [ ! -f "${ENSEMBL_DIR[$species]}/gene_lengths.csv" ]; then
        mkdir -p ${LOG_DIR}/get_gene_lengths
        echo ${GTF_FILE[$species]}
        get_gene_lengths ${GTF_FILE[$species]} > ${ENSEMBL_DIR[$species]}/gene_lengths.csv 2>${LOG_DIR}/get_gene_lengths/${SPECIES[$species]}.log &
    fi

    if [ ! -f "${ENSEMBL_DIR[$species]}/tx2gene.tsv" ]; then
        # Construct transcript->gene mapping file for tximport
        awk '$3=="transcript" {print $14, $10}' ${GTF_FILE[$species]} | sed 's/"//g;s/;//g' > ${ENSEMBL_DIR[$species]}/tx2gene.tsv &
    fi
done

{% if cookiecutter.sargasso == "yes" %}
python3 -m snakemake -s Snakefile.multispecies_analysis bams -j $NUM_TOTAL_THREADS
python3 -m snakemake -s Snakefile.multispecies_analysis multiqc -j $NUM_TOTAL_THREADS
{% else %}
python3 -m snakemake -s Snakefile.singlespecies_analysis bams -j $NUM_TOTAL_THREADS
python3 -m snakemake -s Snakefile.singlespecies_analysis multiqc -j $NUM_TOTAL_THREADS
{% endif %}

{% if cookiecutter.qSVA == "yes" %}
python3 -m snakemake -s Snakefile.common all_qsva
{% endif %}

#### we check if all the sample are the same strandness settings
#### https://github.com/sidbdri/cookiecutter-de_analysis_skeleton/issues/127
strandedness=`head -1 strand.txt`
strandness_qc=`cat ${MAIN_DIR}/multiqc_data/multiqc_picard_RnaSeqMetrics.txt | tail -n +2 | \
                awk -v s=$strandedness -F '\t' '{ if (s==0 && ($19<=45 || $19>=55))
                                                    print $1"\t"$19;
                                                  else if (s==1 && $19>=5)
                                                    print $1"\t"$19;
                                                  else if (s==2 && $19<=95)
                                                    print $1"\t"$19};'`
if [ ! `echo -n "${strandness_qc}" | wc -l` = 0 ]; then
    echo "Error: The following sample may have the wrong strandedness setting:"
    echo "${strandness_qc}"
    exit
fi


exit;

{% for s in cookiecutter.species.split(' ') %}
Rscript diff_expr_{{ s }}.R
{% endfor %}
