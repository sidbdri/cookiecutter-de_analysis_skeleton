#!/bin/bash

function cleanup {
#   echo "Killing all sub-processes..."
   kill -- -$$
}

trap exit INT
trap cleanup EXIT

export WORKON_HOME=${HOME}/{{cookiecutter.virtualenv_home}}
export PROJECT_HOME=${HOME}/{{cookiecutter.projects_base}}
export VIRTUALENVWRAPPER_PYTHON=/usr/bin/python3.8
source /usr/local/bin/virtualenvwrapper.sh

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
[[ ! -s strand.txt ]] && echo "Error: strand.txt not found. Cannot perform strandness check." && exit
[[ ! -s multiqc_data/multiqc_picard_RnaSeqMetrics.txt ]] && echo "Error: multiqc_picard_RnaSeqMetrics.txt not found. Cannot perform strandness check." && exit

strandedness=`head -1 strand.txt`
strandness_qc=`cat multiqc_data/multiqc_picard_RnaSeqMetrics.txt | \
               awk -F '\t' 'NR==1 { for (i=1; i<=NF; i++) { f[$i] = i}}
                            NR>1  { print $(f["Sample"])"\t"$(f["PCT_R2_TRANSCRIPT_STRAND_READS"])}' | \
               awk -v s=$strandedness -F '\t' '{ if (s==0 && ($2<=45 || $2>=55))
                                                    print $1"\t"$2;
                                                 else if (s==1 && $2>=5)
                                                    print $1"\t"$2;
                                                 else if (s==2 && $2<=95)
                                                    print $1"\t"$2};'`
if [ ! `echo -n "${strandness_qc}" | wc -l` = 0 ]; then
    echo "Error: The following sample may have the wrong strandedness setting:"
    echo "${strandness_qc}"
    exit
fi

exit;

{% for s in cookiecutter.species.split(' ') %}
Rscript diff_expr_{{ s }}.R
{% endfor %}

{% for s in cookiecutter.species.split(' ') %}
#Rscript wgcna_{{ s }}.R
{% endfor %}

# Generate Shiny app
echo Generating Shiny app
bash ./generate_shiny.sh
