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

strandedness=`head -1 strand.txt`

python strandness_qc.py
# Check if the Python script returned an error
[[ $? -ne 0 ]] && echo "Error: strandness_qc.py failed. Exiting." && exit 1

# NB. remove this exit line before committing project to a GitHub repo
exit;

{% for s in cookiecutter.species.split(' ') %}
Rscript diff_expr_{{ s }}.R | tee ${LOG_DIR}/diff_expr_{{ s }}.log
Rscript diff_expr_tx_{{ s }}.R | tee ${LOG_DIR}/diff_expr_tx_{{ s }}.log
{% endfor %}

{% for s in cookiecutter.species.split(' ') %}
#Rscript rMATS_{{ s }}.R | tee ${LOG_DIR}/rMATS_{{ s }}.log
{% endfor %}

{% for s in cookiecutter.species.split(' ') %}
#Rscript wgcna_{{ s }}.R  | tee ${LOG_DIR}/wgcna_{{ s }}.log
{% endfor %}

#### Create deliverable results in gzipped tar archive
dated_results_name=$(date '+%Y%m%d').{{cookiecutter.project_name}}.results
mkdir -p ${dated_results_name}

if [ -f "${RESULTS_DIR}/sessionInfo.txt" ]; then
        cp ${RESULTS_DIR}/sessionInfo.txt ${dated_results_name}/
elif [ -f "${MAIN_DIR}/sessionInfo.txt" ]; then
        cp ${MAIN_DIR}/sessionInfo.txt ${dated_results_name}/
fi

if [ -f "${RESULTS_DIR}/multiqc_report.html" ]; then
        cp ${RESULTS_DIR}/multiqc_report.html ${dated_results_name}/
elif [ -f "${MAIN_DIR}/multiqc_report.html" ]; then
        cp ${MAIN_DIR}/multiqc_report.html ${dated_results_name}/
fi

cp -r ${RESULTS_DIR}/differential_expression* ${dated_results_name}/
cp -r ${RESULTS_DIR}/read_counts ${dated_results_name}/

find ${dated_results_name} -name "*count*.csv" -exec rm {} \;
find ${dated_results_name} -name "*genes_in_sets*.csv" -exec rm {} \;

tar cfz results/${dated_results_name}.tar.gz -P ${dated_results_name}
rm -r ${dated_results_name}

# Generate Shiny app
echo Generating Shiny app
bash ./generate_shiny.sh

## one might have used/installed R packages in the project that are not in the docker image,
## for example, one installed a package called 'packageA' when developing the project,
## this package will be install under the project directory ./renv/library/,
## but will not automatically be recorded in the renv.lock file.
##
## Thus before commit project to github for a docker run, one needs to make sure that
## any R package that is used in the project is recorded in the renv.lock file.
## This can be done following the below script:
# bash ./update_r_packages.sh