#!/bin/bash

export WORKON_HOME=${HOME}/{{cookiecutter.virtualenv_home}}
export PROJECT_HOME=${HOME}/{{cookiecutter.projects_base}}
source /usr/local/bin/virtualenvwrapper.sh

mkproject -f {{cookiecutter.project_name}}

# If this is the first time the script has been run, We find out the current
# master hash of transcript-utils and sidbdri-utils and put them in the
# config.sh, if they are not set. We also set up the directory as a git
# repository.
function findHashFromBranchName {
    org=$1
    repo_name=$2
    commit_hash=$3
    echo $(curl --silent -H "Accept: application/vnd.github.VERSION.sha" \
    https://api.github.com/repos/${org}/${repo_name}/commits/${commit_hash})
}

if grep -Fq "unknown_hash" config.sh; then
    echo "# We replace the unknown hash with the master branch hash"
    sed -i "s/unknown_hash_transcript-utils/$(findHashFromBranchName "sidbdri" "transcript-utils" "master")/" config.sh
    sed -i "s/unknown_hash_sidbdri-utils/$(findHashFromBranchName "sidbdri" "sidbdri-utils" "master")/" config.sh
    {% if cookiecutter.sargasso == "yes" %}
    sed -i "s/unknown_hash_sargasso/$(findHashFromBranchName "statbio" "Sargasso" "master")/" config.sh
    {% endif %}

    git init
    mv gitignore .gitignore
fi

source config.sh

pip install scipy
pip install multiqc
pip install git+https://github.com/sidbdri/transcript-utils.git@${transcript_utils_hash}
{% if cookiecutter.sargasso == "yes" %}
pip install git+https://github.com/statbio/Sargasso.git@${sargasso_hash}
{% endif %}

## Clone sidbdri-utils package
rm -rf sidbdri-utils && git clone https://github.com/sidbdri/sidbdri-utils.git
(cd sidbdri-utils && git checkout -q ${sidbdri_utils_hash})
source sidbdri-utils/includes.sh

DATA_DIR=data
RNASEQ_DIR=${DATA_DIR}/rnaseq
PICARD_DATA=${DATA_DIR}/picard
PICARD=/opt/picard-tools-{{cookiecutter.picard_version}}/picard.jar

mkdir -p ${RNASEQ_DIR}
for sample in {{cookiecutter.rnaseq_samples}}; do
    ln -s {{cookiecutter.rnaseq_samples_dir}}/$sample ${RNASEQ_DIR}/$sample
done

{% for s in cookiecutter.species.split(' ') %}
ENSEMBL_DIR=${DATA_DIR}/{{ s }}_ensembl_{{cookiecutter.ensembl_version}}
GENOME_DATA_DIR=/srv/data/genome/{{ s }}/ensembl-{{cookiecutter.ensembl_version}}
REF_FLAT=${PICARD_DATA}/{{ s }}/{{cookiecutter.rff_files[s]}}
GTF_FILE=${GENOME_DATA_DIR}/{{cookiecutter.gtf_files[s]}}

mkdir -p ${ENSEMBL_DIR}
ln -s ${GENOME_DATA_DIR}/STAR_indices ${ENSEMBL_DIR}
ln -s ${GENOME_DATA_DIR}/{{cookiecutter.gtf_files[s]}} ${ENSEMBL_DIR}
ln -s ${GENOME_DATA_DIR}/SALMON_indices ${ENSEMBL_DIR}
ln -s ${GENOME_DATA_DIR}/KALLISTO_indices ${ENSEMBL_DIR}
ln -s ${GENOME_DATA_DIR}/*orthologs.tsv  ${ENSEMBL_DIR}

## Instead of using the genes.tsv from ${GENOME_DATA_DIR}, we download it from Ensembl
download_gene_tb {{ s }} {{cookiecutter.ensembl_version}} > ${ENSEMBL_DIR}/genes.tsv

# Generating refFlat file for Picard RNA-seq metrics
generate_picard_refFlat ${PICARD_DATA}/{{ s }} {{ s }} {{cookiecutter.ensembl_version}} ${GTF_FILE} &
{% endfor %}

{% if "human" not in cookiecutter.species.split(' ') %}
HUMAN_ENSEMBL_DIR=${DATA_DIR}/human_ensembl_{{cookiecutter.ensembl_version}}
mkdir -p ${HUMAN_ENSEMBL_DIR}

## Instead of using the genes.tsv from ${GENOME_DATA_DIR}, we download it from Ensembl
download_gene_tb human {{cookiecutter.ensembl_version}} > ${HUMAN_ENSEMBL_DIR}/genes.tsv
{% endif %}

ln -s /srv/data/genome/human/msigdb ${DATA_DIR}
