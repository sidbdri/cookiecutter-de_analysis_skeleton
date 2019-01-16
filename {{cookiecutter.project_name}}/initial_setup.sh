#!/bin/bash

export WORKON_HOME=${HOME}/{{cookiecutter.virtualenv_home}}
export PROJECT_HOME=${HOME}/{{cookiecutter.projects_base}}
source /usr/local/bin/virtualenvwrapper.sh

mkproject -f {{cookiecutter.project_name}}

pip install scipy
pip install multiqc
pip install git+https://github.com/sidbdri/transcript-utils.git
{% if cookiecutter.sargasso == "yes" %}
pip install git+https://github.com/statbio/Sargasso.git@master
{% endif %}

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
ln -s ${GENOME_DATA_DIR}/genes.tsv ${ENSEMBL_DIR}
ln -s ${GENOME_DATA_DIR}/SALMON_indices ${ENSEMBL_DIR}
ln -s ${GENOME_DATA_DIR}/KALLISTO_indices ${ENSEMBL_DIR}
ln -s ${GENOME_DATA_DIR}/*orthologs.tsv  ${ENSEMBL_DIR}

# Generating refFlat file for Picard RNA-seq metrics
mkdir -p ${PICARD_DATA}/{{ s }}
gtfToGenePred -genePredExt -geneNameAsName2 ${GTF_FILE} ${PICARD_DATA}/{{ s }}/refFlat.tmp.txt
paste <(cut -f 12 ${PICARD_DATA}/{{ s }}/refFlat.tmp.txt) <(cut -f 1-10 ${PICARD_DATA}/{{ s }}/refFlat.tmp.txt) > ${REF_FLAT}
rm ${PICARD_DATA}/{{ s }}/refFlat.tmp.txt

#{% if s!="human" %}
#ln -s ${GENOME_DATA_DIR}/human_orthologs.tsv ${ENSEMBL_DIR}
#{% endif %}
#{% endfor %}

{% if "human" not in cookiecutter.species.split(' ') %}
HUMAN_ENSEMBL_DIR=${DATA_DIR}/human_ensembl_{{cookiecutter.ensembl_version}}
mkdir -p ${HUMAN_ENSEMBL_DIR}
ln -s /srv/data/genome/human/ensembl-{{cookiecutter.ensembl_version}}/genes.tsv ${HUMAN_ENSEMBL_DIR}
{% endif %}

ln -s /srv/data/genome/human/msigdb ${DATA_DIR}

git init
mv gitignore .gitignore
