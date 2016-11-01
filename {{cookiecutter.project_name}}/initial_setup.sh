#!/bin/bash

export WORKON_HOME={{cookiecutter.virtualenv_home}}
export PROJECT_HOME={{cookiecutter.projects_base}}
source /usr/local/bin/virtualenvwrapper.sh

mkproject -f {{cookiecutter.project_name}}

pip install multiqc
pip install git+https://github.com/lweasel/transcript_utils.git

RNASEQ_DIR=data/rnaseq
GENOME_DATA_DIR=/srv/data/genome/{{cookiecutter.species}}/ensembl-{{cookiecutter.ensembl_version}}
ENSEMBL_DIR=data/{{cookiecutter.species}}_ensembl_{{cookiecutter.ensembl_version}}

mkdir -p ${RNASEQ_DIR}

for sample in {{cookiecutter.rnaseq_samples}}; do
    ln -s {{cookiecutter.rnaseq_samples_dir}}/$sample ${RNASEQ_DIR}/$sample
done 

mkdir -p ${ENSEMBL_DIR}

ln -s ${GENOME_DATA_DIR}/STAR_indices/{{cookiecutter.assembly_name}} ${ENSEMBL_DIR}
ln -s ${GENOME_DATA_DIR}/{{cookiecutter.gtf_file}} ${ENSEMBL_DIR}
ln -s ${GENOME_DATA_DIR}/genes.tsv ${ENSEMBL_DIR}
ln -s ${GENOME_DATA_DIR}/{{cookiecutter.salmon_index}} ${ENSEMBL_DIR}
ln -s ${GENOME_DATA_DIR}/{{cookiecutter.kallisto_index}} ${ENSEMBL_DIR}

git init
mv gitignore .gitignore
