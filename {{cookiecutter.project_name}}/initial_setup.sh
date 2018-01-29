#!/bin/bash

export WORKON_HOME={{cookiecutter.virtualenv_home}}
export PROJECT_HOME={{cookiecutter.projects_base}}
source /usr/local/bin/virtualenvwrapper.sh

mkproject -f {{cookiecutter.project_name}}

pip install multiqc
pip install git+https://github.com/sidbdri/transcript-utils.git

DATA_DIR=data
RNASEQ_DIR=${DATA_DIR}/rnaseq
GENOME_DATA_DIR=/srv/data/genome/{{cookiecutter.species}}/ensembl-{{cookiecutter.ensembl_version}}
ENSEMBL_DIR=${DATA_DIR}/{{cookiecutter.species}}_ensembl_{{cookiecutter.ensembl_version}}

mkdir -p ${RNASEQ_DIR}

for sample in {{cookiecutter.rnaseq_samples}}; do
    ln -s {{cookiecutter.rnaseq_samples_dir}}/$sample ${RNASEQ_DIR}/$sample
done 

mkdir -p ${ENSEMBL_DIR}

ln -s ${GENOME_DATA_DIR}/STAR_indices/{{cookiecutter.assembly_names[cookiecutter.species]}} ${ENSEMBL_DIR}
ln -s ${GENOME_DATA_DIR}/{{cookiecutter.gtf_files[cookiecutter.species]}} ${ENSEMBL_DIR}
ln -s ${GENOME_DATA_DIR}/genes.tsv ${ENSEMBL_DIR}
ln -s ${GENOME_DATA_DIR}/{{cookiecutter.salmon_index}} ${ENSEMBL_DIR}
ln -s ${GENOME_DATA_DIR}/{{cookiecutter.kallisto_index}} ${ENSEMBL_DIR}

ln -s /srv/data/genome/human/msigdb ${DATA_DIR}

if [[ "{{cookiecutter.species}}" != "human" ]]; then
    ln -s ${GENOME_DATA_DIR}/human_orthologs.tsv ${ENSEMBL_DIR}

    HUMAN_ENSEMBL_DIR=${DATA_DIR}/human_ensembl_{{cookiecutter.ensembl_version}}
    mkdir -p ${HUMAN_ENSEMBL_DIR}
    ln -s /srv/data/genome/human/ensembl-{{cookiecutter.ensembl_version}}/genes.tsv ${HUMAN_ENSEMBL_DIR}  
fi

git init
mv gitignore .gitignore
