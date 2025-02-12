#!/bin/bash

export WORKON_HOME=${HOME}/{{cookiecutter.virtualenv_home}}
export PROJECT_HOME=${HOME}/{{cookiecutter.projects_base}}
source /usr/local/bin/virtualenvwrapper.sh

## we do some init setup such as recording versions
[[ -f "./project_setup.sh" ]] && bash ./project_setup.sh

mkproject -f -p python3.8 {{cookiecutter.project_name}}

source config.sh

pip install --requirement requirements.txt
pip install git+https://github.com/sidbdri/transcript-utils.git@${transcript_utils_hash}
{% if cookiecutter.sargasso == "yes" %}
pip install git+https://github.com/biomedicalinformaticsgroup/Sargasso.git@${sargasso_hash}
{% endif %}

## Clone sidbdri-utils package
rm -rf sidbdri-utils && git clone https://github.com/sidbdri/sidbdri-utils.git
(cd sidbdri-utils && git checkout -q ${sidbdri_utils_hash})
source sidbdri-utils/includes.sh

DATA_DIR=data
RNASEQ_DIR=${DATA_DIR}/rnaseq
PICARD_DATA=${DATA_DIR}/picard
PICARD=/opt/picard-tools-{{cookiecutter.picard_version}}/picard.jar

# Check sample names for special characters
# https://github.com/sidbdri/cookiecutter-de_analysis_skeleton/issues/98
if [[ `echo {{cookiecutter.rnaseq_samples}} | grep -P '[\t\n.]'` != ''  ]];then
    echo "Error: Please make sure sample names don't contain special characters."
    exit 1
fi

mkdir -p ${RNASEQ_DIR}
for sample in {{cookiecutter.rnaseq_samples}}; do
    ln -s -T {{cookiecutter.rnaseq_samples_dir}}/$sample ${RNASEQ_DIR}/$sample
done

# we make sure the sample are in fastq.gz format
# wrap it under if block to hopefully gain some compatibility for old project
if [[ -f 'sidbdri-utils/scripts/genozip_to_fastqgz.sh' ]]; then
    bash sidbdri-utils/scripts/genozip_to_fastqgz.sh -d {{cookiecutter.rnaseq_samples_dir}}
fi

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
ln -s ${GENOME_DATA_DIR}/BOWTIE2_indices ${ENSEMBL_DIR}
ln -s ${GENOME_DATA_DIR}/*orthologs.tsv  ${ENSEMBL_DIR}
ln -s ${GENOME_DATA_DIR}/msigdb ${ENSEMBL_DIR}

## Instead of using the genes.tsv from ${GENOME_DATA_DIR}, we download it from Ensembl
for i in 1 2 3 4 5; do
  download_gene_tb {{ s }} {{cookiecutter.ensembl_version}} > ${ENSEMBL_DIR}/genes.tsv
  if [[ $(wc -l < ${ENSEMBL_DIR}/genes.tsv) -gt 1000 ]]; then
    break
  fi
  echo "error downloading ${ENSEMBL_DIR}/genes.tsv, retry "${i}" out of 5 times."
  sleep 1;
done

if [[ $(wc -l < ${ENSEMBL_DIR}/genes.tsv) -le 1000 ]]; then
  echo "Error: ${ENSEMBL_DIR}/genes.tsv download failed. Exiting!"
  exit 1
fi

# Generating refFlat file for Picard RNA-seq metrics
generate_picard_refFlat ${PICARD_DATA}/{{ s }} {{ s }} {{cookiecutter.ensembl_version}} ${GTF_FILE} &

# Link gene_length.csv if exist in genome folder
[[ -f "${GENOME_DATA_DIR}/gene_lengths.csv" ]] && ln -s ${GENOME_DATA_DIR}/gene_lengths.csv ${ENSEMBL_DIR}
[[ -f "${GENOME_DATA_DIR}/tx2gene.tsv" ]] && ln -s ${GENOME_DATA_DIR}/tx2gene.tsv ${ENSEMBL_DIR}
{% endfor %}


{% if "human" not in cookiecutter.species.split(' ') %}
HUMAN_ENSEMBL_DIR=${DATA_DIR}/human_ensembl_{{cookiecutter.ensembl_version}}
mkdir -p ${HUMAN_ENSEMBL_DIR}

## Instead of using the genes.tsv from ${GENOME_DATA_DIR}, we download it from Ensembl
download_gene_tb human {{cookiecutter.ensembl_version}} > ${HUMAN_ENSEMBL_DIR}/genes.tsv
{% endif %}

## if we have a lock file but the project library is not set up, we set it up
echo "checking renv setup..."
[ ! -d "renv/library" ] && [ -f "renv.lock" ] && echo "setting up renv" &&  R -e 'renv::restore(prompt = FALSE)' && echo "new renv initialized base on .lock file" || echo "renv has already been setup. do nothing."