#!/bin/bash

set -o nounset
set -o errexit
set -o xtrace

MAIN_DIR=${HOME}/{{cookiecutter.projects_base}}/{{cookiecutter.project_name}}
RESULTS_DIR=${MAIN_DIR}/results

##### Record software version information

mkdir -p ${RESULTS_DIR}

README=${RESULTS_DIR}/README
echo "Ensembl version: {{cookiecutter.ensembl_version}}" > ${README}
echo "Git commit: $(git log --pretty=format:'%H' -n 1)" >> ${README}
echo "Software versions:" >> ${README}
fastqc --version >> ${README}
STAR --version >> ${README}
echo $(featureCounts -v 2>&1) >> ${README}
echo "Salmon $(salmon --version 2>&1)" >> ${README}
kallisto version >> ${README}
multiqc --version >> ${README}
get_gene_lengths --version >> ${README}
Rscript load_packages.R 2>&1 | sed -n '/R version/,$p' >> ${README}
