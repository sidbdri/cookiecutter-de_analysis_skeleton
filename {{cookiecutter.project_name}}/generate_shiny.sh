#!/usr/bin/env bash

# Clone github repository with app scripts
git clone https://github.com/sidbdri/shiny-rna-seq

# Create Data dir
mkdir -p shiny-rna-seq/Data/{{cookiecutter.projects_base}}

# Copy DE files
cp results/differential_expression/de_gene/de_summary* shiny-rna-seq/Data/{{cookiecutter.projects_base}}/
cp results/differential_expression/de_gene/deseq2_results* shiny-rna-seq/Data/{{cookiecutter.projects_base}}/

# Copy GO and GSA files
cp -R results/differential_expression/gene_set_tests shiny-rna-seq/Data/{{cookiecutter.projects_base}}/gene_set_tests
cp -R results/differential_expression/enrichment_tests shiny-rna-seq/Data/{{cookiecutter.projects_base}}/enrichment_tests

# Deleting big redundant files
find shiny-rna-seq/Data/{{cookiecutter.projects_base}}/gene_set_tests/ -name "*genes_in_sets*" -exec rm {} \;
