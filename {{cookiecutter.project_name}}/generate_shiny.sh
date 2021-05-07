#!/usr/bin/env bash

# Clone github repository with app scripts
git clone https://github.com/sidbdri/shiny-rna-seq

# Create Data dir
PROJECT=`basename "$PWD"`
mkdir -p shiny-rna-seq/Data/$PROJECT

# Copy DE files
cp results/differential_expression/de_gene/de_summary* shiny-rna-seq/Data/$PROJECT/
cp results/differential_expression/de_gene/deseq2_results* shiny-rna-seq/Data/$PROJECT/

# Copy GO and GSA files
cp -R results/differential_expression/gene_set_tests shiny-rna-seq/Data/$PROJECT/gene_set_tests
cp -R results/differential_expression/enrichment_tests shiny-rna-seq/Data/$PROJECT/enrichment_tests

# Deleting big redundant files
for file in `find shiny-rna-seq/Data/microglia_stroke_ckirby/gene_set_tests/mouse/ | grep 'genes_in_sets'`
do
	rm $file
done

