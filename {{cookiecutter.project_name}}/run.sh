
NUM_TOTAL_THREADS={{cookiecutter.number_total_threads}}

{% if cookiecutter.sargasso == "yes" %}
snakemake -s Snakefile.multispecies_analysis bams -j $NUM_TOTAL_THREADS
snakemake -s Snakefile.multispecies_analysis multiqc -j $NUM_TOTAL_THREADS
{% else %}
snakemake -s Snakefile.singlespecies_analysis  bams -j $NUM_TOTAL_THREADS
snakemake -s Snakefile.singlespecies_analysis  multiqc -j $NUM_TOTAL_THREADS
{% endif %}


