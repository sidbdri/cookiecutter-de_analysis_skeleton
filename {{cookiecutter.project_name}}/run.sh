
NUM_THREADS_PER_SAMPLE={{cookiecutter.number_threads_per_sample}}
NUM_TOTAL_THREADS={{cookiecutter.number_total_threads}}
NUM_PARALLEL_JOBS=`expr ${NUM_TOTAL_THREADS} / ${NUM_THREADS_PER_SAMPLE}`

{% if cookiecutter.sargasso == "yes" %}
snakemake -s Snakefile.multispecies_analysis bams -j $NUM_PARALLEL_JOBS
snakemake -s Snakefile.multispecies_analysis multiqc -j $NUM_PARALLEL_JOBS
{% else %}
snakemake -s Snakefile.singlespecies_analysis  bams -j $NUM_PARALLEL_JOBS
snakemake -s Snakefile.singlespecies_analysis  multiqc -j $NUM_PARALLEL_JOBS
{% endif %}


