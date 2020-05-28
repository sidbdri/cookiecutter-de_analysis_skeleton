
{% if cookiecutter.sargasso == "yes" %}
snakemake -s Snakefile.run_analysis sargasso 
{% endif %}


snakemake -s Snakefile.run_analysis bams
snakemake -s Snakefile.run_analysis  multiqc

