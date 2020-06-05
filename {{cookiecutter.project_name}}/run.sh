
{% if cookiecutter.sargasso == "yes" %}
snakemake -s Snakefile.run_analysis bams_mixed_species
{% else %}
snakemake -s Snakefile.run_analysis bams
{% endif %}
snakemake -s Snakefile.run_analysis  multiqc

