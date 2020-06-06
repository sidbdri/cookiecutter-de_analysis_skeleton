
{% if cookiecutter.sargasso == "yes" %}
snakemake -s Snakefile.multispecies_analysis multiqc
{% else %}
snakemake -s Snakefile.singlespecies_analysis  multiqc
{% endif %}


