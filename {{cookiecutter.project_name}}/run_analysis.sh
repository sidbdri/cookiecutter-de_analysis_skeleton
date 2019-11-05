
from glob import glob
from subprocess import call, check_output
import yaml
import os


MAIN_DIR = '~/{{cookiecutter.projects_base}}/{{cookiecutter.project_name}}'
DATA_DIR = os.path.join(MAIN_DIR, 'data')
PICARD_DIR = os.path.join(DATA_DIR, 'picard')
SPECIES=[]
ENSEMBL_DIR=[]
STAR_INDEX=[]
GTF_FILE=[]
SPECIES_PARA=[]
REF_FLAT=[]

{% for s in cookiecutter.species.split(' ') %}
SPECIES+=({{ s }})
ENSEMBL_DIR+=(DATA_DIR/{{ s }}_ensembl_{{cookiecutter.ensembl_version}})
STAR_INDEX+=(DATA_DIR/{{ s }}_ensembl_{{cookiecutter.ensembl_version}}/STAR_indices/{{cookiecutter.assembly_names[s]}}_{{cookiecutter.star_version}})
BOWTIE2_INDEX+=(DATA_DIR/{{ s }}_ensembl_{{cookiecutter.ensembl_version}}/BOWTIE2_indices/{{cookiecutter.assembly_names[s]}}_{{cookiecutter.bowtie2_version}})
GTF_FILE+=(DATA_DIR/{{ s }}_ensembl_{{cookiecutter.ensembl_version}}/{{cookiecutter.gtf_files[s]}})
REF_FLAT+=(PICARD_DATA/{{ s }}/{{cookiecutter.rff_files[s]}})
SALMON_INDEX+=(DATA_DIR/{{ s }}_ensembl_{{cookiecutter.ensembl_version}}/SALMON_indices/{{cookiecutter.assembly_names[s]}}_{{cookiecutter.salmon_version}})
KALLISTO_INDEX+=(DATA_DIR/{{ s }}_ensembl_{{cookiecutter.ensembl_version}}/KALLISTO_indices/{{cookiecutter.assembly_names[s]}}_{{cookiecutter.kallisto_version}})
{% if cookiecutter.data_type == "rnaseq" %}
SPECIES_PARA+=("{{ s }} ${DATA_DIR}/{{ s }}_ensembl_{{cookiecutter.ensembl_version}}/STAR_indices/{{cookiecutter.assembly_names[s]}}_{{cookiecutter.star_version}}")
{% else %}
SPECIES_PARA+=("{{ s }} ${DATA_DIR}/{{ s }}_ensembl_{{cookiecutter.ensembl_version}}/BOWTIE2_indices/{{cookiecutter.assembly_names[s]}}_{{cookiecutter.bowtie2_version}}")
{% endif %}
{% endfor %}

print(SPECIES_PARA)


#STAR_EXECUTABLE=STAR{{cookiecutter.star_version}}
#BOWTIE2_EXECUTABLE=bowtie2-{{cookiecutter.bowtie2_version}}
#SALMON_EXECUTABLE=salmon{{cookiecutter.salmon_version}}
#KALLISTO_EXECUTABLE=kallisto{{cookiecutter.kallisto_version}}
#FASTQC_EXECUTABLE=fastqc{{cookiecutter.fastqc_version}}
#
#{% if cookiecutter.data_type == "rnaseq" %}
#MAPPER_EXECUTABLE=${STAR_EXECUTABLE}
#{% else %}
#MAPPER_EXECUTABLE=${BOWTIE2_EXECUTABLE}
#{% endif %}

#NUM_THREADS_PER_SAMPLE={{cookiecutter.number_threads_per_sample}}
#NUM_TOTAL_THREADS={{cookiecutter.number_total_threads}}
#NUM_PARALLEL_JOBS=$(awk '{print int($1/$2)}' <<< "${NUM_TOTAL_THREADS} ${NUM_THREADS_PER_SAMPLE}")
#THREAD_USING=0
#MEM_USING=0

#qSVA="{{cookiecutter.qSVA}}"
#SAMPLES="{{cookiecutter.rnaseq_samples}}"
#SAMPLE_STRAND_TEST=${SAMPLES%% *}
#PAIRED_END_READ="{{cookiecutter.paired_end_read}}"
#USE_SARGASSO="{{cookiecutter.sargasso}}"


