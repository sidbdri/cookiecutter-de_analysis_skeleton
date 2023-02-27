import os

HOME_DIR = os.path.expanduser("~")
MAIN_DIR = os.path.join(HOME_DIR, '{{cookiecutter.projects_base}}/{{cookiecutter.project_name}}')
DATA_DIR = os.path.join(MAIN_DIR, 'data')
RNASEQ_DIR=os.path.join(DATA_DIR, 'rnaseq')
PICARD_DIR = os.path.join(DATA_DIR, 'picard')

NUM_THREADS_PER_SAMPLE={{cookiecutter.number_threads_per_sample}}
NUM_TOTAL_THREADS={{cookiecutter.number_total_threads}}

STAR_EXECUTABLE="STAR{{cookiecutter.star_version}}"
BOWTIE2_EXECUTABLE="bowtie2-{{cookiecutter.bowtie2_version}}"
FASTQC_EXECUTABLE="fastqc{{cookiecutter.fastqc_version}}"
FEATURECOUNTS_EXECUTABLE="featureCounts{{cookiecutter.featurecounts_version}}"
PICARD_EXECUTABLE="/opt/picard-tools-{{cookiecutter.picard_version}}/picard.jar"
SALMON_EXECUTABLE="salmon{{cookiecutter.salmon_version}}"

{% if cookiecutter.paired_end_read == "yes" %}
FEATURECOUNTS_PAIRED_END_FLAG = '-p --countReadPairs'
{% else %}
FEATURECOUNTS_PAIRED_END_FLAG = ''
{% endif %}

SPECIES=[]
ENSEMBL_DIR=[]
MAPPER_INDEX=[]
BOWTIE_INDEX=[]
SALMON_INDEX={}
GTF_FILE=[]
REF_FLAT=[]

DATA_TYPE = "{{cookiecutter.data_type}}"
if DATA_TYPE == "rnaseq":
     MAPPER_EXECUTABLE=STAR_EXECUTABLE
else:
    MAPPER_EXECUTABLE=BOWTIE2_EXECUTABLE

{% for s in cookiecutter.species.split(' ') %}
SPECIES.append("{{ s }}")
ENSEMBL_DIR.append("%s/{{ s }}_ensembl_{{cookiecutter.ensembl_version}}" % DATA_DIR)
if DATA_TYPE == "rnaseq":
    MAPPER_INDEX.append("%s/{{ s }}_ensembl_{{cookiecutter.ensembl_version}}/STAR_indices/{{cookiecutter.assembly_names[s]}}_{{cookiecutter.star_version}}" % DATA_DIR)
else:
    MAPPER_INDEX.append("%s/{{ s }}_ensembl_{{cookiecutter.ensembl_version}}/BOWTIE2_indices/{{cookiecutter.assembly_names[s]}}_{{cookiecutter.bowtie2_version}}" % DATA_DIR)
GTF_FILE.append("%s/{{ s }}_ensembl_{{cookiecutter.ensembl_version}}/{{cookiecutter.gtf_files[s]}}" % DATA_DIR)
REF_FLAT.append("%s/{{ s }}/{{cookiecutter.rff_files[s]}}" % PICARD_DIR)
SALMON_INDEX["{{ s }}"]=("%s/{{ s }}_ensembl_{{cookiecutter.ensembl_version}}/SALMON_indices/{{cookiecutter.assembly_names[s]}}_{{cookiecutter.salmon_version}}" % DATA_DIR)
{% endfor %}


STRATEGY="{{cookiecutter.strategy}}"


IS_PAIRED_END="{{cookiecutter.paired_end_read}}"
gtf_dict = dict(zip(SPECIES, GTF_FILE))
SAMPLES="{{cookiecutter.rnaseq_samples}}"
SAMPLES=SAMPLES.split()

READ1_SUFFIX = ''
READ2_SUFFIX = ''
read_identifiers =  ["{{ cookiecutter.read1_identifier }}", "{{ cookiecutter.read2_identifier }}"]

if IS_PAIRED_END == "yes":
    READ1_SUFFIX = '*'
    READ1_SUFFIX = ''
else:
    READ1_SUFFIX = read_identifiers[0]
    READ2_SUFFIX = read_identifiers[1]


