from snake_variables import *
from subprocess import call, check_output
from glob import glob
import os


def pick_first_sample():
    test_sample = SAMPLES[0]
    test_species = SPECIES[0]
    test_gtf = GTF_FILE[0]
    test_bam = os.path.join("results/final_bams", "%s.%s.bam" % (test_sample, test_species))
    return test_bam, test_gtf


def run_featurecounts(strand_setting):
    read_counts_folder = "results/read_counts"
    if not os.path.exists(read_counts_folder):
        os.makedirs(read_counts_folder)
    test_bam, test_gtf = pick_first_sample()
    outfile = os.path.join(read_counts_folder, "counts_temp.%s.out" % strand_setting)
    call("%s %s -T 4 -a %s -o %s -s %s %s" % (FEATURECOUNTS_EXECUTABLE, FEATURECOUNTS_PAIRD_END_FLAG, test_gtf, outfile, strand_setting, test_bam), shell=True)
    outfile_old = ".".join([outfile, "summary"])
    outfile_new = ".".join([outfile, "testsummary"])
    call("mv %s %s" % (outfile_old, outfile_new), shell=True)
    return outfile_new


def get_assigned_reads(counts_file):
    counts = check_output("grep Assigned %s | awk '{print $2}'" % counts_file, universal_newlines=True, shell=True)
    counts = [int(i) for i in counts.split('\n') if i][0]
    return counts


def calculate_strandedness(strand_dict):
    strand_zero = strand_dict["0"]
    strand_one = strand_dict["1"]
    strand_two = strand_dict["2"]
    if all(v == 0 for v in list(strand_dict.values())):
        raise Exception('No mapped reads detected')
    difference = (float(strand_one) - float(strand_two))/(float(strand_one) + float(strand_two))
    if abs(difference) < 0.75:
        return "0"
    elif difference >= 0.75:
        return "1"
    elif difference <= -0.75:
        return "2"
    else:
        raise Exception('Can not calculate strandedness from values: one=%s, two=%s, three=%s' % (strand_zero, strand_one, strand_two))

def get_strand():
    if not os.path.exists("strand.txt"):
        return None
    with open("strand.txt", "r") as openfile:
        strand=openfile.read().rstrip()
        return strand

def set_strand(strandedness):
    with open("strand.txt", "w") as openfile:
        openfile.write(strandedness)

def strand_test(picard=False):
    STRAND=get_strand()
    picard_strand_value = {"0": "NONE", "1": "FIRST_READ_TRANSCRIPTION_STRAND", "2": "SECOND_READ_TRANSCRIPTION_STRAND"}
    if STRAND is None:
        strand_assigned_reads = {"0": 0, "1": 0, "2": 0}
        for i in strand_assigned_reads:
            counts_file = run_featurecounts(i)
            assigned = get_assigned_reads(counts_file)
            strand_assigned_reads[i] = assigned
        strandedness = calculate_strandedness(strand_assigned_reads)
        set_strand(strandedness)
    else:
        strandedness = STRAND
    if picard is True:
        return picard_strand_value[strandedness]
    return strandedness


def retrieve_fastqs(sample):
    fastq = glob("data/rnaseq/%s/*{{ cookiecutter.fastq_suffix }}" % sample)
    return fastq


def species_index():
    assert len(SPECIES) == len(MAPPER_INDEX)
    species_index_pairs = []
    for i in range(len(SPECIES)):
        species_index_pair = " ".join([SPECIES[i], MAPPER_INDEX[i]])
        species_index_pairs.append(species_index_pair)
    species_index_pairs = " ".join(species_index_pairs)
    return species_index_pairs

def sample_tsv():
    outfile_name = 'sample.tsv'
    with open(outfile_name, 'w') as outfile:
        tsv_lines = []
        for sample in SAMPLES:
            fastqs = []
            fastq_f = glob(os.path.join(RNASEQ_DIR, sample, ''.join(['*', READ1_SUFFIX, '.{{ cookiecutter.fastq_suffix }}'])))
            fastq_f = ",".join(fastq_f)
            fastqs.append(fastq_f)
            if IS_PAIRED_END == 'yes':
                fastq_r = glob(os.path.join(RNASEQ_DIR, sample, ''.join(['*', READ2_SUFFIX, '.{{ cookiecutter.fastq_suffix }}'])))
                fastq_r = ",".join(fastq_r)
                fastqs.append(fastq_r)
            fastqs = " ".join(fastqs)
            line = " ".join([sample, fastqs])
            tsv_lines.append(line)
        outfile.write("\n".join(tsv_lines))
    return outfile_name

