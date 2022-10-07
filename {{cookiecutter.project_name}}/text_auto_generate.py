import subprocess
import math

def millify(n):
    millnames = ['',' Thousand',' Million',' Billion',' Trillion']
    n = float(n)
    millidx = max(0,min(len(millnames)-1,
                        int(math.floor(0 if n == 0 else math.log10(abs(n))/3))))

    return '{:.0f}{}'.format(n / 10**(3 * millidx), millnames[millidx])

def get_average_read_count(samples):
    number_of_reads = []

    for sample in samples:
        with open(f"results/mapped_reads/{sample}.log.out") as iH:
            for line in iH:
                if line.strip().startswith("Number of input reads"):
                    number_of_reads.append(int(line.strip().split("|")[-1].strip()))

    return millify(round(sum(number_of_reads)/len(number_of_reads), -6))

def get_average_read_length(samples):
    length_of_reads = []

    for sample in samples:
        with open(f"results/mapped_reads/{sample}.log.out") as iH:
            for line in iH:
                if line.strip().startswith("Average input read length"):
                    length_of_reads.append(int(line.strip().split("|")[-1].strip()))

    return int(sum(length_of_reads)/len(length_of_reads))

def r_versions(to_search):
    packages_dict = {}
    packages = []
    r_script_output = subprocess.run(["Rscript", "load_packages.R"], stderr=subprocess.DEVNULL, stdout=subprocess.PIPE)
    for line in r_script_output.stdout.decode().split("\n"):
        if line.lstrip().startswith("["):
            packages += line.split("]")[-1].split()

        for search in to_search:
            for package in packages:
                if package.startswith(search):
                    packages_dict[search] = package.split("_")[-1]
    
    return packages_dict

packages_dict = r_versions(["DESeq2", "limma", "topGO"])

rnaseq_samples="{{cookiecutter.rnaseq_samples}}"
rnaseq_samples=rnaseq_samples.split()
#rnaseq_samples=["SWAP12HR1","SWAP12HR2","SWAP12HR3","SWAP12HR4","SWAPdark1","SWAPdark2","SWAPdark3","SWAPdark4","SWAPdark5","SWAPdark6","SWAPlight1","SWAPlight2","SWAPlight3","SWAPlight4","SWAPlight5","SWAPlight6","WT12HR1","WT12HR2","WT12HR3","WT12HR4","WTdark1","WTdark2","WTdark3","WTdark4","WTdark5","WTlight1","WTlight2","WTlight3","WTlight4","WTlight5","WTlight6"]

is_paired_end="{{cookiecutter.paired_end_read}}"
#is_paired_end="yes"
if is_paired_end=="yes":
    read_type="paired-end"
else:
    read_type="single-end"

species = "{{cookiecutter.species}}"
species = species.split()
#species = ["mouse","human"]
assembly = {}
for s in species:
    if s=="human":
        assembly["human"]="{{cookiecutter.human_assembly}}"
    elif s=="mouse":
        assembly["mouse"]="{{cookiecutter.mouse_assembly}}"
    elif s=="rat":
        assembly["rat"]="{{cookiecutter.rat_assembly}}"

assembly_string = ""
for key in assembly:
    assembly_string = assembly_string + key + "(" + assembly[key] + "), "

assembly_string=assembly_string.rstrip(assembly_string[-2])

ensembl_version="{{cookiecutter.ensembl_version}}"
msigdb_version="{{cookiecutter.msigdb_version}}"
star_version="{{cookiecutter.star_version}}"
featurecounts_version="{{cookiecutter.featurecounts_version}}"

print ("RNA sequencing was performed using FILL_THIS_IN library preparation along with next-generation sequencing on the FILL_THIS_IN platform; sequencing was carried out by FILL_THIS_IN. samples were sequenced to a depth of approximately", get_average_read_count(rnaseq_samples), get_average_read_length(rnaseq_samples), "base pair,", read_type, "reads. The reads were mapped to the primary assembly of the", assembly_string, "reference genome contained in Ensembl release", ensembl_version, ",using the STAR RNA-seq aligner, version", star_version, "[1]. Tables of per-gene read counts were generated from the mapped reads with featureCounts, version", featurecounts_version, "[2]. Differential gene expression was performed in R using DESeq2, version", packages_dict["DESeq2"], "[3]. Gene ontology enrichment analysis was performed using topGO, version", packages_dict["topGO"], "[4]. Gene set testing was then performed using Camera [5] from the R package limma, version", packages_dict["limma"], "[6], using gene sets from the Molecular Signatures Database, version", msigdb_version," (https://www.gsea-msigdb.org/gsea/msigdb/).")

print ()
print ("[1] Dobin, A. et al. STAR: ultrafast universal RNA-seq aligner. Bioinformatics 29, 15–21 (2013).")
print ("[2] Liao, Y., Smyth, G. K. & Shi, W. featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics 30, 923–930 (2014).")
print ("[3] Love, M. I., Huber, W. & Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology 15, 550 (2014).")
print ("[4] Alexa, A., Rahnenfuhrer, J. & Lengauer, T. Improved scoring of functional groups from gene expression data by decorrelating GO graph structure. Bioinformatics 22, 1600–1607 (2006).")
print ("[5] Wu, D. & Smyth, G. K. Camera: a competitive gene set test accounting for inter-gene correlation. Nucleic Acids Research 40, e133 (2012).")
print ("[6] Ritchie et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Research 43, e47 (2015).")
print ()


