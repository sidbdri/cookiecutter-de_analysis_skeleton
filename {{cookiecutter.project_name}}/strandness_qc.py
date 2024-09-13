# https://github.com/sidbdri/cookiecutter-de_analysis_skeleton/issues/226
import json

with open('strand.txt', 'r') as f:
    strandedness = int(f.readline().strip())

with open('./multiqc_data/multiqc_data.json', 'r') as file:
    data = json.load(file)

strandness_qc = []

# Loop through the dict to find PCT_R2_TRANSCRIPT_STRAND_READS
for key, value in data['report_general_stats_data'][1].items():
    if "PCT_R2_TRANSCRIPT_STRAND_READS" in value:
        # print(f"{key}: {value['PCT_R2_TRANSCRIPT_STRAND_READS']}")
        pct_r2 = float(value['PCT_R2_TRANSCRIPT_STRAND_READS'])
        # print(pct_r2)
        if strandedness == 0 and (pct_r2 <= 45 or pct_r2 >= 55):
            strandness_qc.append(f"{key}\t{pct_r2}")
        elif strandedness == 1 and pct_r2 >= 5:
            strandness_qc.append(f"{key}\t{pct_r2}")
        elif strandedness == 2 and pct_r2 <= 95:
            strandness_qc.append(f"{key}\t{pct_r2}")

# Check if there's any strandedness QC issue
if len(strandness_qc) > 0:
    print("Error: The following sample may have the wrong strandedness {} based on the PCT_R2_TRANSCRIPT_STRAND_READS:".format(strandedness))
    print("\n".join(strandness_qc))
    exit(1)



