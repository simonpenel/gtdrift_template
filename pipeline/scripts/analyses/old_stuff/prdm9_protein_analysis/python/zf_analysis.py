import pandas as pd
import os

# This script is used to count the number of proteins with 5 or more Zinc finger domains for every organism

full_data = []

zf_domain_summary_files = snakemake.input
output_file = snakemake.output[0]
prefix = snakemake.params.path
prefix_length = len(prefix)

for zf_domain_summary_file in zf_domain_summary_files:
    accession_number = zf_domain_summary_file[(prefix_length + 1) :].split("/")[0]
    with open(zf_domain_summary_file) as reader:
        lines = reader.readlines()[1:]
        i = 0
        prot_count = 0
        while i < len(lines):
            domains_nb = int(lines[i].split('\t')[10])
            if domains_nb >= 5:
                prot_count += 1
            i += domains_nb
        full_data.append([accession_number, prot_count])

zf_data = pd.DataFrame(full_data, columns=['Accession', '5+ ZF'])
zf_data.to_csv(output_file, sep = ';', index= False)
