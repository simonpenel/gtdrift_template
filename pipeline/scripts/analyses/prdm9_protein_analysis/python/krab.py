import pandas as pd
import os

# This script creates a dataframe with the list of proteins presenting a KRAB domain for every organism

krab_tab_files = snakemake.input
output_file = snakemake.output[0]
prefix = snakemake.params.path
print("prefix = "+prefix)
prefix_length = len(prefix)

full_data = []

for krab_tab_file in krab_tab_files:
    accession_number = krab_tab_file[(prefix_length + 1) :]
    with open(krab_tab_file) as reader:
        prot_list = []
        for line in reader.readlines():
            prot_name = line.split('\t')[0].split(' ')[0]
            if prot_name not in prot_list:
                prot_list.append(prot_name)
        full_data.append([accession_number, prot_list, len(prot_list)])

zf_data = pd.DataFrame(full_data, columns=['Accession', 'Protein List', 'KRAB nb'])
zf_data.to_csv(output_file, sep = ';', index= False)
