import pandas as pd
import os

# This script creates a dataframe with the list of proteins presenting a KRAB domain for every organism

accession = snakemake.input

full_data = []

for accession_number in accession:
    with open(accession_number) as reader:
        prot_list = []
        for line in reader.readlines():
            prot_name = line.split('\t')[0].split(' ')[0]
            if prot_name not in prot_list:
                prot_list.append(prot_name)
        full_data.append([accession_number, prot_list, len(prot_list)])

zf_data = pd.DataFrame(full_data, columns=['Accession', 'Protein List', 'KRAB nb'])
zf_data.to_csv(snakemake.output, sep = ';', index= False)
