import pandas as pd
import os
import json
import argparse

parser = argparse.ArgumentParser(description='Reads hmm_search processed files and domains processed files and creates an overview table in the csv format')
parser.add_argument('-i', '--input_dir', type=str, required=True, help='Input dir path')
parser.add_argument('-o', '--output_file', type=str, required=True, help='Processed file path')

args = parser.parse_args()
output_file = args.output_file
input_dir = args.input_dir

f = open('config.json')
data = json.load(f)
accession = data['assembly_list']
f.close()
# This script creates a dataframe with the list of proteins presenting a KRAB domain for every organism

full_data = []

for accession_number in accession:
    with open(f"{input_dir}genome_assembly/{accession_number}/analyses/prdm9_prot/hmm_search/tbl/KRAB_tabulated") as reader:
        prot_list = []
        for line in reader.readlines():
            prot_name = line.split('\t')[0].split(' ')[0]
            if prot_name not in prot_list:
                prot_list.append(prot_name)
        full_data.append([accession_number, prot_list, len(prot_list)])

zf_data = pd.DataFrame(full_data, columns=['Accession', 'Protein List', 'KRAB nb'])
zf_data.to_csv(f'{output_file}', sep = ';', index= False)
