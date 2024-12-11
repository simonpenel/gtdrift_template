import pandas as pd
import os
import json
import argparse

parser = argparse.ArgumentParser(description='This script is used to count the number of proteins with 5 or more Zinc finger domains for every organism')

parser.add_argument('-i', '--input_dir', type=str, required=True, help='Input dir path')
parser.add_argument('-o', '--output_file', type=str, required=True, help='Processed file path')

args = parser.parse_args()

input_dir = args.input_dir
output_file = args.output_file



f = open('config.json')
data = json.load(f)
accession = data['assembly_list']
f.close()
# This script is used to count the number of proteins with 5 or more Zinc finger domains for every organism

full_data = []

for accession_number in accession:
    with open(f"{input_dir}genome_assembly/{accession_number}/analyses/prdm9_prot/hmm_search/domtbl/ZF_domains_summary") as reader:
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
zf_data.to_csv(f'{output_file}', sep = ';', index= False)
