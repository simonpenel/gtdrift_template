import pandas as pd
import os
import json
import argparse

parser = argparse.ArgumentParser(description='This script takes the previously created krab_data dataframe, that contains the list of every protein presenting a krab domain for every organism, and checks for every entries in all of the lists whether these proteins also carry a zinc finger domain or not')
parser.add_argument('-i', '--input_dir', type=str, required=True, help='Input dir path')
#parser.add_argument('-o', '--output_dir', type=str, required=True, help='Processed file path')
parser.add_argument('-o', '--output_file', type=str, required=True, help='Output file path')
parser.add_argument('-k', '--krab_data', type=str, required=True, help='krab_data.csv path')

args = parser.parse_args()

input_dir = args.input_dir
#output_dir = args.output_dir
output_file = args.output_file
krab_data_file = args.krab_data

f = open('config.json')
data = json.load(f)
accession = data['assembly_list']
f.close()
# This script takes the previously created krab_data dataframe, that contains the list of every protein presenting a krab domain for every organism,
# and checks for every entries in all of the lists whether these proteins also carry a zinc finger domain or not.

full_data = []

#krab_data = pd.read_csv(f'{output_dir}analyses_summaries/table_results/krab_data.csv', sep=';')
krab_data = pd.read_csv(f'{krab_data_file}', sep=';')
for accession_number in accession:
    with open(f"{input_dir}genome_assembly/{accession_number}/analyses/prdm9_prot/hmm_search/tbl/ZF_tabulated") as reader:
        prot_list = []
        for line in reader.readlines():
            prot_name = line.split('\t')[0].split(' ')[0]
            if prot_name in krab_data.loc[krab_data['Accession'] == accession_number, 'Protein List'].iloc[0]: 
                prot_list.append(prot_name)
        full_data.append([accession_number, prot_list, len(prot_list)])

zf_data = pd.DataFrame(full_data, columns=['Accession', 'KRAB+ZF protein list', 'KRAB+ZF nb'])
zf_data.to_csv(f'{output_dir}analyses_summaries/table_results/krabzf_data.csv', sep= ';', index = False)
