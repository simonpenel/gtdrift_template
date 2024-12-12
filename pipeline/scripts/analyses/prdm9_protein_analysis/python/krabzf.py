import pandas as pd
import os

krab_data_file = snakemake.input.krab
zf_tab_files = snakemake.input.zf_tabs
output_file = snakemake.output[0]
prefix = snakemake.params.path
prefix_length = len(prefix)


# This script takes the previously created krab_data dataframe, that contains the list of every protein presenting a krab domain for every organism,
# and checks for every entries in all of the lists whether these proteins also carry a zinc finger domain or not.

full_data = []

#krab_data = pd.read_csv(f'{output_dir}analyses_summaries/table_results/krab_data.csv', sep=';')
krab_data = pd.read_csv(krab_data_file, sep=';')
for zf_tab_file in zf_tab_files:
    print("open " + zf_tab_file)
    accession_number = zf_tab_file[(prefix_length + 1) :].split("/")[0]
    with open(zf_tab_file) as reader:
        prot_list = []
        for line in reader.readlines():
            prot_name = line.split('\t')[0].split(' ')[0]
            if prot_name in krab_data.loc[krab_data['Accession'] == accession_number, 'Protein List'].iloc[0]: 
                prot_list.append(prot_name)
        full_data.append([accession_number, prot_list, len(prot_list)])

zf_data = pd.DataFrame(full_data, columns=['Accession', 'KRAB+ZF protein list', 'KRAB+ZF nb'])
zf_data.to_csv(output_file, sep= ';', index = False)
