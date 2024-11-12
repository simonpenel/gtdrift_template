import pandas as pd
import os
import json

# This script takes the previously created krab_data dataframe, that contains the list of every protein presenting a krab domain for every organism,
# and checks for every entries in all of the lists whether these proteins also carry a zinc finger domain or not.

with open('../environment_path.json') as f:
    d = json.load(f)
    
#accession = [elt for elt in os.listdir('results/') if elt.startswith('GC')]
full_data = []

accession = []
with open('data/resources/organisms_data') as reader:
    for line in reader.readlines()[1:]:
        #taxid = int(line.split('\t')[1])
        assembly = line.split('\t')[2]
        accession.append(assembly)
        #dico[taxid] = assembly
        
krab_data = pd.read_csv(d["pathGTDriftGlobalResults"]+'prdm9_genomic_protein_analysis/summarized_results/krab_data.csv', sep=';')
for accession_number in accession:
    with open(f"results/{accession_number}/Step4_Hmm/tbl/ZF_processed") as reader:
        prot_list = []
        for line in reader.readlines():
            prot_name = line.split('\t')[0].split(' ')[0]
            if prot_name in krab_data.loc[krab_data['Accession'] == accession_number, 'Protein List'].iloc[0]: 
                prot_list.append(prot_name)
        full_data.append([accession_number, prot_list, len(prot_list)])

zf_data = pd.DataFrame(full_data, columns=['Accession', 'KRAB+ZF protein list', 'KRAB+ZF nb'])
zf_data.to_csv(d["pathGTDriftGlobalResults"]+'prdm9_genomic_protein_analysis/summarized_results/krabzf_data.csv', sep= ';', index = False)
