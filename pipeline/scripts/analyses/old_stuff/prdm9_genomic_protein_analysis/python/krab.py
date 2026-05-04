import pandas as pd
import os
import json

# This script creates a dataframe with the list of proteins presenting a KRAB domain for every organism

#accession = [elt for elt in os.listdir('results/') if elt.startswith('GC')]
full_data = []

with open('../environment_path.json') as f:
    d = json.load(f)

accession = []
with open('data/resources/organisms_data') as reader:
    for line in reader.readlines()[1:]:
        #taxid = int(line.split('\t')[1])
        assembly = line.split('\t')[2]
        accession.append(assembly)
        #dico[taxid] = assembly
        
for accession_number in accession:
    with open(f"results/{accession_number}/Step4_Hmm/tbl/KRAB_processed") as reader:
        prot_list = []
        for line in reader.readlines():
            prot_name = line.split('\t')[0].split(' ')[0]
            if prot_name not in prot_list:
                prot_list.append(prot_name)
        full_data.append([accession_number, prot_list, len(prot_list)])

zf_data = pd.DataFrame(full_data, columns=['Accession', 'Protein List', 'KRAB nb'])
zf_data.to_csv(d["pathGTDriftGlobalResults"]+'prdm9_genomic_protein_analysis/summarized_results/krab_data.csv', sep = ';', index= False)
