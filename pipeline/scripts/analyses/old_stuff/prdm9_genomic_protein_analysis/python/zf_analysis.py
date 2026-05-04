import pandas as pd
import os
import json

# This script is used to count the number of proteins with 5 or more Zinc finger domains for every organism

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

for accession_number in accession:
    with open(f"results/{accession_number}/Step4_Hmm/domtbl/ZF_domains_summary") as reader:
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
zf_data.to_csv(d["pathGTDriftGlobalResults"]+'prdm9_genomic_protein_analysis/summarized_results/zf_count.csv', sep = ';', index= False)
