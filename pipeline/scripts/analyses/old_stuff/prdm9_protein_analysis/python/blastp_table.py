from ete3 import NCBITaxa
import pandas as pd
import os
import sys
import argparse

parser = argparse.ArgumentParser(description='Creates an output file named blastp_results.csv')
parser.add_argument('-i', '--input_dir', type=str, required=True, help='Input dir path')
args = parser.parse_args()
input_dir = args.input_dir


"""
This script creates an output file named blastp_results.csv. It is a tabular formating of the previously created blastp_summary.txt file. Those two outputs are created separately due to
snakemake iteration constraints, and optimisation constraints.
"""

ncbi = NCBITaxa()
data = []

with open(f'{input_dir}../data_results_per_assembly/organisms_data') as reader:
    dico = {}
    for line in reader.readlines()[1:]:
        taxid = int(line.split('\t')[1])
        species = line.split('\t')[0]
        dico[taxid] = species
        
with open(f'{input_dir}analyses_summaries/BLASTP_results/blastp_summary.txt', 'r+') as reader:
    for line in reader.readlines():
        line = line.strip()
        # Get taxid and taxonomy data
        if line.startswith('>') == True:
            taxid = line.lstrip('>')
            if int(taxid) in dico:
                species = dico[int(taxid)]
            else:
                species = "Not found"
                print("Taxid "+taxid+" not found")
            #try:
            #    species = dico[int(taxid)]
                #lineage = ncbi.get_lineage(taxid)
                #taxonomy = ncbi.get_taxid_translator(lineage)
                #species = ncbi.get_taxid_translator([taxid])[int(taxid)]
            #except ValueError as err:
            #    print("Error : ",format(err))
                #lineage = []
                #taxonomy = {}
             #   species = "Not found"
            # superorder = taxonomy[lineage[18]] # this line was used when working on insects. It can be modified to keep the chosen tanonomy level (change lineage index) or simply ignored.
        # get proteic domains data
        elif line.startswith('<') == True:
            set = int(line.split('\t')[1])
            krab = int(line.split('\t')[2])
            ssxrd = int(line.split('\t')[3])
            zf = int(line.split('\t')[4])
        # get BLASTP results
        else:
            line = line.split('\t')
            # if the best match was not PRDM9
            if len(line) == 2:
                prot_id = line[0]
                bit_score = 0
                ratio = 0
                match = line[-1]
            # if it was
            else:                
                prot_id = line[0]
                bit_score = line[1]
                ratio = line[2]
                match = line[-1]
            data.append([taxid, species, prot_id, bit_score, ratio, set, krab, ssxrd, zf, match]) # You can add the superorder or any chosen taxonomy level you want to have in the final dataframe

output_df = pd.DataFrame(data, columns=['Taxid', 'Species_name', 'Protein ID', 'Bit score', 'Ratio', 'SET', 'KRAB', 'SSXRD', 'ZF', 'Best_Match']) # Add the corresponding column names at the same position you added the value in the list above.

# Dropping useless old indexes column
output_df = output_df.loc[:, ~output_df.columns.str.contains('^Unnamed')]

# # Lecture des données taxonomiques et ajout d'une colonne de synthèse des noms d'espèces
# taxonomy = pd.read_csv("sorted_taxonomy.csv", sep=";")
# for index, row in taxonomy.iterrows():
#     last_valid_col = taxonomy.iloc[index].last_valid_index()
#     taxonomy.at[index, 'Species_name'] = taxonomy.loc[index, last_valid_col]
# taxonomy = taxonomy[['Accession', 'Species_name']]
# taxonomy.rename(columns={'Species_name': 'Species_name2'}, inplace=True)

# output_df = output_df.sort_values(by=['Superorder', 'Species_name']) 
output_df.to_csv(f'{input_dir}/analyses_summaries/BLASTP_results/blastp_results.csv', sep=';')

# os.system(f"mkdir -p test_seq_pour_align/")
# for index, row in df_fusion.iterrows():
#     accession = row['Accession']
#     protein_id = row['Protein ID']
#     os.system(f"blastdbcmd -db data/ncbi/{accession}/protdb -entry {protein_id} -out test_seq_pour_align/{protein_id}.fa")
        
