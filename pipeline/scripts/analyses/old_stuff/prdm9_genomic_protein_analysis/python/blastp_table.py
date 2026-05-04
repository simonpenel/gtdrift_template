from ete3 import NCBITaxa # type: ignore
import pandas as pd
import os
import sys

"""
This script creates an output file named blastp_results.csv. It is a tabular formating of the previously created blastp_summary.txt file. Those two outputs are created separately due to
snakemake iteration constraints, and optimisation constraints.
"""

ncbi = NCBITaxa()
data = []

with open('results/BLASTP_results/blastp_summary.txt', 'r+') as reader:
    for line in reader.readlines():
        line = line.strip()
        # Get taxid and taxonomy data
        if line.startswith('>') == True:
            taxid = line.lstrip('>')
            try:
                lineage = ncbi.get_lineage(taxid)
                taxonomy = ncbi.get_taxid_translator(lineage)
                species = ncbi.get_taxid_translator([taxid])[int(taxid)]
            except ValueError as err:
                print("Error : ",format(err))
                lineage = []
                taxonomy = {}
                species = ""
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
            if 'ERROR' in prot_id:
                status = 'Error during genewisedb'
            else:
                status = 'Valid'
            
            # getting info for prot to check for pseudogene status
            pseudo = pd.read_csv('results/' + prot_id.split('-')[0] + '/Result_tables/summary_table_' + prot_id.split('-')[0] + '.csv', sep=';')
            psgene = int((pseudo.loc[pseudo['SeqID'] == prot_id, 'Pseudogene (HMMER)'].iloc[0] == 'Yes'))
            trunc  = int((pseudo.loc[pseudo['SeqID'] == prot_id, 'ZF Truncated'].iloc[0] == 'Yes'))
            
            data.append([taxid, species, prot_id, bit_score, ratio, set, krab, ssxrd, zf, match, psgene, trunc, status]) # You can add the superorder or any chosen taxonomy level you want to have in the final dataframe
output_df = pd.DataFrame(data, columns=['Taxid', 'Species_name', 'Protein_ID', 'Bit_score', 'Ratio', 'SET', 'KRAB', 'SSXRD', 'ZF', 'Best_Match', 'Is_Pseudogene', 'ZF_Truncated', 'Status']) # Add the corresponding column names at the same position you added the value in the list above.

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
output_df.to_csv('results/BLASTP_results/blastp_results.csv', sep=';')

# os.system(f"mkdir -p test_seq_pour_align/")
# for index, row in df_fusion.iterrows():
#     accession = row['Accession']
#     protein_id = row['Protein ID']
#     os.system(f"blastdbcmd -db data/ncbi/{accession}/protdb -entry {protein_id} -out test_seq_pour_align/{protein_id}.fa")
        
