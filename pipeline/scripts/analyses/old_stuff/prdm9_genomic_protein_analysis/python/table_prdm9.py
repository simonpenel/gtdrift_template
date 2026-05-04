import pandas as pd
import json
pd.options.mode.copy_on_write = True

"""
This script creates a summary table of domain presence absence in candidate genes for every organism
"""

with open('../environment_path.json') as f:
    d = json.load(f)
    
taxonomy = pd.read_csv('data/resources/sorted_taxonomy.csv', sep= ';')
taxonomy = taxonomy[['Taxid', 'Accession']]
donnees = pd.read_csv("results/BLASTP_results/blastp_results.csv", sep=";")
donnees = donnees.loc[:, ~donnees.columns.str.contains('^Unnamed')]

donnees['PRDM9'] = False
for index, row in donnees.iterrows():
    if row['Best_Match'] == 'PRDM9':
        donnees.at[index, 'PRDM9'] = True 
    else:
        donnees.at[index, 'PRDM9'] = False 

donnees['boolZF'] = donnees['ZF'].apply(lambda x: 0 if x == 0 else 1)
donnees['Complete_PRDM9'] = ((donnees['SET'] + donnees['KRAB'] + donnees['SSXRD'] + donnees['boolZF'] == 4) & (donnees['Is_Pseudogene'] == 0) & (donnees['PRDM9'] == True)).astype(int)
donnees['Incomplete_PRDM9'] = ((donnees['Complete_PRDM9'] == 0) & (donnees['Is_Pseudogene'] == 0) & (donnees['PRDM9'] == True)).astype(int)
donnees['SET+KRAB+SSXRD'] = ((donnees['Complete_PRDM9'] == 0) & (donnees['KRAB'] + donnees['SSXRD'] == 2) & (donnees['Is_Pseudogene'] == 0) & (donnees['PRDM9'] == True)).astype(int)
donnees['SET+KRAB+ZF'] = ((donnees['Complete_PRDM9'] == 0) & (donnees['KRAB'] + donnees['boolZF'] == 2) & (donnees['Is_Pseudogene'] == 0) & (donnees['PRDM9'] == True)).astype(int)
donnees['SET+KRAB'] = ((donnees['Complete_PRDM9'] == 0) & (donnees['SET+KRAB+ZF'] == 0) & (donnees['SET+KRAB+SSXRD'] == 0) & (donnees['KRAB'] == 1) & (donnees['Is_Pseudogene'] == 0) & (donnees['PRDM9'] == True)).astype(int)
donnees['KRAB+ZF'] = ((donnees['Complete_PRDM9'] == 0) & (donnees['SET+KRAB+ZF'] == 0) & (donnees['SET+KRAB+SSXRD'] == 0) & (donnees['KRAB'] == 1) & (donnees['boolZF'] == 1) & (donnees['Is_Pseudogene'] == 0) & (donnees['PRDM9'] == True)).astype(int)
donnees['SET+SSXRD+ZF'] = ((donnees['Complete_PRDM9'] == 0) & (donnees['SSXRD'] + donnees['boolZF'] == 2) & (donnees['Is_Pseudogene'] == 0) & (donnees['PRDM9'] == True)).astype(int)
donnees['SET+SSXRD'] = ((donnees['Complete_PRDM9'] == 0) & (donnees['SET+SSXRD+ZF'] == 0) & (donnees['SET+KRAB+SSXRD'] == 0) & (donnees['SSXRD'] == 1) & (donnees['Is_Pseudogene'] == 0) & (donnees['PRDM9'] == True)).astype(int)
donnees['SSXRD+ZF'] = ((donnees['Complete_PRDM9'] == 0) & (donnees['SET+SSXRD+ZF'] == 0) & (donnees['SET+KRAB+SSXRD'] == 0) & (donnees['SSXRD'] == 1) & (donnees['boolZF'] == 1) & (donnees['Is_Pseudogene'] == 0) & (donnees['PRDM9'] == True)).astype(int)
donnees['SET+ZF'] = ((donnees['Complete_PRDM9'] == 0) & (donnees['SSXRD'] == 0) & (donnees['KRAB'] == 0) & (donnees['boolZF'] == 1) & (donnees['Is_Pseudogene'] == 0) & (donnees['PRDM9'] == True)).astype(int)
donnees['Pseudogene'] = donnees['Is_Pseudogene'][donnees['PRDM9'] == True]
donnees['ZF_Truncated'] = donnees['ZF_Truncated'][donnees['PRDM9'] == True]


#prdm9 = donnees[donnees['PRDM9'] == True]

synth = donnees.groupby(['Taxid', 'Species_name']).agg({
    'Complete_PRDM9': lambda x: list(donnees.loc[x.index, 'Protein_ID'][donnees['Complete_PRDM9'] == 1]),
    'Incomplete_PRDM9': lambda x: list(donnees.loc[x.index, 'Protein_ID'][donnees['Incomplete_PRDM9'] == 1]),
    'SET+KRAB+SSXRD': lambda x: list(donnees.loc[x.index, 'Protein_ID'][donnees['SET+KRAB+SSXRD'] == 1]),
    'SET+KRAB+ZF': lambda x: list(donnees.loc[x.index, 'Protein_ID'][donnees['SET+KRAB+ZF'] == 1]),
    'SET+KRAB': lambda x: list(donnees.loc[x.index, 'Protein_ID'][donnees['SET+KRAB'] == 1]),
    'KRAB+ZF': lambda x: list(donnees.loc[x.index, 'Protein_ID'][donnees['KRAB+ZF'] == 1]),
    'SET+SSXRD+ZF': lambda x: list(donnees.loc[x.index, 'Protein_ID'][donnees['SET+SSXRD+ZF'] == 1]),
    'SET+SSXRD': lambda x: list(donnees.loc[x.index, 'Protein_ID'][donnees['SET+SSXRD'] == 1]),
    'SSXRD+ZF': lambda x: list(donnees.loc[x.index, 'Protein_ID'][donnees['SSXRD+ZF'] == 1]),
    'SET+ZF': lambda x: list(donnees.loc[x.index, 'Protein_ID'][donnees['SET+ZF'] == 1]),
    'Pseudogene': lambda x: list(donnees.loc[x.index, 'Protein_ID'][donnees['Is_Pseudogene'] == 1]),
    'ZF_Truncated': lambda x: list(donnees.loc[x.index, 'Protein_ID'][donnees['ZF_Truncated'] == 1]),
    'Status': lambda x: list(donnees.loc[x.index, 'Protein_ID'][donnees['Status'] == 'Error during genewisedb'])


}).reset_index()

synth['Status'] = synth['Status'].apply(len)
for i in range(len(synth['Status'])):
    if synth['Status'][i] > 0:
        synth['Status'][i] = 'Error during genewisedb'
    else:
        synth['Status'][i] = 'Valid'
synth['Complete_PRDM9 nb'] = synth['Complete_PRDM9'].apply(len)
synth['SET+KRAB+SSXRD nb'] = synth['SET+KRAB+SSXRD'].apply(len)
synth['SET+KRAB+ZF nb'] = synth['SET+KRAB+ZF'].apply(len)
synth['SET+KRAB nb'] = synth['SET+KRAB'].apply(len)
synth['KRAB+ZF nb'] = synth['KRAB+ZF'].apply(len)
synth['SET+SSXRD+ZF nb'] = synth['SET+SSXRD+ZF'].apply(len)
synth['SET+SSXRD nb'] = synth['SET+SSXRD'].apply(len)
synth['SSXRD+ZF nb'] = synth['SSXRD+ZF'].apply(len)
synth['SET+ZF nb'] = synth['SET+ZF'].apply(len)
synth['Pseudogene nb'] = synth['Pseudogene'].apply(len)
synth['Incomplete_PRDM9 nb'] = (synth['SET+KRAB+SSXRD nb'] + synth['SET+KRAB+ZF nb'] + synth['SET+KRAB nb'] + synth['KRAB+ZF nb'] 
                                + synth['SET+SSXRD+ZF nb'] + synth['SET+SSXRD nb'] + synth['SSXRD+ZF nb'] + synth['SET+ZF nb'])

synth['ZF_Truncated nb'] = synth['ZF_Truncated'].apply(len)
nombre_lignes_SET = donnees[donnees['SET'] == 1].groupby(['Taxid', 'Species_name']).size().reset_index(name='Nb_SET')
synth = synth.merge(nombre_lignes_SET, on=['Taxid', 'Species_name'], how='left')
synth['Nb_SET'] = synth['Nb_SET'].fillna(0).astype(int)
synth['Total_PRDM9'] = (synth['Complete_PRDM9 nb'] + synth['Incomplete_PRDM9 nb'] + synth['Pseudogene nb'])

merge = pd.merge(synth, taxonomy, left_on='Taxid', right_on='Taxid', how='inner')
merge = merge.loc[:, ~merge.columns.str.contains('^Unnamed')]
merge = merge[['Taxid', 'Accession', 'Species_name', 'Status',
               'Nb_SET', 'Total_PRDM9', 'Complete_PRDM9 nb', 'Pseudogene nb', 'Incomplete_PRDM9 nb', 'ZF_Truncated nb',
               'Complete_PRDM9', 'Pseudogene', 'Incomplete_PRDM9', 'ZF_Truncated',
               'SET+KRAB+SSXRD nb', 'SET+KRAB+SSXRD', 'SET+KRAB+ZF nb', 'SET+KRAB+ZF', 'SET+KRAB nb', 'SET+KRAB',
               'SET+SSXRD+ZF nb', 'SET+SSXRD+ZF', 'SET+SSXRD nb', 'SET+SSXRD',
               'SSXRD+ZF nb', 'SSXRD+ZF', 'SET+ZF nb', 'SET+ZF']]
merge.to_csv(d["pathGTDriftGlobalResults"]+'prdm9_genomic_protein_analysis/summarized_results/table_prdm9.csv', sep = ';', index=False)
