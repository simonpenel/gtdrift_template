from ete3 import NCBITaxa
import pandas as pd
import os


"""
This script creates a complete taxonomy table using the ete3 python library
"""

def move_last_col(df):
    last_col_name = df.columns[-1]
    last_col = df.pop(last_col_name)
    df.insert(0, last_col_name, last_col)


os.system("""awk '/^>/ {sub(">", "", $1); print $1}' results/BLASTP_results/blastp_summary.txt > data/resources/taxid.txt""")

ncbi = NCBITaxa()

with open('data/resources/taxid.txt') as reader:
    data_list = [elt for elt in reader.readlines()]

tax_data = []
for elt in data_list:
    try:
        lineage = ncbi.get_lineage(elt)
    except ValueError as err:
        print("Error : ",format(err))
        lineage = []
    names = ncbi.get_taxid_translator(lineage)
    organism_data = [int(elt.strip())] + [names[taxid] for taxid in lineage]
    tax_data.append(organism_data)

taille_max = max(tax_data, key=len)
taille_max_value = len(taille_max)
columns = ['Taxid'] + [str(i) for i in range(1, taille_max_value)]

sorted_tax = pd.DataFrame(tax_data, columns=columns)
sorted_tax = sorted_tax.sort_values(by=columns)
for index, row in sorted_tax.iterrows():
    last_valid_col = row.last_valid_index()
    sorted_tax.at[index, 'Species name'] = sorted_tax.loc[index, last_valid_col]

move_last_col(sorted_tax)

with open('data/resources/organisms_data') as reader:
    dico = {}
    for line in reader.readlines()[1:]:
        taxid = int(line.split('\t')[1])
        assembly = line.split('\t')[2]
        dico[taxid] = assembly

for index, row in sorted_tax.iterrows():
    taxid = row['Taxid']
    sorted_tax.at[index, 'Accession'] = dico[taxid]

move_last_col(sorted_tax)
sorted_tax.to_csv('data/resources/sorted_taxonomy.csv', sep = ';')
