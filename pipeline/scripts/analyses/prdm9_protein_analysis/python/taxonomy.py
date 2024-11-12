from ete3 import NCBITaxa
import pandas as pd
import os

import argparse

parser = argparse.ArgumentParser(description='Creates a complete taxonomy table using the ete3 python library')

parser.add_argument('-i', '--input_dir', type=str, required=True, help='Input dir path')
parser.add_argument('-o', '--output_dir', type=str, required=True, help='Output dir path')

args = parser.parse_args()
input_dir = args.input_dir
output_dir = args.output_dir

"""
This script creates a complete taxonomy table using the ete3 python library
"""

def move_last_col(df):
    last_col_name = df.columns[-1]
    last_col = df.pop(last_col_name)
    df.insert(0, last_col_name, last_col)

print(f"Generate "+output_dir+"taxid.txt from "+output_dir+"analyses_summaries/BLASTP_results/blastp_summary.txt")
os.system("awk '/^>/ {sub(\">\", \"\", $1); print $1}' "+ output_dir+"analyses_summaries/BLASTP_results/blastp_summary.txt >  "+output_dir+"taxid.txt")

ncbi = NCBITaxa()

dicospec = {}
with open(f'{input_dir}/organisms_data') as reader:
    #dico = {}
    for line in reader.readlines()[1:]:
        taxid = int(line.split('\t')[1])
        species = line.split('\t')[0]
        dicospec[taxid] = species
                
with open(f'{output_dir}/taxid.txt') as reader:
    data_list = [elt for elt in reader.readlines()]

tax_data = []
for elt in data_list:
    try:
        lineage = ncbi.get_lineage(elt)
    except ValueError as err:
        print("Error  for taxonmy: ",format(err))
        lineage = []
    names = ncbi.get_taxid_translator(lineage)
    taxid = int(elt.strip())
    if taxid in dicospec:
        species = dicospec[taxid]
    else:
        species = "Not found"
        print("Taxid "+str(taxid)+" not found")
    
    #organism_data = [int(elt.strip())] + [names[taxid] for taxid in lineage]
    organism_data = [int(elt.strip())] + [names[taxid] for taxid in lineage] +[species]
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

#with open(f'{input_dir}/ncbi_genome_assembly_taxon.txt') as reader:
with open(f'{input_dir}/organisms_data') as reader:
    dico = {}
    for line in reader.readlines()[1:]:
        taxid = int(line.split('\t')[1])
       #taxid = int(line.split('\t')[6])
        assembly = line.split('\t')[2]
        dico[taxid] = assembly

for index, row in sorted_tax.iterrows():
    taxid = row['Taxid']
#    sorted_tax.at[index, 'Accession'] = dico[taxid]
    if taxid in dico :
        sorted_tax.at[index, 'Assembly Accession'] = dico[taxid]
    else :
        print("Error : taxid is not present "+ str(taxid))
    
move_last_col(sorted_tax)
sorted_tax.to_csv(f'{output_dir}/sorted_taxonomy.csv', sep = ';')
