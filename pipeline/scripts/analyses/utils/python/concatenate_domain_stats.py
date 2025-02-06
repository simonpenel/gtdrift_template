# =======================================
# Join several domain stats outputs. 
# For example
#

import pandas as pd
input_files = snakemake.input.domains
input_files_simple = snakemake.input.domains_simple
input_files += input_files_simple
organisms_file = snakemake.input.organisms_file
output_file = snakemake.output[0]
flag_multiple = False

i = 0;
for file in input_files:
    df = pd.read_csv(file, sep=';', header=0)
    df = df.drop(["Unnamed: 0"],axis=1) 
    if i == 0 :
        # First file
        df_cont = df
        nb_prot = df_cont["Nb proteins"][0]  
        assembly =  df_cont["Assembly"][0]    
    else :
        # Folowing files
        # Remove redundant Taxid Assembly and Species fields
        if df["Nb proteins"][0] != nb_prot :
            flag_multiple = True                        
        df = df.drop(["Nb proteins"],axis=1)              
        # Join the dataframes according to Assemblmy
        newdf = df_cont.join(df.set_index('Assembly'),on='Assembly', how="outer",lsuffix='_caller', rsuffix='_other')
        # Reset index 
        newdf = newdf.reset_index(drop=True)
        df_cont = newdf
    i = i + 1

# Moving Nb proteins columns to the end    
column_nbprot = df_cont.pop("Nb proteins")   
df_cont['Nb proteins']=column_nbprot   
df_orga = pd.read_csv(organisms_file, sep='\t', header=0)
taxid = df_orga.loc[df_orga['Assembly Accession'] == assembly, 'Taxid'].values[0]    
df_cont["Taxid"] = taxid 
species = df_orga.loc[df_orga['Assembly Accession'] == assembly, 'Species Name'].values[0]    
df_cont["Species"] = species  
# Write output
df_cont.to_csv(output_file, sep=';')
