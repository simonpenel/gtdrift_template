# =======================================
# Join several  outputs. 

import pandas as pd
input_files = snakemake.input
output_file = snakemake.output[0]

i = 0;
for file in input_files:    
    df = pd.read_csv(file, sep=';', header=0)
    if i ==  0 :
        df_cont = df
    else :
        df_cont = pd.concat([df_cont, df], ignore_index=True)
    i += 1

#df_cont = df_cont.fillna(0.0)  
#df_cont = df_cont.fillna(0)  
# Write output
# Moving Taxid and Species columns to the end    
column_taxid = df_cont.pop("Taxid")   
column_species = df_cont.pop("Species")   
df_cont['Taxid']=column_taxid
df_cont['Species']=column_species

df_cont.drop(df_cont.columns[df_cont.columns.str.contains('unnamed', case=False)], axis=1, inplace=True)

columns = list(df_cont.columns)

tobeintegers  = ["hits","domains", "Length","start","end","Taxid"] 
for column in columns:
    for test in tobeintegers:
        if test in column:
            print("Change "+column)
            df_cont[column] = df_cont[column].astype('Int64')

tobe0whenNA  = ["hits","domains"] 
for column in columns:
    for test in tobe0whenNA:
        if test in column:
            print("Na is 0 for  "+column)
            df_cont[column] = df_cont[column].fillna(value=0)

if 'Genewise index' in df_cont.columns:
    df_cont['Genewise index'] = df_cont['Genewise index'].astype('Int64')
    df_cont['Protein Length'] = df_cont['Protein Length'].astype('Int64')
    df_cont['Chr Start'] = df_cont['Chr Start'].astype('Int64')
    df_cont['Chr End'] = df_cont['Chr End'].astype('Int64')

df_cont.to_csv(output_file, sep=';',na_rep="NA",index=False)
