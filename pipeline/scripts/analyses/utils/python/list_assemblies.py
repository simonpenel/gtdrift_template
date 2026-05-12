# =======================================
# Join several  outputs. 

import pandas as pd
input_files = snakemake.input
output_file = snakemake.output[0]
# p=snakemake.params
# toto = snakemake.ACCESSNB
# print("paramtres " + str(toto))

i = 0;
for file in input_files:    
    df = pd.read_csv(file, sep='\t', header=0)
    if i ==  0 :
        df_cont = df
    else :
        df_cont = pd.concat([df_cont, df], ignore_index=True)
    i += 1

# #df_cont = df_cont.fillna(0.0)  
# df_cont = df_cont.fillna(0)  
# # Write output
# # Moving Taxid and Species columns to the end    
# column_taxid = df_cont.pop("Taxid")   
# column_species = df_cont.pop("Species")   
# df_cont['Taxid']=column_taxid
# df_cont['Species']=column_species

# df_cont.drop(df_cont.columns[df_cont.columns.str.contains('unnamed', case=False)], axis=1, inplace=True)
df_cont.to_csv(output_file, sep='\t',index=False)
