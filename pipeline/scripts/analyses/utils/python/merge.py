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

df_cont = df_cont.fillna(0.0)  
# Write output
df_cont.to_csv(output_file, sep=';')
