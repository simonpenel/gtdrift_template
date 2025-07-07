# =======================================
# Join several tyrosines outputs. 
# For example
#

import pandas as pd
input_files = snakemake.input
output_file = snakemake.output[0]
flag_multiple = False

i = 0;
for file in input_files:
    df = pd.read_csv(file, sep=';', header=0)
    if i ==  0 :
        df_cont = df
    else :
        df_cont = pd.concat([df_cont, df], ignore_index=True)
    i += 1

#df_cont = df_cont.drop(["0 candidate PRDM9 proteins"],axis=1) 
# Write output
df_cont.to_csv(output_file, sep=';')
