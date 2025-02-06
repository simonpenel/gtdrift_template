# =======================================
# Join several domain stats outputs. 
# For example
#

import pandas as pd
input_files = snakemake.input.all_domains_stats
output_file = snakemake.output[0]
flag_multiple = False

i = 0;
for file in input_files:
    df = pd.read_csv(file, sep=';', header=0)
    df = df.drop(["Unnamed: 0"],axis=1) 
    print(df)
    if i ==  0 :
        df_cont = df
        print(df_cont)
    else :
        #df_cont = pd.merge(df_cont,df)
        df_cont = pd.concat([df_cont, df], ignore_index=True)
        print(df_cont)
    i += 1


# Write output
df_cont.to_csv(output_file, sep=';')
