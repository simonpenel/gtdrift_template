import pandas as pd
input_file = snakemake.input[0]
output_file = snakemake.output[0]

df = pd.read_csv(input_file, sep=';', header=0)
# remove first unamed column 
#df = df.drop(["Unnamed: 0"],axis=1) 
column_names = list(df.columns)
for name in column_names:           
    test = name.split(" ")
    if len(test) >= 2:
        if test[1] == "non-truncated" or  test[1] == "Intron" or  test[1] == "Stop/Frameshift" :
            to_reorder = df.pop(name)
            df[name] = to_reorder
    if name == "ProtRefID" or name == "Pseudogene (Genewise)" or name == "Pseudogene (HMMER)" or name == "Stop/Shift Positions" or name == "Nb Introns" or name == "ZF Truncated":
        to_reorder = df.pop(name)
        df[name] = to_reorder

df.drop(df.columns[df.columns.str.contains('unnamed', case=False)], axis=1, inplace=True)         
#df.to_csv(output_file, sep=';',index = False)
df.to_csv(output_file, sep=';')
