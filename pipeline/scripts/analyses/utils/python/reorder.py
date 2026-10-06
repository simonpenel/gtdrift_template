import pandas as pd
input_file = snakemake.input[0]
output_file = snakemake.output[0]

df = pd.read_csv(input_file, sep=';', header=0,dtype={'Stop/Shift Positions': "str"})
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
columns = list(df.columns)

tobeintegers  = ["hits","domains", "Length","start","end","Taxid","Start","End","index","non-truncated","Intron","Stop/Frameshift"] 
for column in columns:
    for test in tobeintegers:
        if test in column:
            print("Change "+column)
            df[column] = df[column].astype('Int64')

tobe0whenNA  = ["hits","domains"] 
for column in columns:
    for test in tobe0whenNA:
        if test in column:
            print("Na is 0 for  "+column)
            df[column] = df[column].fillna(value=0)

toberemoved  = ["Genewise index","Pseudogene (HMMER)","Pseudogene (Genewise)"]


toberemovedmatch  = ["non-truncated"]
for column in columns:
    for test in toberemovedmatch:
        if test in column:
            toberemoved.append(column)

print("Remove "+str(toberemoved))
df = df.drop(columns=toberemoved)
df.to_csv(output_file, sep=';',index = False,na_rep="NA")
