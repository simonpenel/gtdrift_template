# =======================================
# Join several analysis outputs. 
# For example

#
# =======================================

import pandas as pd
import numpy as np
input_files = snakemake.input
output_file = snakemake.output[0]
accession = snakemake.params.accession

data_list = []
dico_chromo = {} # info sur la position de la sequence dans l'assemblage
i = 0;
assembly = accession
taxid = "undefined"
species = "undefined"
for file in input_files:
    df = pd.read_csv(file, sep=';', header=0,dtype={"Stop/Shift Positions": str})
    #convert_dict = {"Stop/Shift Positions": int}
    #df = df.astype(convert_dict)
    #df["Stop/Shift Positions"] = df["Stop/Shift Positions"].replace(np.nan, 0)
    #df["Stop/Shift Positions"] = df["Stop/Shift Positions"].astype(str)
    # remove first unamed column 
    print("Process read file "+file)
    df["Stop/Shift Positions"] = df["Stop/Shift Positions"].astype(str)
    for index in range(0, len(df.index)):
        print(df["SeqID"][index])
        seqid = df["SeqID"][index]
        chromo =  df["Chromosome"][index]
        chromo_start =  df["Chr Start"][index]
        chromo_end =  df["Chr End"][index]    
        strand = df["Strand"][index]
        prot_len = df["Protein Length"][index]
        prot_id = df["ProtRefID"][index]    
        pseudo = df["Pseudogene (Genewise)"][index] 
        stop_shift = str(df["Stop/Shift Positions"][index])
        stop_shift = df["Stop/Shift Positions"][index]
        if seqid in dico_chromo :
            if   dico_chromo[seqid] != (chromo, chromo_start, chromo_end,strand,prot_len,prot_id,pseudo,stop_shift) :
                print("ERROR:")
                print("Current chromosome info for "+seqid+" :")
                print(dico_chromo[seqid])
                print("New chromosome info for "+seqid+" :")
                print((chromo, chromo_start, chromo_end,strand,prot_len,prot_id,pseudo,stop_shift))
                
                sys.exit('Multiple Chromosome')
        else :
            dico_chromo[seqid] = (chromo, chromo_start, chromo_end,strand,prot_len,prot_id,pseudo,stop_shift) 

for file in input_files:
    df = pd.read_csv(file, sep=';', header=0, dtype={"Stop/Shift Positions": str})
    # remove first unamed column 
    df = df.drop(["Unnamed: 0"],axis=1)
    # get domain name because we need to add it to the fields "Best Match" "Bit Score" and "Score ratio"     
    domain = ""
    column_names = list(df.columns)
    for name in column_names:
        test = name.split(" ")
        if len(test) >= 2:
            if test[1] == "Query" :
                domain = test[0]
                break
    df = df.rename(columns={'Best Match': domain + ' Best Match','Bit Score': domain + ' Bit Score','Score ratio': domain + ' Score ratio'})
    df["Stop/Shift Positions"] = df["Stop/Shift Positions"].replace(np.nan, 0)
    df["Stop/Shift Positions"] = df["Stop/Shift Positions"].astype(str)
    print("Processing "+file) 
    if i == 0 :
        df = df.drop(["Chromosome"],axis=1)
        df = df.drop(["Chr Start"],axis=1)
        df = df.drop(["Chr End"],axis=1)  
        df = df.drop(["Strand"],axis=1)
        df = df.drop(["Protein Length"],axis=1) 
        df = df.drop(["ProtRefID"],axis=1) 
        df = df.drop(["Pseudogene (Genewise)"],axis=1) 
        df = df.drop(["Stop/Shift Positions"],axis=1)   
        # First file
        df_cont = df
        if not df.empty:
            taxid = df_cont["Taxid"][0]
            #assembly = df_cont["Assembly"][0]
            species = df_cont["Species"][0]
    else :
        # Folowing files
        # Remove redundant Taxid Assembly and Species fields
        if not df.empty :
            if taxid == "undefined":
                taxid = df["Taxid"][0]
                #assembly = df["Assembly"][0]
                species = df["Species"][0]
            
            if df["Taxid"][0] != taxid :
                 sys.exit('Multiple Taxids, Species ot Assembly in the input files.')
            if df["Assembly"][0] != assembly :
                 sys.exit('Assembly in the file is different from the parameter')    
            if df["Species"][0] != species :
                 sys.exit('Multiple Taxids, Species ot Assembly in the input files.')                         
        orig = df.dtypes.to_dict()
        df = df.drop(["Taxid"],axis=1)
        df = df.drop(["Assembly"],axis=1)    
        df = df.drop(["Species"],axis=1) 
        df = df.drop(["Chromosome"],axis=1)
        df = df.drop(["Chr Start"],axis=1)
        df = df.drop(["Chr End"],axis=1)  
        df = df.drop(["Strand"],axis=1)
        df = df.drop(["Protein Length"],axis=1) 
        df = df.drop(["ProtRefID"],axis=1) 
        df = df.drop(["Pseudogene (Genewise)"],axis=1) 
        df = df.drop(["Stop/Shift Positions"],axis=1)             
        # Put the types of the new dataframe in a dict
        #orig = df.dtypes.to_dict()
        # Update the types of the current dataframe in the dict
        orig.update(df_cont.dtypes.to_dict())
        #df_cont["Stop/Shift Positions"] = df_cont["Stop/Shift Positions"].astype(str)
        # Join the dataframes according to SeqId
        newdf = df_cont.join(df.set_index('SeqID'),on='SeqID', how="outer",lsuffix='_caller', rsuffix='_other')
        # Reset index 
        newdf = newdf.reset_index(drop=True)
        #newdf["Stop/Shift Positions"] = newdf["Stop/Shift Positions"].astype(str)
        # Set missing data to 0
        #if flag_multiple : 
        #    sys.exit('Multiple Taxids, Species ot Assembly in the input files.')
        # Set missing taxid,species and assembly
        newdf["Taxid"] = newdf["Taxid"].fillna(value=taxid)
        newdf["Assembly"] = newdf["Assembly"].fillna(value=assembly)
        newdf["Species"] = newdf["Species"].fillna(value=species)    
        # Set other missing data to 0
        df_cont = newdf.fillna(0)
        # Apply types   
        df_cont = df_cont.apply(lambda x: x.astype(orig[x.name]))
        #df_cont["Stop/Shift Positions"] = df_cont["Stop/Shift Positions"].astype(str)
    i = i + 1


df_chromo = pd.DataFrame() 
for key in  dico_chromo :
    values = dico_chromo[key]
    new_row = {"SeqID": key, "Chromosome": values[0], "Chr Start": values[1],"Chr End": values[2],
    "Strand" : values[3], "Protein Length" : values[4], "ProtRefID" : values[5],
    "Pseudogene (Genewise)" : values[6], "Stop/Shift Positions" : values[7]}
    print(new_row)
    df_chromo = df_chromo._append(new_row, ignore_index=True)
#print(df_chromo.SeqID)
#print(df_chromo.Chromosome)
print("debug "+ str(len(df_chromo)))
# esquive le cas ou il n'y a pas de candidat
if len(df_chromo) > 0 :
    newdf = df_cont.join(df_chromo.set_index('SeqID'),on='SeqID', how="outer", lsuffix='_caller', rsuffix='_other')
    df_cont= newdf



# Moving Taxid and Species columns to the end    
column_taxid = df_cont.pop("Taxid")   
column_species = df_cont.pop("Species")   
df_cont['Taxid']=column_taxid
df_cont['Species']=column_species
# Write output
df_cont.drop(df_cont.columns[df_cont.columns.str.contains('unnamed', case=False)], axis=1, inplace=True)
df_cont.to_csv(output_file, sep=';')
#df_cont.to_csv("debug.csv", sep=';')
