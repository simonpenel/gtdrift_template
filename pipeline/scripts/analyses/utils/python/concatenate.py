# =======================================
# Join several analysis outputs. 
# For example
#
#
#;SeqID;Assembly;KRAB Query;KRAB E-value;KRAB Score;Nb KRAB domains;KRAB domain start;KRAB domain end;Taxid;Species
#0;XP_009249362.3;GCF_028885655.2;Domain_KRAB_ReferenceAlignment2024;6.6e-29;101.8;1;29;102;9601;Pongo abelii
#1;XP_054411252.1;GCF_028885655.2;Domain_KRAB_ReferenceAlignment2024;8.2e-29;101.5;1;40;113;9601;Pongo abelii
#2;XP_024096381.2;GCF_028885655.2;Domain_KRAB_ReferenceAlignment2024;9.3e-16;59.7;1;26;86;9601;Pongo abelii
#3;XP_054400199.1;GCF_028885655.2;Domain_KRAB_ReferenceAlignment2024;2.9e-15;58.1;1;26;86;9601;Pongo abelii
#
# and
#;SeqID;Assembly;SET Query;SET E-value;SET Score;Nb SET domains;SET domain start;SET domain end;Taxid;Species;Best Match;Bit Score;Score ratio
#0;XP_054411252.1;GCF_028885655.2;Domain_SET_ReferenceAlignment2024;2.5000000000000002e-73;247.7;1;216;388;9601;Pongo abelii;Domain_SET_ReferenceAlignment2024;247.7;1.4301385681293304
#1;XP_009249362.3;GCF_028885655.2;Domain_SET_ReferenceAlignment2024;1.3e-58;199.7;1;205;351;9601;Pongo abelii;Domain_SET_ReferenceAlignment2024;199.7;1.3164139749505603
#2;XP_024111667.1;GCF_028885655.2;Domain_SET_ReferenceAlignment2024;4.7e-44;152.3;1;76;239;9601;Pongo abelii;SET_PRDM11_domain_alignment;362.8;2.3821405121470782
#3;XP_054380416.1;GCF_028885655.2;Domain_SET_ReferenceAlignment2024;4.7e-44;152.3;1;76;239;9601;Pongo abelii;SET_PRDM11_domain_alignment;362.8;2.3821405121470782
#
# =======================================

import pandas as pd
input_files = snakemake.input
output_file = snakemake.output[0]
accession = snakemake.params.accession

data_list = []
i = 0;
assembly = accession
taxid = "undefined"
species = "undefined"
for file in input_files:
    df = pd.read_csv(file, sep=';', header=0)
    # remove first unamed column 
    df = df.drop(["Unnamed: 0"],axis=1) 
    # get the domain name because we need to add it to the fields "Best Match" "Bit Score" and "Score ratio"
    domain = ""
    column_names = list(df.columns)
    for name in column_names:
        test = name.split(" ")
        if len(test) >= 2:
            if test[1] == "Query" :
                domain = test[0]
                break
    df = df.rename(columns={'Best Match': domain + ' Best Match','Bit Score': domain + ' Bit Score','Score ratio': domain + ' Score ratio'})
    print("Processing "+file) 
    if i == 0 :
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
        df = df.drop(["Taxid"],axis=1)
        df = df.drop(["Assembly"],axis=1)    
        df = df.drop(["Species"],axis=1)              
        # Put the types of the new dataframe in a dict
        orig = df.dtypes.to_dict()
        # Update the types of the current dataframe in the dict
        orig.update(df_cont.dtypes.to_dict())
        # Join the dataframes according to SeqId
        newdf = df_cont.join(df.set_index('SeqID'),on='SeqID', how="outer",lsuffix='_caller', rsuffix='_other')
        # Reset index 
        newdf = newdf.reset_index(drop=True)
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
    i = i + 1
# Moving Taxid and Species columns to the end    
column_taxid = df_cont.pop("Taxid")   
column_species = df_cont.pop("Species")   
df_cont['Taxid']=column_taxid
df_cont['Species']=column_species
# Write output
df_cont.to_csv(output_file, sep=';')
