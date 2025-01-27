# =======================================
# Join several analysis outputs. 
# For example
#
#;SeqID;SET Query;SET E-value;SET Score;Nb SET domains;SET domain start;SET domain end;Taxid;Best Match;Bit Score;Score ratio
#0;NP_001363829.1;Domain_SET_ReferenceAlignment2024;1.2e-72;246.2;1;205;377;9606;Domain_SET_ReferenceAlignment2024;246.2;1.4338963308095516
#1;NP_064612.2;Domain_SET_ReferenceAlignment2024;1.2e-72;246.2;1;205;377;9606;Domain_SET_ReferenceAlignment2024;246.2;1.4338963308095516
#2;XP_016878373.1;Domain_SET_ReferenceAlignment2024;4e-71;241.3;1;26;190;9606;Domain_SET_ReferenceAlignment2024;241.3;1.4086398131932283
#3;XP_016878371.1;Domain_SET_ReferenceAlignment2024;7.5e-71;240.4;1;132;296;9606;Domain_SET_ReferenceAlignment2024;240.4;1.4116265413975337
#...
# and
#
#;SeqID;KRAB Query;KRAB E-value;KRAB Score;Nb KRAB domains;KRAB domain start;KRAB domain end;Taxid
#0;NP_001363829.1;Domain_KRAB_ReferenceAlignment2024;1.3e-28;101.6;1;29;102;9606
#1;NP_064612.2;Domain_KRAB_ReferenceAlignment2024;1.3e-28;101.6;1;29;102;9606
#2;XP_011521133.1;Domain_KRAB_ReferenceAlignment2024;2.6e-28;100.7;1;29;102;9606
#3;XP_011521131.1;Domain_KRAB_ReferenceAlignment2024;5.1e-28;99.7;1;29;102;9606
# ...
# =======================================

import pandas as pd
input_files = snakemake.input
output_file = snakemake.output[0]

data_list = []
i = 0;
flag_multiple_taxid = False
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
        taxid = df_cont["Taxid"][0]
    else :
        # Folowing files
        # Remove redundant Taxid field
        if df_cont["Taxid"][0] != taxid :
            flag_multiple_taxid = True
        df = df.drop(["Taxid"],axis=1)
        # Put the types of the new dataframe in a dict
        orig = df.dtypes.to_dict()
        # Update the types of the current dataframe in the dict
        orig.update(df_cont.dtypes.to_dict())
        # Join the dataframes according to SeqId
        newdf = df_cont.join(df.set_index('SeqID'),on='SeqID', how="outer",lsuffix='_caller', rsuffix='_other')
        # Set missing data to 0
        if flag_multiple_taxid : 
            sys.exit('Multiple Taxids in the input files.')
        # Set missing taxid
        newdf["Taxid"] = newdf["Taxid"].fillna(value=taxid)
        # Set other missing data to 0
        df_cont = newdf.fillna(0)
        # Apply types   
        df_cont = df_cont.apply(lambda x: x.astype(orig[x.name]))
    i = i + 1

df_cont.to_csv(output_file, sep=';')
