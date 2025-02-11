# =============================================
# Read analysis output of each domain and 
# select the candidates with all domains (wad). 
# ==============================================

import pandas as pd

input_files_conf = snakemake.input.candidates_confirmed

# -------------------------------------------------------------
# Format:
#SeqID;Assembly;SET Query;SET E-value;SET Score;Nb SET domains;SET domain start;SET domain end;Taxid;Species;Best Match;Bit Score;Score ratio
#XP_030859242.3;GCF_029281585.2;Domain_SET_ReferenceAlignment2024;9.4e-73;245.8;1;205;377;9595;Gorilla gorilla;Domain_SET_ReferenceAlignment2024;245.8;1.4241019698725377
#XP_055223462.2;GCF_029281585.2;Domain_SET_ReferenceAlignment2024;1.0000000000000001e-72;245.7;1;244;416;9595;Gorilla gorilla;Domain_SET_ReferenceAlignment2024;245.7;1.4243478260869564
#XP_030858898.2;GCF_029281585.2;Domain_SET_ReferenceAlignment2024;8.8e-68;229.6;1;205;363;9595;Gorilla gorilla;Domain_SET_ReferenceAlignment2024;229.6;1.359384251036116
# ...
# -------------------------------------------------------------

input_files_simple = snakemake.input.candidates_simple

# -------------
# Format:
#XP_030859242.3
#XP_055223462.2
#XP_030858898.2
#XP_030861841.2
#XP_030861853.1
#...
# -------------

output_file = snakemake.output[0]

# ----------------------------------------------------
# Format:
#;SeqID;Assembly;Taxid;Species
#0;XP_030859242.3;GCF_029281585.2;9595;Gorilla gorilla
#1;XP_055223462.2;GCF_029281585.2;9595;Gorilla gorilla
#...
# ----------------------------------------------------

candidate_list = []
i = 0
# First, read the results for confirmed domains
for file in input_files_conf:
    df = pd.read_csv(file, sep=';', header=0)
    df=df.dropna(how='any') 
    # If there is no data for the domain, create an output file with no data and exit 
    if df.empty :
        data = []
        candidate_list = pd.DataFrame(data, columns=["SeqID","Assembly","Taxid","Species"])
        candidate_list.to_csv(output_file, sep=';')
        exit()
    df = df[["SeqID","Assembly","Taxid","Species"]]    
    if i == 0 :
        # Initialise the output
        candidate_list = df
    else :
        # join the new data to the current output
        candidate_list = candidate_list.join(df.set_index('SeqID'),on='SeqID', how="inner",lsuffix='_caller', rsuffix='_other')
   
    i = i + 1
    
# Then, read the results for the rest of  domains     
for file in input_files_simple:
    df = pd.read_csv(file, sep=';', header=None,names=["SeqID"])
    df=df.dropna(how='any')
    # If there is no data for the domain, create an output file with no data and exit 
    if df.empty :
        data = []
        candidate_list = pd.DataFrame(data, columns=["SeqID","Assembly","Taxid","Species"])
        candidate_list.to_csv(output_file, sep=';')
        exit()
    # If there is no data for the current output, create an output file with no data and exit 
    if candidate_list.empty :
        data = []
        candidate_list = pd.DataFrame(data, columns=["SeqID","Assembly","Taxid","Species"])
        candidate_list.to_csv(output_file, sep=';')
        exit()
    # join the new data to the current output    
    candidate_list = candidate_list.join(df.set_index('SeqID'),on='SeqID', how="inner",lsuffix='_caller', rsuffix='_other')
    
# Write output
candidate_list.to_csv(output_file, sep=';')
