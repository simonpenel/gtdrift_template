import sys
import os
import pandas as pd

"""
This script extracts the sequence selected by hmm search for an organism and runs a blastp analysis against the Human PRDM genes family. If the best match is PRDM9, 
the value is saved and compared to the next best non-PRDM9 match.
The ouput file is named blastp_summary.txt and contains the taxid, the best PRDM match, the presence/absence data for every proteic domain, the bit score of the blastp
if the best match is PRDM9 and the ratio with the second best non-PRDM9 match.
"""
accession_number = snakemake.params.accession
summary_file = snakemake.input.summary
blastdb = snakemake.input.blastdb
prdmdb = snakemake.input.prdmdb

blastp_file = snakemake.output.blastp_file
SET_sequences_dir = snakemake.output.SET_sequences_dir
SET_blastp_dir = snakemake.output.SET_blastp_dir

df = pd.read_csv(summary_file, sep=';')
outputfile = snakemake.output.table

with open(blastp_file, 'w') as writer:
    taxid = None
    string = ''
    #os.system(f"mkdir -p {inputdir}/{accession}/analyses/prdm9_prot/SET_sequences/")
    #os.system(f"mkdir -p {inputdir}/{accession}/analyses/prdm9_prot/SET_blastp/")
    #print("mkdir -p "+SET_sequences_dir)
    os.system("mkdir -p "+SET_sequences_dir)
    #print("mkdir -p "+SET_blastp_dir)  
    os.system("mkdir -p "+SET_blastp_dir)  
    for index, row in df.iterrows():
        taxid = f">{df['Taxid'].iloc[0]}\n"
        set = 0
        krab = 0
        ssxrd = 0
        zf = 0
        if row['Nb SET domains'] != 0: 
            set = 1
        if row['Nb KRAB domains'] != 0: 
            krab = 1
        if row['Nb SSXRD domains'] != 0:
            ssxrd = 1
        if row['Nb ZF domains'] != 0:
            zf = row['Nb ZF domains']
        prot = f"<\t{set}\t{krab}\t{ssxrd}\t{zf}\n"
        #print(f"Run blastdbcmd -db {blastdb}/protdb -entry {row['SeqID']} -range {int(row['SET domain start'])}-{int(row['SET domain end'])} -out {SET_sequences_dir}/{row['SeqID']}.fa")
        if int(row['SET domain end']) > int(row['SET domain start']) :
            ret = os.system(f"blastdbcmd -db {blastdb}/protdb -entry {row['SeqID']} -range {int(row['SET domain start'])}-{int(row['SET domain end'])} -out  {SET_sequences_dir}/{row['SeqID']}.fa")            
            if ret > 0 :
                sys.exit("Error during blastdbcmd")
                
            ret = os.system(f"blastp -db {prdmdb}/prdm_family -outfmt 7 -query  {SET_sequences_dir}/{row['SeqID']}.fa -out  {SET_blastp_dir}/{row['SeqID']}")            
            if ret > 0 :
                sys.exit("Error during blastp")

            with open(f"{SET_blastp_dir}/{row['SeqID']}") as reader:            
                prot_id = row['SeqID']
                lines = reader.readlines()
                prdm_match = lines[5].split()[1].split('_')[0]
                df.at[index, 'Best Match'] = lines[5].split()[1].split('_')[0] # Best match Prdm number
                if prdm_match == 'PRDM9': # if it is prdm9, save the score and compare it to the next non-prdm9 best match
                    df.at[index, 'Bit Score'] = float(lines[5].split()[-1])
                    j = 1
                    while lines[5 + j].split()[1].split('_')[0] == 'PRDM9': 
                        j += 1
                    if lines[5 + j].split(" ")[0] != "#":
                    	df.at[index, 'Score ratio'] = float(lines[5].split()[-1])/float(lines[5 + j].split()[-1])
                    	string += f"{prot}{prot_id}\t{float(lines[5].split()[-1])}\t{float(lines[5].split()[-1])/float(lines[5 + j].split()[-1])}\t{prdm_match}\n"
                    else:
                    	print(f"Only 1 hit for {SET_blastp_dir}/{row['SeqID']}")
                else:
                    string += f"{prot}{prot_id}\t{prdm_match}\n"
        else :
            print(f"blastdbcmd was not processed")        
    df.to_csv(f"{outputfile}", sep=';')
    if taxid == None:
        print("Nothing found in blastp")
    else:
        writer.write(taxid + string)  
