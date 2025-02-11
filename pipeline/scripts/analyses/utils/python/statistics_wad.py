# ========================================
# Create a file describing for a genome
# the number of candidates with all
# domains (wad) and the number of proteins.
# The Taxid and the Species are added too.
# =========================================
import re
import sys
import pandas as pd

res_file = snakemake.input.res
# ---------------------------------------------------
# Format:
#;SeqID;Assembly;Taxid;Species
#0;XP_030859242.3;GCF_029281585.2;9595;Gorilla gorilla
#1;XP_055223462.2;GCF_029281585.2;9595;Gorilla gorilla
# ...
# ---------------------------------------------------

protein_file = snakemake.input.protein
# ---------------------------------------------------
# Fasta file
# ---------------------------------------------------

organisms_file = snakemake.input.organisms_file
# ---------------------------------------------------
# Format:
#Species Name	Taxid	Assembly Accession	Existing Annotation	Existing Protein Sequence	URL
#Grapholita funebrana	568174	GCA_038095595.1	False	False	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/038/095/595/GCA_038095595.1_ASM3809559v1/GCA_038095595.1_ASM3809559v1
#Grapholita dimorpha	934274	GCA_038095585.1	False	False	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/038/095/585/GCA_038095585.1_ASM3809558v1/GCA_038095585.1_ASM3809558v1
#Lyrurus tetrix	1233216	GCA_043882375.1	False	False	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/043/882/375/GCA_043882375.1_ASM4388237v1/GCA_043882375.1_ASM4388237v1
#...
# ---------------------------------------------------

output_file = snakemake.output[0]
# Example:
#;Assembly;Nb complete candidates;Candidates;Nb proteins;Taxid;Species
#0;GCF_029281585.2;2;['XP_030859242.3', 'XP_055223462.2'];80997;9595;Gorilla gorilla
# (only 1 ligne because it is the results for 1 genome)

accession = snakemake.params.accession

# Read the  file
input_df = pd.read_csv(res_file, sep=';', header=0)
# Remove lines with NaN
input_df = input_df.dropna(how='any')
nb_hits = len(input_df)

candidates = input_df["SeqID"].tolist()

nb_prot = 0
for line in open(protein_file, 'r') :
    if re.search(">", line) :
        nb_prot += 1
        
d = {'Assembly': accession,'Nb complete candidates': [nb_hits], 'Candidates':[candidates],'Nb proteins': [nb_prot]}

df = pd.DataFrame(data=d)

# We can not use the info in the res_file because it may happen that the file is empty, then
# there is no info on the species and taxid. However even there is not results, we need to
# generate an output, thus we need this info. So we use the accession and the organisms_data file.

df_orga = pd.read_csv(organisms_file, sep='\t', header=0)
taxid = df_orga.loc[df_orga['Assembly Accession'] == accession, 'Taxid'].values[0]    
df["Taxid"] = taxid 
species = df_orga.loc[df_orga['Assembly Accession'] == accession, 'Species Name'].values[0]    
df["Species"] = species  
df.to_csv(output_file, sep=';')

