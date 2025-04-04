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
import re
import sys
import pandas as pd

res_file = snakemake.input.res
protein_file = snakemake.input.protein
output_file = snakemake.output[0]
domain = snakemake.params.domain
accession = snakemake.params.accession
restype = snakemake.params.type
if restype == "simple" :
    file=open(res_file, 'r')
    nb_hits = len(file.readlines())
    file.close()
elif restype == "confirmed" :
    input_df = pd.read_csv(res_file, sep=';', header=0)
    nb_hits = len(input_df.dropna(how='any')   )
    
else :
    sys.exit("Type inconnu")


nb_prot = 0
for line in open(protein_file, 'r') :
    if re.search(">", line) :
        nb_prot += 1
        
        
#df = pd.read_csv(res_file, sep=';', header=0)
print("RES FILES = ",res_file)
print("PROT FILE = ",protein_file)
print("OUTPUT FILE = ",output_file)
print("DOMAIN = ",domain)
print("Nb hits ",nb_hits)
print("Nb prot ",nb_prot)

if restype == "simple" :
    d = {'Assembly': accession,'Nb hits '+domain: [nb_hits], 'Nb proteins': [nb_prot]}
else :
    d = {'Assembly': accession,'Nb confirmed hits '+domain: [nb_hits], 'Nb proteins': [nb_prot]}

df = pd.DataFrame(data=d)
df.to_csv(output_file, sep=';')

#file=open(output_file, 'w')
#file.write(str(nb_hits))
#file.close()
