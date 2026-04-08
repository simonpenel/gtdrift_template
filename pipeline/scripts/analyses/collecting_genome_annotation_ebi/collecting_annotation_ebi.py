import subprocess, os, time
from html.parser import HTMLParser
#from bs4 import BeautifulSoup
import requests, sys
import numpy as np
import pandas as pd
import json
# URL to download
release_vertebrata = "115"
release_metazoa = "62"
url_index_vertebrata = "https://ftp.ensembl.org/pub/release-" + release_vertebrata + "/"
url_index_metazoa = "https://ftp.ensemblgenomes.ebi.ac.uk/pub/release-" + release_metazoa + "/metazoa/"

# REST request
server = "https://rest.ensembl.org"
ext = "/info/genomes/taxonomy/Primates?"
ext = "/info/genomes/taxonomy/Chordata?"
ext = "/info/genomes/taxonomy/Metazoa?"
#ext = "/info/genomes/taxonomy/Ovis?"
r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
 
if not r.ok:
  r.raise_for_status()
  sys.exit()
 
decoded = r.json()
dico = {}
index = 0
accessions = []
for line in decoded:
    scientific_name = line["scientific_name"]
    assembly_name = line["assembly_name"]
    strain = line["strain"]
    assembly_accession = line["assembly_accession"]
    taxid = line["species_taxonomy_id"]
    division = line["division"]
    print(scientific_name + "\t" + str(taxid) + "\t" + assembly_name + "\t" + str(assembly_accession))
    if division == "EnsemblVertebrates":
        url_index = url_index_vertebrata
        version = release_vertebrata
    elif division == "EnsemblMetazoa":
        url_index = url_index_metazoa
        version = release_metazoa
    else:
        sys.exit("Wrong division")

    url_fasta_genome = url_index  + "fasta/" + line["name"] + "/dna/" + line["url_name"] + "." +  line["assembly_default"] + ".dna.toplevel.fa.gz"
    print(url_fasta_genome)
    url_fasta_cds = url_index  + "fasta/" + line["name"] + "/cds/" + line["url_name"] + "." +  line["assembly_default"] + ".cds.all.fa.gz"
    print(url_fasta_cds)
    url_fasta_pep = url_index + "fasta/" + line["name"] + "/pep/" + line["url_name"] + "." +  line["assembly_default"] + ".pep.all.fa.gz"
    print(url_fasta_pep)    
    url_gff3 = url_index + "gff3/" + line["name"] + "/" + line["url_name"] + "." +  line["assembly_default"] + "." + version + ".gff3.gz"
    print(url_gff3)
    dico[index] = [scientific_name,strain,taxid,assembly_name,assembly_accession,division,url_fasta_cds,url_fasta_pep,url_fasta_genome,url_gff3] 
    accessions.append(assembly_accession)
    index +=1

sum = pd.DataFrame(dico).T
sum.columns = ['Species Name','Strain','Taxid','Assembly Name','Assembly Accession',"Division","URL fasta cds","URL fasta proteins","URL fasta genome","URL gff3 genome"]
print(sum)
sum.to_csv("organisms_data_ebi",sep="\t", index=False )
ref = sum[(sum['Strain'] == 'reference')]
ref2 = sum[(sum['Strain'].isnull())]
ref = pd.concat([ref,ref2]).reset_index(drop=True)
ref.to_csv("organisms_data_ebi.ref_only",sep="\t", index=False )


data = {
    "assembly_list" : accessions,
    "storagetype": "local"
}

json_str = json.dumps(data, indent=4)
with open("sample.json", "w") as f:
    f.write(json_str)
