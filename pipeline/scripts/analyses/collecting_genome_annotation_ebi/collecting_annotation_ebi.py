import subprocess, os, time
from html.parser import HTMLParser
#from bs4 import BeautifulSoup
import requests, sys
import numpy as np
import pandas as pd
import json

# URL to download
url_index = "https://ftp.ensembl.org/pub/" 
server = "https://rest.ensembl.org"
ext = "/info/genomes/taxonomy/Primates?"
version = "115"

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
    assembly_accession = line["assembly_accession"]
    taxid = line["species_taxonomy_id"]
    print(scientific_name + "\t" + str(taxid) + "\t" + assembly_name + "\t" + assembly_accession)

    url_fasta_genome = url_index + "release-" + version + "/" + "fasta/" + line["name"] + "/dna/" + line["url_name"] + "." +  line["assembly_default"] + ".dna.toplevel.fa.gz"
    print(url_fasta_genome)
    url_fasta_cds = url_index + "release-" + version + "/" + "fasta/" + line["name"] + "/cds/" + line["url_name"] + "." +  line["assembly_default"] + ".cds.all.fa.gz"
    print(url_fasta_cds)
    url_fasta_pep = url_index + "release-" + version + "/" + "fasta/" + line["name"] + "/pep/" + line["url_name"] + "." +  line["assembly_default"] + ".pep.all.fa.gz"
    print(url_fasta_pep)    
    url_gff3 = url_index + "release-" + version + "/" + "gff3/" + line["name"] + "/" + line["url_name"] + "." +  line["assembly_default"] + "." + version + ".gff3.gz"
    print(url_gff3)
    dico[index] = [scientific_name,taxid,assembly_name,assembly_accession,url_fasta_cds,url_fasta_pep,url_fasta_genome,] 
    accessions.append(assembly_accession)
    index +=1

sum = pd.DataFrame(dico).T
sum.columns = ['Species Name','Taxid','Assembly Name','Assembly Accession',"URL fasta cds","URL fasta proteins","URL fasta genome"]
print(sum)
sum.to_csv("organisms_data_ebi",sep="\t", index=False )

data = {
    "assembly_list" : accessions,
    "storagetype": "local"
}

json_str = json.dumps(data, indent=4)
with open("sample.json", "w") as f:
    f.write(json_str)
