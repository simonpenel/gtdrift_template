import subprocess, os, time
from html.parser import HTMLParser
from bs4 import BeautifulSoup
import requests, sys



# URL to download
url_index = "https://ftp.ensembl.org/pub/" 
server = "https://rest.ensembl.org"
ext = "/info/genomes/taxonomy/Homo sapiens?"
version = "115"

r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
 
if not r.ok:
  r.raise_for_status()
  sys.exit()
 
decoded = r.json()
for line in decoded:
    print(line)
    url_fasta = url_index + "release-" + version + "/" + "fasta/" + line["name"] + "/dna/" + line["url_name"] + "." +  line["assembly_default"] + ".dna.toplevel.fa.gz"
    print(url_fasta)
    url_gff3 = url_index + "release-" + version + "/" + "gff3/" + line["name"] + "/" + line["url_name"] + "." +  line["assembly_default"] + "." + version + ".gff3.gz"
    print(url_gff3)


