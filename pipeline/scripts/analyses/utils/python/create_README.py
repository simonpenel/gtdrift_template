import pandas as pd
import argparse
from math import isnan
import os
import sys
import json
import warnings
warnings.filterwarnings("error")
universal_col_names = ['Unnamed: 0', 'SeqID', 'Assembly', 'Taxid', 'Species']
dico_col_names_univ = {
       'SeqID': 'Sequence name from the proteome', 
       'Assembly':'Assembly accession number',
       'Taxid':'Taxonomic identifier',
       'Species':'Species name'
       }
col_names_domain_univ = dico_col_names_univ.keys()


dico_col_names_domain_simple = {
       'DOMAIN Query'       :'Sequence used to search the DOMAIN domain with hmmer',
       'DOMAIN E-value'     :'E-value of the DOMAIN domain (hmmer)', 
       'DOMAIN Score'       :'Score of the DOMAIN domain (hmmer)',
       'Nb DOMAIN domains'  :'Number of DOMAIN domains (hmmer) [if many, only the last domain is selected]',
       'DOMAIN domain start':'Start of the DOMAIN domain',
       'DOMAIN domain end'  :'Start of the DOMAIN domain',
       'DOMAIN HMM cov.'    :'Ratio of the hmm part covering the protein over the hmm length  for the DOMAIN domain',
       'DOMAIN HMM cov. pos.'  :'Segments of the hmm  covering the protein for the DOMAIN domain',
       'DOMAIN Prot cov.'   :'Ratio of the protein part covered by  the hmm over the hmm length  for the DOMAIN domain'       }

dico_col_names_domain = {
       'DOMAIN Best Match'  :'Name of the sequence presenting the best reciprocal match for the DOMAIN domain',
       'DOMAIN Bit Score'   :'Score of the best reciprocal mach for the DOMAIN domain',
       'DOMAIN Score ratio' :'Ratio between the scores of the best reciprocal match and the followong match for the DOMAIN domain'
       }
custom_col_names_domain = dico_col_names_domain.keys()
custom_col_names_domain_simple = dico_col_names_domain_simple.keys()

#pd.options.mode.copy_on_write = True
parser = argparse.ArgumentParser(description='Reads overview table in the csv format and returns best candidates for each locus')

parser.add_argument('-i', '--input', type=str, required=True, help='Overview table')
parser.add_argument('-o', '--output', type=str, required=True, help='README file')

args = parser.parse_args()
outfile=open(args.output, 'w')

with open("../environment_path.json", "r") as file:
    environment = json.load(file)
pathResources=environment['pathGTDriftResource']

with open("analyse.json", "r") as file:
    analyse = json.load(file)
domains=analyse["domains"]
domains_simple=analyse["domains_simple"]
resources_dir_name=analyse['resources_dir_name']
data_origin = analyse["domain_aln_data_origin"]
outfile.write("Reference alignments:\n")
outfile.write("====================:\n")
for data in data_origin:
       print(data)
       print(pathResources + resources_dir_name + "reference_alignments/" + data + "/" + data_origin[data] + ".fst")
       outfile.write(f"{data:<20}" + " : " + pathResources + resources_dir_name + "reference_alignments/" + data + "/" + data_origin[data] + ".fst\n")

## Reading overview table for prdm9
table = pd.read_csv(args.input, sep=';', dtype=str, header=0)
fields = list(table.columns)

outfile.write("Field definitions:\n")
outfile.write("=================:\n")
for template in  col_names_domain_univ:
       definition =  dico_col_names_univ[template]
       if template in fields:
              outfile.write(f"{template:<20}" + " : " + definition + "\n")
              fields.remove(template)
       else:  
              print("Error :" + template)
              sys.exit("This field is not in the input file")

for domain in domains_simple:
       for template in    custom_col_names_domain_simple:
              definition =  dico_col_names_domain_simple[template]
              template = template.replace("DOMAIN",domain)
              definition = definition.replace("DOMAIN",domain)
              if template in fields:
                     outfile.write(f"{template:<20}" + " : " + definition + "\n")
                     fields.remove(template)
              else:  
                     print("Error :" + template)
                     sys.exit("This field is not in the input file")


for domain in domains:
       for template in    custom_col_names_domain_simple:
              definition =  dico_col_names_domain_simple[template]
              template = template.replace("DOMAIN",domain)
              definition = definition.replace("DOMAIN",domain)
              if template in fields:
                     outfile.write(f"{template:<20}" + " : " + definition + "\n")
                     fields.remove(template)
              else:
                     print("Error :" + template)
                     sys.exit("This field is not in the input file")
       for template in    custom_col_names_domain:
              definition =  dico_col_names_domain[template]
              template = template.replace("DOMAIN",domain)
              definition = definition.replace("DOMAIN",domain)
              if template in fields:
                     outfile.write(f"{template:<20}" + " : " + definition + "\n")
                     fields.remove(template)
              else:
                     print("Error :" + template)
                     sys.exit("This field is not in the input file")


template = 'Unnamed: 0'
if template in fields:
        fields.remove(template)
if len(fields) > 0 :
       print("Missing information for :")
       print(fields)
       sys.exit("Missing information for some fields")

outfile.write("Python packages:\n")
outfile.write("===============:\n")
with open("pyproject.toml", "r") as file:
       s=file.read()
outfile.write(s)
outfile.close()
