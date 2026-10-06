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
        'SeqID': 'ID of the GeneWise predicted protein', 
        'Assembly':'Genome assembly accession number',
        'Taxid':'Taxonomic identifier',
        'Species':'Species name',
        'Chromosome':'Chromosome or scaffold name',
        'Chr Start':'Start position of the predicted gene in the chromosome', 
        'Chr End':'End position of the predicted gene in the chromosome', 
        'Strand':'Chromosome strand', 
        'Protein Length':'Length of the predicted protein', 
        'ProtRefID':'ID of the representative protein used for the GeneWise prediction', 
        'Nb Introns':'Number of introns in the predicted gene', 
        'Nb Stop/Frameshift':'Number of stops and frameshifts in the predicted gene', 
        'Stop/Shift Positions':'Positions of frameshifts and stop codons reported by GeneWise (in protein coordinates)',
        'ZF Truncated':'Is the zinc finger part truncated?'
       }
col_names_domain_univ = dico_col_names_univ.keys()

dico_col_names_domain_simple = {
       'DOMAIN Query'       :'Reference protein alignment used to search for DOMAIN domains with HMMsearch',
       'DOMAIN E-value'     :'E-value of the best DOMAIN hit (HMMsearch)', 
       'DOMAIN Score'       :'HMMsearch score of the best DOMAIN domain (sum of hit scores in case of merging of split hits)',
       'Nb DOMAIN hits'     :'Number of hits on the DOMAIN domains found by HMMsearch',
       'Nb DOMAIN domains'  :'Number of DOMAIN domains found (after merging of split hits)',
       'DOMAIN domain start':'Start of the best DOMAIN domain in the protein',
       'DOMAIN domain end'  :'End of the best DOMAIN domain in the protein',
       'DOMAIN coverage'    :'Start and end positions of the segments of the reference DOMAIN domain alignment that are covered by at least one hit of the best domain',
       'DOMAIN position'    :'Start and end positions of all hits of the best domain in the reference DOMAIN domain alignment',
       'DOMAIN Prot Length' :'Length of the DOMAIN domain (segments) detected in that protein',  
       'DOMAIN Length'      :'Length of the reference DOMAIN domain',
       'DOMAIN Intron' :'Number of introns in DOMAIN domain',     
       }

dico_col_names_domain_combined = {
       'DOMAIN E-value'     :'Combined E-value of DOMAIN domains (HMMsearch)', 
       'DOMAIN Score'       :'Combined score of DOMAIN domains (HMMsearch)',
       'DOMAIN domain start':'Start of the first DOMAIN domain in the protein',
       'DOMAIN domain end'  :'End of the last DOMAIN domain in the protein',
       'DOMAIN coverage'    :'Fraction of the DOMAIN reference alignment covered by at least one hit',
       'DOMAIN position'    :'Start and end positions of the hits in the reference DOMAIN domain',
       'DOMAIN Prot Length' :'Cumulated length of DOMAIN domains detected in that protein', 
       }

dico_col_names_domain = {
       'DOMAIN paralog Match'   :'ID of the paralog showing the highest similarity to the DOMAIN domain of this protein',
       'DOMAIN paralog Score'   :'Score of the best match among paralogs',
       'DOMAIN paralog Ratio'   :'Ratio (score of the best match)/(score of the 2nd best match)' 
       }
custom_col_names_domain = dico_col_names_domain.keys()
custom_col_names_domains_combined = dico_col_names_domain_combined.keys()
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
data_origin = analyse["domain_references"]
exons  = analyse["exons"]
outfile.write("Reference alignments:\n")
outfile.write("====================\n")
for data in data_origin:
       print(data)
       print(pathResources + resources_dir_name + "reference_alignments/" + data + "/" + data_origin[data] + ".fst")
       outfile.write(f"{data:<20}" + " " + pathResources + resources_dir_name + "reference_alignments/" + data + "/" + data_origin[data] + ".fst\n")
outfile.write("Exons:\n")
outfile.write("=====\n")     
for exon in exons:
       print(pathResources + "ref_align/Prdm9_Metazoa_Reference_alignment/exon_peptides/" + exon + ".fst")
       outfile.write(f"{exon:<20}" + " " + pathResources + "ref_align/Prdm9_Metazoa_Reference_alignment/exon_peptides/" + exon + ".fst\n")
## Reading overview table for prdm9
table = pd.read_csv(args.input, sep=';', dtype=str, header=0)
fields = list(table.columns)

outfile.write("Field definitions:\n")
outfile.write("=================\n")
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
outfile.write("===============\n")
with open("pyproject.toml", "r") as file:
       s=file.read()
outfile.write(s)
outfile.close()
