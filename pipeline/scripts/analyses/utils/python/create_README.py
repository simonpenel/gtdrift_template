import pandas as pd
import argparse
from math import isnan
import os
import sys
import json
import warnings
import glob
warnings.filterwarnings("error")
intro = '''
For each species analyzed, we screened its annotated proteome with HMMsearch to identify proteins showing similarity with at least one domain of ANALYSE_NAME (E-value<=1e-3).
For each hit, we record its score and E-value, its position in the protein, and its position in the reference domain alignment used by HMMsearch. 
If a domain is split into several hits (e.g. if the protein has a large insertion or deletion within this domain), then the hits corresponding to the different parts of the domain are merged, and their scores and E-values are combined.
If a domain is duplicated in a protein, then only the best hit is selected, except for the following domains for which we combine all copies : COMBINED_DOMAINS.
'''

intro_paralog = '''
To be able to distinguish ANALYSE_NAME DOMAIN domains from DOMAIN domains present in other paralogs of the  PARALOGY_DOMAIN_NAME protein family, we extracted each DOMAIN domain identified, and compared them with HMMscan against a database of DOMAIN domain reference alignments from all PARALOGY_DOMAIN_NAME paralogs.
We recorded the best hit, and the ratio of the  best score over the 2nd best score (to assess to which extent the first hit is better than the following ones).
'''

intro_conclusion = '''
The file 'candidate_homologs.ANALYSE_NAME.proteome.csv' contains a table of the identified homologs, along with information on their domain content, in CSV format (separated by ;).

'''


universal_col_names = ['Unnamed: 0', 'SeqID', 'Assembly', 'Taxid', 'Species']
dico_col_names_univ = {
       'SeqID': 'Protein ID (from the annotated proteome)', 
       'Assembly':'Genome assembly accession number',
       'Taxid':'Taxonomic identifier',
       'Species':'Species name',
       'Chromosome': "Chromosome or scaffold name (according to gff)",
       'Chr Start': "Start position of the gene (according to gff)", 
       'Chr End': "End position of the gene (according to gff)", 
       'Strand': "Chromosome strand (according to gff)", 
       'Protein Length': "Length of the protein", 
       'Pseudo': "[0/1] 1 if the annotated gene is reported as being a pseudogene (according to gff)"
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
       'DOMAIN Length'      :'Length of the reference DOMAIN domain'    
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

with open("../environment_path.json", "r") as file:
    environment = json.load(file)
pathResources=environment['pathGTDriftResource']

with open("analyse.json", "r") as file:
    analyse = json.load(file)

analyse_name = args.input.split('.')[1]
domains_to_merge = str(analyse["domains_to_merge"])
intro = intro.replace("ANALYSE_NAME",analyse_name)
intro = intro.replace("COMBINED_DOMAINS",domains_to_merge)
intro_paralog  = intro_paralog.replace("ANALYSE_NAME",analyse_name)
intro_conclusion = intro_conclusion.replace("ANALYSE_NAME",analyse_name)

outfile=open(args.output, 'w')


domains=analyse["domains"]
domains_simple=analyse["domains_simple"]
resources_dir_name=analyse['resources_dir_name']
data_origin = analyse["domain_references"]
paralogy_domain_names = analyse["paralogy_domain_names"]
outfile.write(intro)

for domain in domains:
       paralogy_domain_name = paralogy_domain_names[domain]
       intro_paralog_tmp  = intro_paralog.replace("PARALOGY_DOMAIN_NAME",paralogy_domain_name)
       intro_paralog_tmp  = intro_paralog_tmp.replace("DOMAIN",domain)
       outfile.write(intro_paralog_tmp)
outfile.write(intro_conclusion)       

outfile.write("Reference alignments  used for protein domain searches with HMMsearch:\n")
outfile.write("=====================================================================:\n")
for data in data_origin:
       print(data)
       print(pathResources + resources_dir_name + "reference_alignments/" + data + "/" + data_origin[data] + ".fst")
       outfile.write(f"{data:<20}" + " : " + pathResources + resources_dir_name + "reference_alignments/" + data + "/" + data_origin[data] + ".fst\n")
outfile.write("\n")
# Recupoer ls findos sur la paralogi
for domain in domains: 
       outfile.write("Alignments used to build the HMM database of " + paralogy_domain_names[domain] + " family database for " + domain +" domain :\n")
       outfile.write("===========================================================================================================================:\n")
       print(pathResources + resources_dir_name + "reference_alignments/" + domain + "/ .fst")
       file_list = glob.glob(pathResources + resources_dir_name + "reference_alignments/" + domain + "/*fst")
       for buf in file_list:
              outfile.write(buf+"\n")      
       print(file_list)

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
              print("debug template " + template)
              definition =  dico_col_names_domain_simple[template]
              print("debug definition " + definition)
              if domain in domains_to_merge:
                     if template in dico_col_names_domain_combined: 
                            definition =  dico_col_names_domain_combined[template] 
              print("debug definition " + definition)      
              template = template.replace("DOMAIN",domain)
              definition = definition.replace("DOMAIN",domain)

              if template in fields:
                     outfile.write(f"{template:<20}" + " : " + definition + "\n")
                     fields.remove(template)
              else:  
                     print("Error :" + template)
                     #sys.exit("This field is not in the input file")


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
