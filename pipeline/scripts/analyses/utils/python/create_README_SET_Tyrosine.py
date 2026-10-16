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
Identifies in the protein the three positions in which it should be the three 
catalytic tyrosines of the SET domain of the PRDM9 gene and checks if there are
tyrosines there (1) or not (0).
Reference (query sequence), in this case the SET domain of the PRDM9 human gene:
"CEMCQNFFIDSCAAHGPPTFVKDSAVDKGHPNRSALSLPPGLRIGPSGIPQAGLGVWNEASDLPLGLHFGPYEGRITEDEEAANNGYSWLITKGRNCYEYVDGKDKSWANWMRYVNCARDDEEQNLVAFQYHRQIFYRTCRVIRPGCELLVWYGDEYGQELGIKWGSKWKKELMAGR"

'''


dico_col_names_univ = {
    "Subject_ID":"Protein identifiant",
    "Y276":"1 if a Tyrosine is presenet in position 276",
    "Y341":"1 if a Tyrosine is presenet in position 341",
    "Y357":"1 if a Tyrosine is presenet in position 357"
    }
col_names_domain_univ = dico_col_names_univ.keys()


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



outfile=open(args.output, 'w')



outfile.write(intro)


## Reading overview table for prdm9
table = pd.read_csv(args.input, sep=',', dtype=str, header=0)
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
