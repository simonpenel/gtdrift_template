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
For each species analyzed, blah blah ....
'''


intro_conclusion = '''
The file 'zf_dna_diversity_summary.csv' contains a table of blah blah..., in CSV format (separated by ;).

'''



dico_col_names_univ = {
    'SeqID': "Protein ID",
    'Status': "0 = the translated DNA sequence has the same length as the annotated protein sequence, 1 = Frameshift: the translated DNA sequence is shorter than the  annotated protein sequence.",
    'Contig': "Chromosome or scaffold name",
    'ZF_Dataset': "The dataset used: it can be : all the ZFs, the ZFs in the longest ZF array ot the ZFs in the largest ZF cluster",
    'ID_set': "The name of the array, or of the cluster according to the dataset used",
    'Nb_ZF': "Number of ZF domains",
    'Nb_ZF_28': "Number of 28 aa ZF domains",
    'Nb_ZF_29': "Number of 29 aa ZF domain",
    'Nb_Arrays': "Nb of arrays, an array is defined as a set of successive ZF domains",
    'Nb_Clusters': "Nb of clusters of ZF domains after a clustering based on similarity at 3rd codon positions",
    'Nb_singletons': "Nb of clusters of ZF domains after a clustering based on similarity at 3rd codon positions",
    'Mean_divS': "Mean synonym divergence (divergence codon - divergence aa)",
    'Mean_ratio_AS': "Mean ratio (aa / mean syno)",
    'Ratio_AS-1': "Mean ratio in position -1",
    'Ratio_AS2': "Mean ratio in position 2",
    'Ratio_AS3': "Mean ratio in position 3",
    'Ratio_AS6': "Mean ratio in position 6",
    'Mean_ratio_AS-1236': "",
    'Mean_ratio_AS_non-1236': "",
    'ZFD': "",
    'Nb ZF hits in annotated protein': ""
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
