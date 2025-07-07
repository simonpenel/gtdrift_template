import pandas as pd
import argparse
from math import isnan

#pd.options.mode.copy_on_write = True
parser = argparse.ArgumentParser(description='Reads overview table in the csv format and returns best candidates for each locus')

parser.add_argument('-p', '--protein', type=str, required=True, help=' ')
parser.add_argument('-g', '--genomic', type=str, required=True, help=' ')
parser.add_argument('-o', '--output', type=str, required=True, help='Processed file path')

args = parser.parse_args()

## Reading protein
table_protein = pd.read_csv(args.protein, sep=';', header=0)

## Sorting by chromosome then by start position
sorted_protein = table_protein.sort_values(by = ['Chromosome','Chr Start'])

## Getting rid of extra index column and resetting index
sorted_protein.drop(columns='Unnamed: 0', inplace=True)
sorted_protein.reset_index(drop=True, inplace=True)

## Reading genomic
table_genomic = pd.read_csv(args.genomic, sep=';', header=0)

## Sorting by chromosome then by start position
sorted_genomic = table_genomic.sort_values(by = ['Chromosome','Chr Start'])

## Getting rid of extra index column and resetting index
sorted_genomic.drop(columns='Unnamed: 0', inplace=True)
sorted_genomic.reset_index(drop=True, inplace=True)

## Overlap function for comparison between two loci
def overlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))


#Comparison between loci i and j
i = 0
id = 0
find_in_prot = []
find_in_genewise = []
df = pd.DataFrame() 
while i < len(sorted_protein.index):
    chromosome = sorted_protein.loc[i].Chromosome
    start = sorted_protein.loc[i]['Chr Start']
    end = sorted_protein.loc[i]['Chr End']
    j = 0
    while j < len(sorted_genomic.index):
        chromosome_g = sorted_genomic.loc[j].Chromosome
        start_g = sorted_genomic.loc[j]['Chr Start']
        end_g = sorted_genomic.loc[j]['Chr End']
        if chromosome_g == chromosome :
            over = overlap([start_g,end_g],[start,end])
            if over > 0 :
                print("Fusion for "+ sorted_protein.loc[i].SeqID)
                print("       and "+ sorted_genomic.loc[j].SeqID)
                print("Protein : start = " + str(start) + ", end = " + str(end) + ",Chromo = "+ chromosome)
                print("Genomic : start = " + str(start_g) + ", end = " + str(end_g) + ",Chromo = "+ chromosome_g)
                if not sorted_genomic.loc[j].SeqID in find_in_prot:
                    find_in_prot.append(sorted_genomic.loc[j].SeqID)
                if not sorted_protein.loc[i].SeqID in find_in_genewise:
                    find_in_genewise.append(sorted_protein.loc[i].SeqID)    
                id = id + 1
                fusion_row = {"UniqueID": id}
                fusion_row.update(sorted_protein.loc[i]) 
                df = df._append(fusion_row, ignore_index=True)
                fusion_row = {"UniqueID": id}
                fusion_row.update(sorted_genomic.loc[j]) 
                df = df._append(fusion_row, ignore_index=True)
        j = j + 1
    i = i + 1
df.to_csv(args.output, sep=';')
fbase=args.output.removesuffix(".csv")
for id in find_in_prot:
    sorted_genomic = sorted_genomic.drop(sorted_genomic[sorted_genomic.SeqID == id].index)
for id in find_in_genewise:
    sorted_protein = sorted_protein.drop(sorted_protein[sorted_protein.SeqID == id].index)    
sorted_protein.to_csv(fbase+".protein_only.csv", sep=';')
sorted_genomic.to_csv(fbase+".genewise_only.csv", sep=';')
 
            
  
                



