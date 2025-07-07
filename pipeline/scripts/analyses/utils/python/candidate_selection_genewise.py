import pandas as pd
import argparse
from math import isnan

#pd.options.mode.copy_on_write = True
parser = argparse.ArgumentParser(description='Reads overview table in the csv format and returns best candidates for each locus')

parser.add_argument('-i', '--input', type=str, required=True, help='Overview table for prdm9')
parser.add_argument('-o', '--output', type=str, required=True, help='Processed file path')

args = parser.parse_args()

## Reading overview table for prdm9
df = pd.read_csv(args.input, sep=';', header=0)

column_list = df.columns.tolist()

print(column_list)
for column in column_list :
    split_col = column.split()
    print(split_col)
    if len(split_col) >1 :
        if split_col[1] == "Query" :
            print(column)
            df = df[df[column] != "0"]
            #df = df.query('English_Score >= 90')
    if len(split_col) >2 :
        if split_col[1] == "Best" and split_col[2] == "Match" :
            print(column)
            query = split_col[0] + " Query"
            df = df[df[column] == df[query]]

df = df.drop(["Unnamed: 0"],axis=1) 
df.to_csv(args.output, sep=';')
