import pandas as pd
import argparse
from math import isnan



## Reading overview table for prdm9
df = pd.read_csv(snakemake.input[0], sep=';', header=0)

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
df.to_csv(snakemake.output[0], sep=';')
