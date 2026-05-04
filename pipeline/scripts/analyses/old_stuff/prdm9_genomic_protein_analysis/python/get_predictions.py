import argparse
import pandas as pd

parser = argparse.ArgumentParser()

parser.add_argument('-i', '--input', type=str, required=True, help='table summary path')
parser.add_argument('-o', '--output', type=str, required=True, help='output file path')

args = parser.parse_args()

data = pd.read_csv(args.input, sep=';')
data = data.loc[data['Best Match'] == 'PRDM9']
data = data.reset_index(drop=True)

with open(args.output, 'w') as writer:
    for i in range(len(data.index+2)):
        writer.write(data.loc[data.index == i]['SeqID'].iloc[0] + '\n')