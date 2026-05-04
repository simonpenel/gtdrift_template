import argparse

"""
This script rewrites the result files of hmm search using tabulation as separator
"""

parser = argparse.ArgumentParser(description='Reads hmm_search tabular output file and writes another file with tab separator instead of whitespaces')

parser.add_argument('-i', '--input', type=str, required=True, help='HMM_search file path')
parser.add_argument('-o', '--output', type=str, required=True, help='New file path')
args = parser.parse_args()

with open(args.input) as reader, open(args.output, 'w') as writer:
    for line in reader.readlines():
        if line.startswith('#'):
            del(line)
        else:
            for elt in line.split(maxsplit=18):
                writer.write(f"{elt.strip()}\t")
            writer.write('\n')
