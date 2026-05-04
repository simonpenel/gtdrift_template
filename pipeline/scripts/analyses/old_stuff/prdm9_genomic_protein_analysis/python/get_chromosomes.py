import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-i', '--input', type=str, required=True, help='loci file path')
parser.add_argument('-o', '--output', type=str, required=True, help='output file path')

args = parser.parse_args()

## Script for preparing blastdbcmd entry_batch file based on candidate loci
with open(args.input, 'r') as reader:
    with open(args.output, 'w') as writer:
        lines = reader.readlines()
        for line in lines[1:]:
            out = line.split('\t')
            if out[4] == '+':
                sign = 'plus'
            else:
                sign = 'minus'
            begin = int(out[2])-50000
            end   = int(out[3])+50000
            if begin < 1:
                begin = 1
            writer.write(f'{out[1]} {str(begin)}-{str(end)} {sign}\n')