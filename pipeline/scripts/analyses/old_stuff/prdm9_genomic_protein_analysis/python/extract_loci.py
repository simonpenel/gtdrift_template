import argparse
import os
 
parser = argparse.ArgumentParser()

parser.add_argument('-i', '--input', type=str, required=True, help='candidate loci file path')
parser.add_argument('-o', '--output', type=str, required=True, help='output file path')
parser.add_argument('-db', '--database', type=str, required=True, help='database file path')

args = parser.parse_args()

with open(args.input, 'r') as reader, open(args.output, 'w') as writer:
    lines = reader.readlines()[1:]
    with open(args.input+'.batch', 'w') as batch:
        for line in lines:
            spl = line.split('\t')
            batch.write(spl[5]+'\n')
    ret = os.system(f'blastdbcmd -db {args.database} -entry_batch {args.input+".batch"} > {args.output}')
    if ret > 0:
        print('Something strange happened during blastdbcmd, check protein ids')
    os.system(f'rm {args.input+".batch"}')