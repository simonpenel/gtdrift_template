import argparse
import os

parser = argparse.ArgumentParser(description="Execute genewisedb on sequences: return special format if fails")

parser.add_argument('-i', '--input', type=str, required=True, help='Input sequence file path')
parser.add_argument('-r', '--reference', type=str, required=True, help='Reference sequences file path')
parser.add_argument('-n', '--number', type=str, required=True, help='Incrementation number')
parser.add_argument('-a', '--align', type=str, required=True, help='Max alignment number')
parser.add_argument('-o', '--output', type=str, required=True, help='Output file path')

args = parser.parse_args()

ret = os.system(f'genewisedb {args.reference} {args.input} -prodb -dnas -genes -pseudo -cdna -pep -quiet -init local -subs 1e-6 -indel 1e-6 -pretty -aln {args.align} > {args.output}')
if ret > 0:
    with open(args.output, 'w') as writer:
        writer.write('ERROR DETECTED DURING GENEWISEDB')