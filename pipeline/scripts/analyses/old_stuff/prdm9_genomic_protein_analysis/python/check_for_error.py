import argparse
import os

parser = argparse.ArgumentParser()

parser.add_argument('-i', '--input', type=str, required=True, help='genewise output file path')
parser.add_argument('-o', '--output', type=str, required=True, help='summary output file path')

args = parser.parse_args()

with open(args.input, 'r') as reader:
    l = reader.readlines()
    if '>fasta ERROR DETECTED DURING GENEWISEDB\n' in l:
        with open(args.output+'/errcheck.error_detected', 'w') as output:
            output.write('ERROR DETECTED DURING GENEWISEDB')
    else:
        with open(args.output+'/errcheck.no_error_detected', 'w') as output:
            output.write('No error was detected during genewisedb')
