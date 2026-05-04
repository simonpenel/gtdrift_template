import argparse

parser = argparse.ArgumentParser(description="Check if the file was ended abruptly")

parser.add_argument('-i','--input', type=str, required=True, help="input file path")

args = parser.parse_args()

check = 0
with open(args.input, 'r') as reader:
    l = reader.readlines()
    if not '\n' in l[-1]:
        check = 1

if check == 1:
    with open(args.input, 'w') as writer:
        writer.write('ERROR DETECTED DURING GENEWISEDB')
