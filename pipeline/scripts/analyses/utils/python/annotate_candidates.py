import argparse

parser = argparse.ArgumentParser(description='appends unique id to each candidate')

parser.add_argument('-i', '--input', type=str, required=True, help='candidate loci file path')
parser.add_argument('-o', '--output', type=str, required=True, help='Processed file path')

args = parser.parse_args()

with open(args.input, 'r') as reader, open(args.output, 'w') as writer:
    i = 1
    l = reader.readline()
    while l:
        if '>' in l:
            writer.write(f"{l.split(' ')[0]}:Gene_{str(i)} {' '.join(l.split(' ')[1:])}")
            i += 1
        else:
            writer.write(l)
        l = reader.readline()