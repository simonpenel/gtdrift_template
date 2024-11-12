import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-i', '--input', type=str, required=True, help='candidates fasta path')
parser.add_argument('-o', '--output', type=str, required=True, help='output directory path')

args = parser.parse_args()

with open(args.input, 'r') as reader:
    l = reader.readline()
    first_line = l
    title = None
    while l:
        if '>' in l:
            if l != first_line:
                chr = title.split(' ')[0].split(':')[0][1:]
                pos = title.split(' ')[0].split(':')[1]
                with open(args.output + '/' + chr + '@' + pos + '.fna', 'w') as writer:
                    writer.write(title + '\n' + seq)
            title = l
            l = reader.readline()
            seq = l
        else:
            seq += l
        l = reader.readline()

    if title is not None :
        chr = title.split(' ')[0].split(':')[0][1:]
        pos = title.split(' ')[0].split(':')[1]
        with open(args.output + '/' + chr + '@' + pos + '.fna', 'w') as writer:
            writer.write(title + '\n' + seq)
