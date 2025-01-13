"""
This script rewrites the result files of hmm search using tabulation as separator
"""
with open(snakemake.input[0]) as reader, open(snakemake.output[0], 'w') as writer:
    for line in reader.readlines():
        if line.startswith('#'):
            del(line)
        else:
            for elt in line.split(maxsplit=18):
                writer.write(f"{elt.strip()}\t")
            writer.write('\n')
