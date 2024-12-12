import pandas as pd
import os
print("Input per domain file "+snakemake.input.per_domain[0])
print("Output tabulated domain file "+snakemake.output.tabulated_per_domain)
print("Output summary domain file "+snakemake.output.domain_summary[0])

with open(snakemake.input.per_domain[0]) as reader, open(snakemake.output.per_domain_tabulated[0], 'w') as writer:
    for line in reader.readlines():
        if line.startswith('#'):
            del(line)
        else:
            for elt in line.split(maxsplit=23):
                writer.write(f"{elt.strip()}\t")
            writer.write('\n')

with open(snakemake.output.per_domain_tabulated[0]) as reader, open(snakemake.output.domain_summary[0], 'w') as writer:
    seq_ID = ''
    newline = ''
    for line in reader.readlines(): 
        current_seq_ID = line.split('\t')[0]
        if current_seq_ID != seq_ID:
            seq_ID = current_seq_ID
            writer.write(newline + '\n')
            newline = ''
            for elt in line.split(maxsplit=23):
                newline += f"{elt.strip()}\t"
        else:
            # overlapping zinc finger domains are merged to create one big domain with multiple repetitions.
            print("\nDEBUG\n\n COMPARE "+snakemake.output.per_domain_tabulated[0]+ " WITH ZF_domains_processed\n\n")
            if snakemake.output.tabulated[0] == 'ZF_domains_processed':
                line_data = line.split(maxsplit=23)
                newline_data = newline.split('\t')
                evalue = line_data[12]
                start = line_data[17]
                end = line_data[18]
                newline_data[12] = str(min(float(evalue), float(newline_data[12])))
                newline_data[17] = str(min(int(start), int(newline_data[17])))
                newline_data[18] = str(max(int(end), int(newline_data[18])))
                newline = '\t'.join(newline_data)
            else:
                writer.write(newline + '\n')
                newline = ''
                for elt in line.split(maxsplit=23):
                    newline += f"{elt.strip()}\t"
                
    writer.write(newline + '\n')
        



        
