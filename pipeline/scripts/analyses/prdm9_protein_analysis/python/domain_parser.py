import pandas as pd
import os
print("Input per domain file "+snakemake.input.per_domain)
print("Output tabulated domain file "+snakemake.output.tabulated_per_domain)
print("Output summary domain file "+snakemake.output.domain_summary)
per_domain_file = snakemake.input.per_domain
tabulated_per_domain_file = snakemake.output.tabulated_per_domain
summary_per_domain_file = snakemake.output.domain_summary

with open(per_domain_file) as reader, open(tabulated_per_domain_file, 'w') as writer:
    for line in reader.readlines():
        if line.startswith('#'):
            del(line)
        else:
            for elt in line.split(maxsplit=23):
                writer.write(f"{elt.strip()}\t")
            writer.write('\n')

with open(tabulated_per_domain_file) as reader, open(summary_per_domain_file, 'w') as writer:
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
            print("\n\n\nDEBUG\n\n COMPARE "+tabulated_per_domain_file+ " WITH ZF_domains_processed\n\n\n\n")
            if tabulated_per_domain_file == 'ZF_domains_processed':
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
        



        
