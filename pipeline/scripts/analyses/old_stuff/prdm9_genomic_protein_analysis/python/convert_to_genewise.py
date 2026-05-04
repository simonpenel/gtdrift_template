import argparse

parser = argparse.ArgumentParser(description="Convert annotated protein outputs to genewise output for concatenation")

parser.add_argument('-i', '--input', type=str, required=True, help='Annotated protein output file path')
parser.add_argument('-s', '--sequences', type=str, required=True, help='Annotated protein sequences file path')
parser.add_argument('-n', '--new', type=str, required=True, help='Annotated sequences output file path')
parser.add_argument('-a', '--accession', type=str, required=True, help='Accession number')
parser.add_argument('-o', '--output', type=str, required=True, help='Processed file path')

args = parser.parse_args()

with open(args.input, 'r') as reader, open(args.sequences, 'r') as sequences, open(args.output, 'w') as writer, open(args.new, 'w') as final:
    rl = reader.readlines()
    sl = sequences.readline()
    while sl:
        if '>' in sl:
            seq = ''
            prot_id = sl.split(' ')[0][1:]
            for line in rl:
                if line.split('\t')[5] == prot_id:
                    line_buffer = line.split('\t')
            sl = sequences.readline()
            while sl and not '>' in sl:
                seq += sl.strip()
                sl = sequences.readline()
            
            loc_id  = line_buffer[0]
            chrom   = line_buffer[1]
            start   = line_buffer[2]
            end     = line_buffer[3]
            strand  = line_buffer[4]
            len_seq = len(seq)
            ref     = line_buffer[7]
            #nb stop
            #stop
            loc_start = line_buffer[8].strip()

            introns = line_buffer[6][1:-1].replace("\'", "").split(', ')
            if not introns:
                introns = []
            else:
                if 'NA' in introns:
                    introns.remove('NA')
                int_adj  = 0
                int_buff = []
                new_int  = []
                for elt in introns:
                    int_start = elt.split('-')[0]
                    int_end   = elt.split('-')[1]
                    new_int.append([int_start, int_end])
                
                new_int.sort(key=lambda x: x[0])
                for elt in new_int:

                    new_start = int(elt[0]) - int(loc_start) - int_adj
                    int_adj  += int(elt[1]) - int(elt[0])
                    int_buff.append(str(round(new_start / 3)))

                introns = int_buff


            int_len = len(introns)
            #pseudo

            writer.write(f'{args.accession}-{chrom}-{loc_id}\t{chrom}\t{start}\t{end}\t{strand}\t{len_seq}\t{prot_id}-{ref}\t0\tNA\t{int_len}\t{";".join(introns)}\n')
            final.write(f'>{args.accession}-{chrom}-{loc_id} {chrom},NA,{";".join(introns)},{len_seq},{start}-{end},{strand},{prot_id}-{ref},No\n')
            final.write(f'{seq}\n')