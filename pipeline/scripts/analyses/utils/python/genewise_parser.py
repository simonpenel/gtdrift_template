import argparse

"""
This script outputs two files: a text file for stats on each
"""

parser = argparse.ArgumentParser()

parser.add_argument('-i', '--input', type=str, required=True, help='genewise output file path')
parser.add_argument('-a', '--accession', type=str, required=True, help='accession number')
parser.add_argument('-o', '--output', type=str, required=True, help='summary output file path')
parser.add_argument('-f', '--fasta', type=str, required=True, help='output fasta file path')

args = parser.parse_args()

## Reading the genewise file
with open(args.input, 'r') as reader:
    line = reader.readline()
    output = {}
    overall_num = 1
    index_genewise = 0
    while line:

        if 'ERROR DETECTED DURING GENEWISEDB' in line:
            with open(args.output, 'w') as writer, open(args.fasta, 'w') as faawriter:
                faawriter.write('>fasta ERROR DETECTED DURING GENEWISEDB\nX\n')
                writer.write('ERROR DETECTED DURING GENEWISEDB')
            quit()

        if '>Results' in line:
            introns      = []
            intronstart  = 0
            intronend    = 0
            pseudoStatus = "No"
            seq          = ''
            stop_shift   = []
            QueryID      = line.strip().split(' ')[2]
            TargetID     = line.strip().split(' ')[4].split(':')[0]
            adj          = line.strip().split(' ')[4].split(':')[1]
            GeneID       = line.strip().split(' ')[4].split(':')[2]
            titlenum     = 1
            index_genewise += 1

            ## Sometimes contains a c: we remove it
            if 'c' in adj:
                adj = adj.replace('c','')

## We make sure the chromosomal positions are in ascending order
            if int(adj.split('-')[0]) > int(adj.split('-')[1]):
                Adjust = [adj.split('-')[1], adj.split('-')]
                strand = '-'
            else:
                Adjust = adj.split('-')
                strand = '+'

## Getting locus title and positions
# The positions are relative to the nucleic sequence, not the chromosome sequence
        elif line.startswith('Gene') and not 'Paras' in line:
            if len(line.split(' ')) == 2:
                title = GeneID+'_'+str(titlenum)+'_'+str(overall_num)
                titlenum    += 1
                print("Parsing " + title)
            elif len(line.split(' ')) > 2:
                TargetRange = [line.strip().split(' ')[1], line.strip().split(' ')[2]]
                if 'pseudogene' in line:
                    pseudoStatus = "Yes"
## Must account for target strand: negative if the range is in decreasing order

                if int(TargetRange[0]) < int(TargetRange[1]):
                    AdjustedRange = TargetRange
                else:
                    AdjustedRange = [TargetRange[1], TargetRange[0]]

## The true positions are the range + the chromosome start
                AdjustedRange[0] = str(int(AdjustedRange[0]) + int(Adjust[0]))
                AdjustedRange[1] = str(int(AdjustedRange[1]) + int(Adjust[0]))
                print(title + " range = " + str(AdjustedRange[0]) + '-' + str(AdjustedRange[1]))

## Reading exons to get introns by deduction
        elif line.strip().startswith('Exon'):
            if intronstart == 0:
                intronstart = line.strip().split(' ')[2]
                alignstart  = line.strip().split(' ')[1]
            elif intronstart != 0:
                intronend = line.strip().split(' ')[1]
                introns.append([intronstart, intronend])
                print("Intron at positions: " + str(intronstart) + '-' + str(intronend))
                intronstart = line.strip().split(' ')[2]

## Getting sequence
## -Modified by Simon :remove_cr
        elif '>' in line and '.pep' in line:
            line = reader.readline()
            print("Reading sequence for " + title)
            while not '/' in line:
                seq += line
                line = reader.readline()
            remove_cr = 0    
            for i in range(len(seq)):
                if seq[i] == '\n':
                    remove_cr += 1
                if seq[i] == 'X':
                    stop_shift.append(i - remove_cr + 1)
                    print("Stop codon or frameshift detected in " + title)
            output[args.accession+'-'+TargetID+'-'+title] = [TargetID, AdjustedRange, strand, QueryID, stop_shift, introns, seq, alignstart, pseudoStatus,index_genewise]
            print(title + " finished parsing")
            overall_num += 1
        line = reader.readline()
        
## Writing on both the fasta output and the text output
with open(args.output, 'w') as writer, open(args.fasta, 'w') as faawriter:

## We go through all elements in our dictionary

    for elt in output:
## Accounting for empty lists by adding NA value
        if output[elt][4] == []:
            output[elt][4].append('NA')
        if output[elt][5] == []:
            output[elt][5].append('NA')

## Correcting intron position based on protein positions
        if output[elt][5] != ['NA']:
            k = 0
            intronadjust = 0
            while k < len(output[elt][5]):
                buffer = output[elt][5][k]
                output[elt][5][k] = str(round((int(output[elt][5][k][0]) - (int(output[elt][7]) + intronadjust))/3))
                intronadjust += (int(buffer[1]) - int(buffer[0]))
                k += 1


## Writing fasta
## Order: prot title, chr, stops-n-shifts, introns, prot length, chromosome range, strand, rep id, pseudostatus
        intronlist = []
        for intron in output[elt][5]:
            intronlist.append(intron)
        shiftlist = []
        for shift in output[elt][4]:
            shiftlist.append(str(shift))
        protlength = len(output[elt][6].replace('\n', ''))
        faawriter.write(f">{elt} {output[elt][0].split(':')[0]},{';'.join(shiftlist)},")
        faawriter.write(f"{';'.join(intronlist)},{str(protlength)},{str(output[elt][1][0])}-{str(output[elt][1][1])},")
        faawriter.write(f"{output[elt][2]},{output[elt][3]},{output[elt][8]},{output[elt][9]}\n")
        faawriter.write(output[elt][6])

## Writing text
## Order: prot title, chr, start, end, strand, prot length, rep id, nb stops-n-shifts, stops-n-shifts, nb introns, introns, pseudostatus
        shift_no_na = [x for x in output[elt][4] if x != 'NA']
        intron_no_na = [x for x in output[elt][5] if x != 'NA']
        writer.write(f"{elt}\t{output[elt][0].split(':')[0]}\t{output[elt][1][0]}\t{output[elt][1][1]}\t{output[elt][2]}\t{str(protlength)}")
        writer.write(f"\t{output[elt][3]}\t{str(len(shift_no_na))}\t{';'.join(shiftlist)}\t{str(len(intron_no_na))}\t{';'.join(intronlist)}\t{output[elt][8]}\t{output[elt][9]}\n")
