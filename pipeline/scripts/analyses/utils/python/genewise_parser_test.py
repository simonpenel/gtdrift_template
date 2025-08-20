import argparse
from Bio.Seq import Seq
import re

"""
This script outputs two files: a text file for stats on each
"""

parser = argparse.ArgumentParser()

parser.add_argument('-i', '--input', type=str, required=True, help='genewise output file path')
parser.add_argument('-a', '--accession', type=str, required=True, help='accession number')
parser.add_argument('-o', '--output', type=str, required=True, help='summary output file path')
parser.add_argument('-f', '--fasta', type=str, required=True, help='output fasta file path')

args = parser.parse_args()

flog  = open("test.log","w")
# output file
f = open("test.output", "w")
f.write("SeqID;Contig;mrna;Status;Nb ZF full;Pattern;Pattern num;Match num;Tandem num;ZF num;ZF name;Start in prot;End in prot;Length;uniformised ZF string;original SF string;Contig;mrna;dna sequence;dna sequence reading strand;dna sequence length\n")


# define patterns
pattern1 = r'..C..C.{12}H...H.{5}'
pattern2 = r'..C..C.{12}H....H.{5}'

# max covering for  covering sf removing
max_covering = 1


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
            dnaseq       = ''
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

        elif '>' in line and '.sp' in line:
                    line = reader.readline()
                    print("Reading dna sequence for " + title)
                    while not '/' in line:
                        dnaseq += line.rstrip()
                        line = reader.readline()
                    output[args.accession+'-'+TargetID+'-'+title] = [TargetID, AdjustedRange, strand, QueryID, stop_shift, introns, seq, alignstart, pseudoStatus,index_genewise,dnaseq]

## Getting sequence
## -Modified by Simon :remove_cr
        elif '>' in line and '.pep' in line:
            line = reader.readline()
            print("Reading sequence for " + title)
            while not '/' in line:
                seq += line.rstrip()
                line = reader.readline()
            remove_cr = 0    
            for i in range(len(seq)):
                if seq[i] == '\n':
                    remove_cr += 1
                if seq[i] == 'X':
                    stop_shift.append(i - remove_cr + 1)
                    print("Stop codon or frameshift detected in " + title)
            #output[args.accession+'-'+TargetID+'-'+title] = [TargetID, AdjustedRange, strand, QueryID, stop_shift, introns, seq, alignstart, pseudoStatus,index_genewise,dnaseq]
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
        protein_seq = output[elt][6]
        print("Protein seq. =           " + output[elt][6])
        print("Dna seq. =               " + output[elt][10])
        dna_seq = output[elt][10]
        bio_dna_seq = Seq(dna_seq)
        translated_dna = bio_dna_seq.translate()
        print("Transl. dna seq. =       " + translated_dna)

        status = "Ok"
        # Create the  part of the translated protein which is correct (i.e util frame shift)
        correct_protein  = []
        flag_fs = False
        for idx in range(0, len(translated_dna)) :
            if translated_dna[idx] ==  protein_seq[idx] :
                correct_protein.append(translated_dna[idx])
            else :
                print("Frameshift!")
                flag_fs = True
                status = "frameshift"
                break
        correct_protein = "".join(correct_protein)
        print("Frameshift =             " + str(flag_fs))       
        print("Corr. Transl. dna seq. = " + correct_protein)
        # Search for zf
        # search patterns

        flog.write("Pattern search in full protein :\n")
        pattern_nb = 1
        list_of_matches = []
        for pattern in [pattern1,pattern2]:
            flog.write("Pattern "+pattern+":\n")
            matches_test = re.finditer(pattern, protein_seq)
            for match in matches_test:
                list_of_matches.append([pattern,match])
        sorted_list_of_matches = sorted(list_of_matches, key=lambda element: element[1].span()[0])   # sort
        flog.write(str(len(sorted_list_of_matches))+ " matches.\n")
        nb_zf_full = len(sorted_list_of_matches)

        flog.write("Pattern search in translated protein :\n")
        pattern_nb = 1
        list_of_matches = []
        for pattern in [pattern1,pattern2]:
            flog.write("Pattern "+pattern+":\n")
            matches_test = re.finditer(pattern, correct_protein)
            for match in matches_test:
                list_of_matches.append([pattern,match])
        sorted_list_of_matches = sorted(list_of_matches, key=lambda element: element[1].span()[0])   # sort
        flog.write(str(len(sorted_list_of_matches))+ " matches.\n")

        # Check for supeprosition
        if len(sorted_list_of_matches) > 0 :
            element = sorted_list_of_matches[0]
            match= element[1]
            cur_start = match.span()[0]
            cur_end = match.span()[1]
            flog.write("Check Match 0 : "+str(match.span()[0])+" - "+str(match.span()[1])+"\n")
            superposed_to_remove = []
            for iel in range(1,len(sorted_list_of_matches)) : 
                element = sorted_list_of_matches[iel]
                match= element[1]
                flog.write("Check Match "+str(iel)+" : "+str(match.span()[0])+" - "+str(match.span()[1])+"\n")
                if match.span()[0] <  cur_end :
                    flog.write(" Warning : superposition\n")
                    if match.span()[0] > cur_start:
                        flog.write("Current match "+str(iel-1)+" starts earlier\n")
                        covering  = cur_end -  match.span()[0]
                        flog.write("Covering = " + str(covering) + "\n")
                        if covering > max_covering :
                            flog.write("Remove new match "+str(iel) +"\n")
                            superposed_to_remove.append(sorted_list_of_matches[iel])
                        else :
                            flog.write("Covering is not  > "+str(max_covering) +"\n")

                    if match.span()[0] == cur_start:
                        flog.write("Check on end\n")
                        if match.span()[1] < cur_end:
                            flog.write("Remove current match "+str(iel-1) +"\n")
                            superposed_to_remove.append(sorted_list_of_matches[iel-1])
                        if match.span()[1] > cur_end:
                            flog.write("Remove new match "+str(iel) +"\n")
                            superposed_to_remove.append(sorted_list_of_matches[iel])
                        if match.span()[1] == cur_end: 
                            sys.exit("Unexpected situation (1)")   

                    if match.span()[0] < cur_start:  
                        sys.exit("Unexpected situation (2)")
                cur_start = match.span()[0]
                cur_end = match.span()[1]

            if len(superposed_to_remove) > 0 :
                flog.write("Superposed zinc finger to be removed:\n")
                for to_remove in superposed_to_remove:
                    flog.write(str(to_remove[1].span()) + "\n")
                    sorted_list_of_matches.remove(to_remove)

                flog.write("New zinc finger list\n")
                for element in sorted_list_of_matches :  
                    flog.write("\tMatch "+str(element[1].span()[0])+ " - "+ str(element[1].span()[1])+"\n")        


        tandem = 0 # will increment for each group of tandem zincfingers
        match_nb = 1 # num of the match in set
        match_tandem_nb = 1 # num of the match in set of tandem zincfingers
        flag_match_ok = True
        flog.write("\nProcessing matches\n")
        for element in sorted_list_of_matches:
            flog.write("Match "+str(element)+"\n")
            if flag_match_ok == False:
                flog.write("Match was flaged as erroneous,  sequence flaged as erroneous\n")
                status = "match transl. problem"
                flag_sequence_ok = False
                continue
            match  = element[1]
            pattern = element[0]
        #for match in matches:
            print("debug match_nb " + str(match_nb))
            if match_nb == 1 :
                current_tandem = match.span()[1]
            else :
                flog.write("Compare current end  "+str(match.span()[0])+" to previous start "+str(current_tandem)+"\n")    
            if match.span()[0] > current_tandem:
                tandem +=1
                #sys.exit("lol")
                match_tandem_nb = 1
            # check for superposiition    
            if match_nb > 2 :
                if match.span()[0] < current_tandem:
                    flog.write("\nWARNING: zinc finger superposition. Covering = " + str( current_tandem - match.span()[0] ) + "\n")
                    if current_tandem - match.span()[0] > max_covering :
                        sys.exit("ERROR: zinc finger superposition")
                    
            flog.write("Match "+str(match_nb)+" "+str(match_tandem_nb)+" "+str(tandem)+" "+str(match.span()[0]) + "-"+str(match.span()[1]))
            current_tandem = match.span()[1]
            # get the matching part of the sequence
            zf_length = match.span()[1] - match.span()[0]
            zf = correct_protein[match.span()[0] : match.span()[1]]
            flog.write("("+str(zf_length)+")\n")
            flog.write("Original sequence : "+str(match.group())+"\n")
            modified = False
            # modify the matching part if needed
            if zf_length == 28 :
                modified = True
                zf = correct_protein[match.span()[0] : match.span()[0]+19] + "X" + correct_protein[match.span()[0] + 19 : match.span()[1]]
                #zf = sequence[match.span()[0] : match.span()[0]+19] + sequence[match.span()[0] + 20 : match.span()[1]]
            flog.write("Modified sequence : "+zf+"\n")
            print("Processing match "+str(match_nb))

            shift_s = 19
            st = match.span()[0]
            en = match.span()[1]
            frame = 0
            raw_seq_extract = ""
            # build the dna sequence of the matching part of the protein
            if modified == False :
                for pos_prot in range(st, en):
                    raw_seq_extract += dna_seq[(pos_prot)*3+0 + frame]
                    raw_seq_extract += dna_seq[(pos_prot)*3+1 + frame]
                    raw_seq_extract += dna_seq[(pos_prot)*3+2 + frame]
            else :
                for pos_prot in range(st, st + shift_s):
                    raw_seq_extract += dna_seq[(pos_prot)*3+0 + frame]
                    raw_seq_extract += dna_seq[(pos_prot)*3+1 + frame]
                    raw_seq_extract += dna_seq[(pos_prot)*3+2 + frame]
                #for pos_prot in range(st + shift_s + 1, en):
                raw_seq_extract += "NNN"

                for pos_prot in range(st + shift_s , en):
                    raw_seq_extract += dna_seq[(pos_prot)*3+0 + frame]
                    raw_seq_extract += dna_seq[(pos_prot)*3+1 + frame]
                    raw_seq_extract += dna_seq[(pos_prot)*3+2 + frame]

            flog.write("Prot:   ")
            for aa in zf:
                flog.write(aa+"  ")
            flog.write("\n")
            flog.write("DNA:    ")
            flog.write(raw_seq_extract+"\n")
            bioseq_dna = Seq(raw_seq_extract)
            bioseq_prot = bioseq_dna.translate()
            flog.write("Trans.: ")
            for aa in bioseq_prot:
                flog.write(aa+"  ")
            flog.write("\n")
            if bioseq_prot != zf :
                sys.exit("translation error")

            flog.write("\n\nCheck OK: Translated sequence and protein sequence are identical.\n\n")
            zfname = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcefghijklmnopqrsruvwxyz"[tandem]+str(match_tandem_nb)
            if modified :
                zfname += "_28"
            else :
                zfname += "_29"
            
            id = elt
            contig = id.split("-")[1]
            mrna = id.split("-")[2]
            f.write(id+";"+contig+";"+mrna+";"+status+";"+str(nb_zf_full)+";"+pattern+";"+str(pattern_nb)+";"+str(match_nb)+";"+str(tandem)+";"+str(match_tandem_nb)+";"+zfname+";"+str(match.span()[0])+";"+str(match.span()[1])+";"+str(zf_length)+";"+zf+";"+str(match.group())+";"+contig+";"+mrna+";"+raw_seq_extract+";"+raw_seq_extract+";"+str(len(raw_seq_extract))+"\n")

            match_nb += 1
            match_tandem_nb += 1
        pattern_nb += 1


## Writing text
## Order: prot title, chr, start, end, strand, prot length, rep id, nb stops-n-shifts, stops-n-shifts, nb introns, introns, pseudostatus
#        shift_no_na = [x for x in output[elt][4] if x != 'NA']
#        intron_no_na = [x for x in output[elt][5] if x != 'NA']
#        writer.write(f"{elt}\t{output[elt][0].split(':')[0]}\t{output[elt][1][0]}\t{output[elt][1][1]}\t{output[elt][2]}\t{str(protlength)}")
#        writer.write(f"\t{output[elt][3]}\t{str(len(shift_no_na))}\t{';'.join(shiftlist)}\t{str(len(intron_no_na))}\t{';'.join(intronlist)}\t{output[elt][8]}\t{output[elt][9]}\n")
