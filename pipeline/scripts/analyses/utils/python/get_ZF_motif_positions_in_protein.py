# version 1.1
# ============
# les sequences de assemblies.v1.1.json sont correctement traduites, mais il y a des pb avec les matches
# par ex. avec GCF_028858775.2
# version 1.2
# ============
# Ok pour GCF_028858775.2 (
# Strand : - 
# First CDS frame: 2
# Last CDS frame: 0)
#
# Ok pour  assemblies.v1.1.json
#
# Ok pour assemblies.json.GCA_018342105.1
# ---
# Strand = +
# First CDS frame: 0
# Last CDS frame: 1
# => warning sur la longeur / 3
# --
# Strand = -
# First CDS frame: 0
# Last CDS frame: 0
# => warning sur la longeur / 3
# --
# Strand = -
# First CDS frame: 0
# Last CDS frame: 0
# => warning sur la longeur / 3
# --
# Strand = -
# First CDS frame: 1
# Last CDS frame: 0
# => ok
#
# assemblies.json.GCA_013235015.1
# no exons
#
# assemblies.json.GCF_001890085.2
# Problemes:
# Stop car on demande d'arreter 
#                if ratio < 0.2 :
#                 sys.exit("debug")  
# pour debuguer


# version 1.4
# ============
# assemblies.json.GCF_001890085.2  : pb de match

# pb avec la sequence XP_019490596.1 qui contient un frameshift

# ok  avec GCF_033118175.1 qui contient un frameshift

# version 1.5
# ============
# ok avec assemblies.json.GCF_001890085.2  :qui contient un frameshift


# ok  avec GCF_033118175.1 qui contient un frameshift

# pb avec la sequence CAI9150971.1 dans GCA_949782905.1
# Frame first cds = 1
# Frame last cds = 0
# Strand -

# ok with 
# +  0 2 
# -  2 0 
# -  0 0 
# -  1 0

import argparse
import time
import os
import sys
import pandas as pd
from Bio import SeqIO
import re
import numpy as np
#from BCBio import GFF
#from Bio import SeqIO
from Bio.Seq import Seq
parser = argparse.ArgumentParser()

parser.add_argument('-i', '--input', type=str, required=True, help='protein in fasta format')
parser.add_argument('-d', '--dna', type=str,required=True, help='path to fasta dna file')
parser.add_argument('-g', '--gff', type=str,required=True, help='path to gff file')
parser.add_argument('-o', '--output', type=str, required=True, help='output file path')
parser.add_argument('-l', '--log', type=str, required=True, help='log file path')
parser.add_argument('-w', '--warnings', type=str, required=True, help='warnings file path')

args = parser.parse_args()

#-------------------------------------------------------------------------
# Fonction diff_seq
# get the difference between 2 sequences. ratio = 1 if sequences are identical
#-------------------------------------------------------------------------
def diff_seq(s1:str,s2:str):
    score = 0
    l1 = len(s1)
    l2 = len(s2)
    l = min(l1,l2)
    for i in range(0,l) :
        if s1[i] == s2[i] :
            score +=1
    ratio =  (2 * score / (l1 + l2))
    return ratio

#-------------------------------------------------------------------------
# Fonction processExonsGff
# fill 3 dictionnaries :
# dico_exons_pos: (contig, mrna) -> liste des (start,end,strand) des exons
# dico_cds_pos:   (contig, mrna) -> liste des (start,end,strand) des cds
# dico_prot:      protein name   -> liste des (contig,mrna)
#
# Notes on GFF3 Fields:
#
# Fields must be tab-separated. Also, all but the final field in each feature line must contain a value; "empty" columns should be denoted with a '.'
#
#     seqname - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix. Important note: the seqname must be one used within Ensembl, i.e. a standard chromosome name or an Ensembl identifier such as a scaffold ID, without any additional content such as species or assembly. See the example GFF output below.
#     source - name of the program that generated this feature, or the data source (database or project name)
#     feature - feature type name, e.g. Gene, Variation, Similarity
#     start - Start position* of the feature, with sequence numbering starting at 1.
#     end - End position* of the feature, with sequence numbering starting at 1.
#     score - A floating point value.
#     strand - defined as + (forward) or - (reverse).
#     frame - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
#     attribute - A semicolon-separated list of tag-value pairs, providing additional information about each feature.
#
# *- Both, the start and end position are included. For example, setting start-end to 1-2 describes two bases, the first and second base in the sequence.
#
#-------------------------------------------------------------------------
def processExonsGff(gff:str):
    with open(gff, 'r') as reader:
        print("Processing gff... (This may take several minutes)")
        print("Reading gff...")
        for line in reader:
            if line.startswith('#'):
                continue
            split_line = line.split('\t')
            # Processing exons
            if split_line[2] == 'exon':
                contig=split_line[0]
                start=split_line[3]
                end=split_line[4]
                strand=split_line[6]
                exon_info = split_line[8]
                exon_infos = exon_info.split(';')
                parent = "none"
                for info in exon_infos:
                    infos = info.split("=")
                    if infos[0] == "Parent":
                        parent = infos[1]
                if (contig,parent) in dico_exons_pos:
                    dico_exons_pos[(contig,parent)].append([start,end,strand])
                else :
                    dico_exons_pos[(contig,parent)] = []
                    dico_exons_pos[(contig,parent)].append([start,end,strand])

            #processing cds
            if split_line[2] == 'CDS':
                contig=split_line[0]
                start=split_line[3]
                end=split_line[4]
                strand=split_line[6]
                frame=split_line[7]
                cds_gene_id = "none"
                prot_name = "none"
                cds_info   = split_line[8]
                cds_infos = cds_info.split(';')
                parent = "none"
                for info in cds_infos:
                    infos = info.split("=")
                    if infos[0] == "Parent":
                        parent = infos[1]
                flag_gene = 0
                flag_prot = 0
                xrefs = cds_info.split(';')[2].split(",")
                for xref in xrefs:
                    ref=xref.split(':')
                    if ref[0] == "Dbxref=GeneID" or ref[0] == "GeneID" :
                        cds_gene_id = ref[1]
                        flag_gene = 1
                    if ref[0] == "GenBank"  or ref[0]== "Dbxref=GenBank":
                        prot_name = ref[1]
                        flag_prot = 1
                    if ref[0] == "Genbank"  or ref[0] == "Dbxref=Genbank":
                        prot_name = ref[1]
                        flag_prot = 1
                    if ref[0] == "NCBI_GP"  or ref[0] == "Dbxref=NCBI_GP":
                        prot_name = ref[1]
                        flag_prot = 1
                if flag_gene == 0 :
                    xrefs = cds_info.split(';')
                    for xref in xrefs:
                        ref=xref.split('=')
                        if ref[0] == "locus_tag" :
                            cds_gene_id = ref[1]
                            flag_gene = 1
                if cds_gene_id == "none":
                    fwarn.write("WARNING: Unable to get gene id : ")
                    fwarn.write(cds_info)
                if prot_name == "none":
                    fwarn.write("WARNING: Unable to get protein name : ")
                    fwarn.write(cds_info)
                if prot_name != "none" and cds_gene_id != "none":
                    if prot_name in dico_prot:
                        dico_prot[prot_name].append([contig,parent])
                    else :
                        dico_prot[prot_name] = []
                        dico_prot[prot_name].append([contig,parent])
                    if (contig,parent) in dico_cds_pos:
                        dico_cds_pos[(contig,parent)].append([start,end,strand,frame])
                    else :
                        dico_cds_pos[(contig,parent)] = []
                        dico_cds_pos[(contig,parent)].append([start,end,strand,frame])

# output file
f = open(args.output, "w")
# output file
ferr = open(args.output+".ERROR", "w")
# log file
flog = open(args.log, "w")
f.write("SeqID;Contig;mrna;Status;Nb matches in protein;Pattern;Pattern num;Match num;Tandem num;ZF num;ZF name;Start in prot;End in prot;Length;uniformised ZF string;original SF string;Contig;mrna;dna sequence;dna sequence reading strand;dna sequence length\n")
# file of warning
fwarn = open(args.warnings, "w")

# load dna sequence
print("Load dna sequence...")
dico_genome = SeqIO.index(args.dna, 'fasta')
print("Ok.")

# initialise dictionnaries
dico_exons_pos = {}     # Dico ["nom Contig","nom mrna"] => liste [start exon, end exon, strand]
dico_cds_pos = {}		# Dico ["nom Contig","nom mrna"] => liste [start cds, end cds, strand]
dico_prot = {}  		# Dico "nom protein" => liste ["nom Contig","nom mrna"]

processExonsGff(args.gff)

# define patterns
pattern1 = r'..C..C.{12}H...H.{5}'
pattern2 = r'..C..C.{12}H....H.{5}'

# max covering for  covering sf removing
max_covering = 1

# dictionnary dico_sequence:  protein name -> list of ((contig, mrna), sequence_positions)
dico_sequence = {}

# Processing protein sequences
for seq_record in SeqIO.parse(args.input, "fasta"):
    print("SEQUENCE "+seq_record.id)
    # get the protein sequence
    sequence = str(seq_record.seq)
    flog.write("\n===============================\n")
    flog.write(seq_record.id)
    flog.write("\n===============================\n")
    flog.write("\n")
    flog.write("Sequence:\n")
    flog.write(sequence)
    flog.write("\n")
    protein_length = len(sequence)
    flog.write("Protein length from fasta file is "+str(protein_length)+"\n")
    dico_sequence[seq_record.id] = []
    partial_start = False
    partial_end = False
    # flag_identical: True if the translation of dna into protein is correct.
    # if true, the matching sequences should be true too.
    # (used to check the modified pattern results)
    flag_identical = False
    flag_sequence_ok = True
    # info on frame of first and last cds. Useful if sequence is partial
    frame_first_cds = 0
    frame_last_cds = 0
    if seq_record.id in dico_prot:
        status = "Ok"
        contig_mrna_dico  = {}
        # get the uniques couples (contig, mrna) associated to the protein
        for contig_mrna in dico_prot[seq_record.id]:
            if not (contig_mrna[0],contig_mrna[1]) in contig_mrna_dico:
                contig_mrna_dico[(contig_mrna[0],contig_mrna[1])] = "ok"
        unique_contig_mrna = contig_mrna_dico.keys()
        num_cds = 1
        # loop over (contig, mrna) couples associated to the protein
        for contig_mrna in unique_contig_mrna:
            # initialise the array of exons positions in the contig
            sequence_pos = []
            # if (contig, mrna) is in the dictionnary
            if (contig_mrna[0],contig_mrna[1]) in dico_exons_pos:
                flog.write("Transcrit "+str(num_cds) + ": " +contig_mrna[0]+" "+contig_mrna[1]+"\n")
                # get the cds associated to the couple (contig, mrna)
                cds_feat = dico_cds_pos[(contig_mrna[0],contig_mrna[1])]
                # get the strand of the first cds (should be the same for all)
                cds_strand = cds_feat[0][2]
                flog.write("Strand = "+cds_strand+"\n")
                flog.write("CDS   : ")
                # write cds info on log and check if all the cds have the same strand
                for cds in cds_feat:
                    start = int(cds[0])
                    end = int(cds[1])
                    frame = int(cds[3])
                    flog.write("["+str(start)+"-"+str(end)+"] ("+str(frame)+") ")
                    if cds[2] != cds_strand :
                        sys.exit("Error: different strands")
                flog.write("\n")
                # get the range of the protein according to the cds ranges,
                # and the first and last cds
                if  cds_strand  == "+" :
                    first_cds = cds_feat[0]
                    last_cds = cds_feat[len(cds_feat)-1]
                    start_prot = int(first_cds[0])
                    end_prot = int(last_cds[1])
                    frame_first_cds = int(first_cds[3])
                    frame_last_cds = int(last_cds[3])
                elif cds_strand  == "-" :
                    last_cds = cds_feat[0]
                    first_cds = cds_feat[len(cds_feat)-1]
                    start_prot = int(first_cds[0])
                    end_prot = int(last_cds[1])
                    frame_first_cds = int(first_cds[3])
                    frame_last_cds = int(last_cds[3])
                else :
                     sys.exit("Unknown strand")
                flog.write("Protein range (according to CDS) = "+str(start_prot)+"-"+str(end_prot)+"\n")

                # get the exons associated to the couple (contig, mrna)
                exon_features = dico_exons_pos[(contig_mrna[0],contig_mrna[1])]

                # processing partial sequences:
                # ----------------------------
                # check if the first cds covers entirely the first exon
                # brin direct
                #            | ex1   | int1 | ex2 | int
                #        ...123456789012345678901234567890
                #complete       |cds1|      | cds2|
                #               ^
                #               5
                #partial     |    cds1|      | cds2|
                # frame 0    ^
                #            2
                # frame 1     ^
                #             3
                # frame 2      ^
                #              4

                # get first and last exons
                if  cds_strand  == "+" :
                    first_exon = exon_features[0]
                    last_exon = exon_features[len(exon_features)-1]
                else :
                    last_exon = exon_features[0]
                    first_exon = exon_features[len(exon_features)-1]


                flog.write("Protein range  = "+str(start_prot)+"-"+str(end_prot)+"\n")
                flog.write("First CDS frame: "+str(frame_first_cds)+"\n")
                flog.write("Last CDS frame: "+str(frame_last_cds)+"\n")
                start_prot = start_prot + frame_first_cds
                end_prot = end_prot - frame_last_cds

                # if cds_strand  == "+" :
                #     start_prot = start_prot + frame_first_cds
                #     end_prot = end_prot - frame_last_cds
                # else :
                #     start_prot = start_prot - frame_first_cds
                #     end_prot = end_prot + frame_last_cds
                flog.write("Protein range (after check for frame) = "+str(start_prot)+"-"+str(end_prot)+"\n")
                # if cds_strand  == "+" :
                #     end_prot = end_prot - 3 # Suppresion du codon stop
                # flog.write("Protein range (after stop removing) = "+str(start_prot)+"-"+str(end_prot)+"\n")    
                # else :
                #     start_prot = start_prot + 3 # Suppresion du codon stop

                # if  cds_strand  == "+" :
                #     start_prot = start_prot + frame_first_cds
                #     end_prot = 
                # else :
                #     end_prot = end_prot - frame_last_cds
                    
                flog.write("Protein range (after check for frame) = "+str(start_prot)+"-"+str(end_prot)+"\n")
                flog.write("exons : ")

                # Build the array sequence_pos containing the positions of
                # the exons in the dna sequence

                # for exon_feat in exon_features:
                #     flog.write("debug pos\n")  
                #     debug = 0
                #     for pos in range(int(exon_feat[0]), int(exon_feat[1]) +1):
                #         debug += 1
                #     flog.write("debug pos + ["+str(debug)+"]\n")       
                #     flog.write("debug pos\n")  
                #     debug = 0
                #     for pos in range(int(exon_feat[1]) , int(exon_feat[0]) -1 , -1):
                #         debug += 1  
                #     flog.write("debug pos - ["+str(debug)+"]\n")        
                
                if cds_strand == "+" :
                    # direct strand
                    num_exon = 1
                    for exon_feat in exon_features:
                        start = int(exon_feat[0])
                        end = int(exon_feat[1])
                        flog.write("["+str(start)+"-"+str(end)+"]")
                        # get the postion of exon from the cds start to the cds end
                        for pos in range(start, end + 1):
                            if pos >= start_prot and pos < end_prot :
                                sequence_pos.append(pos-1) # pos -1 car l'indexation commence à 0 dans le fichier qui contient l'adn
                        num_exon += 1
                else :
                    # reverse strand
                    reverse = list(reversed(exon_features))
                    num_exon = 1
                    nb_exons = len(exon_features)
                    for exon_feat in exon_features:
                        start = int(exon_feat[1])
                        end = int(exon_feat[0])
                        flog.write("["+str(start)+"-"+str(end)+"]")
                        # get the postion of exon from the cds start to the cds end
                        #flog.write("\ndebug range "+str(start)+ " "+str(end)+"-1 \n")
                        for pos in range(start , end -1 , -1):
                            #flog.write("pos "+str(pos) +" ?  >= "+str(start_prot)+" <= "+str(end_prot)+"\n")
                            # attention a la fin:
                            # on doit etre superieur >= 3start_prot + 3
                            if pos > start_prot + 2 - frame_first_cds  and pos <= end_prot :
                                #flog.write("ADD "+ str(pos-1)+" : "+ debug_raw_seq[pos-1]+"\n")
                                sequence_pos.append(pos - 1) # pos -1 a cause de l'indexation qui commence a 0
                        num_exon +=1
                    reverse_seq=list(reversed(sequence_pos))
                    sequence_pos = reverse_seq
                flog.write("\n")
            else:
                print("ERROR: No exons!")
                #ferr.write("No exons for "+seq_record.id+".\n")
                #continue
                #ferr.close()
                flog.write("Error: No exon in the sequence\n")
                print("No exon in the sequence, jump to  the next sequence")
                flog.write("Sequence flaged as erreoneous.\n")
                flag_sequence_ok = False
                status = "no exon"
                continue

            # Check if dna sequence is ok
            flog.write("DNA sequence length :"+str(len(sequence_pos))+"\n")
            intdiv = int(len(sequence_pos) / 3)
            if intdiv * 3 == len(sequence_pos):
                flog.write("OK : DNA sequence length is multiple of 3\n")
                flog.write("Strand "+cds_strand+"\n")
                flog.write("First CDS frame: "+str(frame_first_cds)+"\n")
                flog.write("Last CDS frame: "+str(frame_last_cds)+"\n")
            else :
                flog.write("WARNING : DNA sequence length is not multiple of 3\n")
                flog.write("Strand "+cds_strand+"\n")
                flog.write("First CDS frame: "+str(frame_first_cds)+"\n")
                flog.write("Last CDS frame: "+str(frame_last_cds)+"\n")
                #sys.exit("DNA sequence length is not multiple of 3 "+cds_strand+" "+str(frame_first_cds)+" "+str(frame_last_cds))

            length_diff_codon = 0
            if len(sequence_pos) == protein_length * 3 :
                flog.write("OK : DNA sequence length is 3 x "+str(protein_length)+"\n")
            else :
                flog.write("WARNING : DNA sequence length is not  3 x "+str(protein_length)+"\n")
                length_diff_codon  = len(sequence_pos)  -  protein_length * 3
            flog.write("Length diff codon : "+str(length_diff_codon)+"\n")    
            # get the dna sequence of the contig
            genome_seq  = dico_genome[contig_mrna[0]]
            raw_seq = genome_seq.seq
            # build the coding sequence of the protein (for information)
            dna_seq = ""
            for i in sequence_pos:
                dna_seq += raw_seq[i]
            flog.write("DNA sequence:\n")
            flog.write(dna_seq)
            flog.write("\n")
            bio_dna_seq = Seq(dna_seq)
            if  cds_strand  == "+" :
                trans_dna_seq = bio_dna_seq.translate()
            else :
                trans_dna_seq = bio_dna_seq.reverse_complement().translate()
            flog.write("Translated DNA sequence:\n")
            flog.write(str(trans_dna_seq))
            flog.write("\n")
            trans_dna_seq_nostop = str(trans_dna_seq)
            #if partial_end == False :
            #    trans_dna_seq_nostop = str(trans_dna_seq)[:-1]
            length_diff = 0 
            if trans_dna_seq_nostop != sequence :
                flog.write("\n\n********\nWarning: translated sequence and protein sequence are different.\n")
                flog.write("Translated sequence:\n")
                flog.write(str(trans_dna_seq_nostop))
                flog.write("\n")
                flog.write(str(len(str(trans_dna_seq_nostop))))
                flog.write("\n")
                flog.write("Original protein sequence:\n")
                flog.write(str(sequence))
                flog.write("\n")
                flog.write(str(len(str(sequence))))
                flog.write("\n")
                length_diff = len(str(sequence)) - len(str(trans_dna_seq_nostop))
                
                ratio = diff_seq(trans_dna_seq_nostop,sequence)
                flog.write("Match ratio  :"+str(ratio))
                flog.write("\n********\n\n")
                if ratio < 0.98 :
                    flog.write("Problem: translated sequence and protein sequence are too different\n")
                    #sequence_pos_translatable  = []
                    correct_protein  = []
                    flag_fs = False
                    for idx in range(0, len(trans_dna_seq_nostop)) :
                        if trans_dna_seq_nostop[idx] ==  sequence[idx] :
                            #sequence_pos_translatable.append(sequence_pos[idx])
                            correct_protein.append(sequence[idx])
                        else :
                            #print("Frameshift!")
                            flag_fs = True
                            status = "frameshift"
                            break
                    correct_protein = "".join(correct_protein)
                    flog.write("Frameshift =             " + str(flag_fs)+"\n")       
                    flog.write("Corr. Transl. dna seq. :\n")   
                    flog.write(correct_protein)   
                    flog.write("\n")   
                else :
                    flog.write("\n\nCheck OK: Translated sequence and protein sequence are 90percent identical.\n\n")
                    correct_protein = sequence
                    full_sequence = sequence  
            else :
                flog.write("\n\nCheck OK: Translated sequence and protein sequence are identical.\n\n")
                flag_identical = True   
                correct_protein = sequence
                full_sequence = sequence

            flog.write("Sequence length diff  = "+str(length_diff)+"\n")
            # magouille
            full_sequence = sequence
            sequence =   correct_protein      

                        
            dico_sequence[seq_record.id].append((contig_mrna,sequence_pos))
            test = dico_sequence[seq_record.id] 
            num_cds +=1

    else:
        print("Sequence inconnue.")
        sys.exit("Sequence inconnue.")

    if flag_sequence_ok == False :
        flog.write("Jump to next sequence\n.")
        f.write(seq_record.id+";"+contig_mrna[0]+";"+contig_mrna[1]+";"+status+";;;;;;;;;;;;;;;;;\n")
        #f.write(seq_record.id+";OK;"+pattern+";"+str(pattern_nb)+";"+str(match_nb)+";"+str(tandem)+";"+str(match_tandem_nb)+";"+zfname+";"+str(match.span()[0])+";"+str(match.span()[1])+";"+str(zf_length)+";"+zf+";"+str(match.group())+";"+seq_genomic[0][0]+";"+seq_genomic[0][1]+";"+raw_seq_extract+";"+str(compl)+";"+str(len(raw_seq_extract))+"\n")
        continue
    
     # search patterns in original protein
    flog.write("Pattern search in original protein:\n")
    pattern_nb = 1
    list_of_matches = []
    for pattern in [pattern1,pattern2]:
        flog.write("Pattern "+pattern+":\n")
        matches_test = re.finditer(pattern, full_sequence)
        for match in matches_test:
            list_of_matches.append([pattern,match])
    sorted_list_of_matches = sorted(list_of_matches, key=lambda element: element[1].span()[0])   # sort
    flog.write(str(len(sorted_list_of_matches))+ " matches.\n")
    nb_zf_full = len(sorted_list_of_matches);

    # search patterns in translated protein
    flog.write("Pattern search in correct protein:\n")
    pattern_nb = 1
    list_of_matches = []
    for pattern in [pattern1,pattern2]:
        flog.write("Pattern "+pattern+":\n")
        matches_test = re.finditer(pattern, sequence)
        for match in matches_test:
            list_of_matches.append([pattern,match])
    sorted_list_of_matches = sorted(list_of_matches, key=lambda element: element[1].span()[0])   # sort
    flog.write(str(len(sorted_list_of_matches))+ " matches.\n")
    # Check if nb of matches is the same than in full protein
    if nb_zf_full != len(sorted_list_of_matches):
        flog.write("Warning: number of matches is different:  frameshift\n")
        status = "frameshift"
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

    # for pattern in [pattern1,pattern2]:
    #
    #     flog.write("\nPattern "+pattern+":\n")
    #     matches = re.finditer(pattern, sequence)

    tandem = 0 # will increment for each group of tandem zincfingers
    match_nb = 1 # num of the match in set
    match_tandem_nb = 1 # num of the match in set of tandem zincfingers
    flag_match_ok = True
    flog.write("\nProcessing " + str(len(sorted_list_of_matches))+" matches\n")
    if len(sorted_list_of_matches) == 0 :
        for seq_genomic in dico_sequence[seq_record.id]:
            f.write(seq_record.id+";"+seq_genomic[0][0]+";"+seq_genomic[0][1]+";"+status+";"+str(nb_zf_full)+";;;;;;;;;;;;;;;;\n")
    for element in sorted_list_of_matches:
        flog.write("Match "+str(element)+"\n")
        if flag_match_ok == False:
            flog.write("Match was flaged as erroneous.\n")
            status = "match transl. problem"
            #flag_sequence_ok = False
            continue
        match  = element[1]
        pattern = element[0]
    #for match in matches:
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
        zf = sequence[match.span()[0] : match.span()[1]]
        flog.write("("+str(zf_length)+")\n")
        flog.write("Original sequence : "+str(match.group())+"\n")
        modified = False
        # modify the matching part if needed
        if zf_length == 28 :
            modified = True
            zf = sequence[match.span()[0] : match.span()[0]+19] + "X" + sequence[match.span()[0] + 19 : match.span()[1]]
            #zf = sequence[match.span()[0] : match.span()[0]+19] + sequence[match.span()[0] + 20 : match.span()[1]]
        flog.write("Modified sequence : "+zf+"\n")
        print("Processing match "+str(match_nb))
        print(seq_record.id+";"+str(match_nb)+";"+str(match.span()[0])+";"+str(match.span()[1])+";"+str(zf_length)+";"+zf+";"+repr(match.group())+"\n")
        for seq_genomic in dico_sequence[seq_record.id]:
            flog.write("For contig = "+seq_genomic[0][0]+"; mrna = "+ seq_genomic[0][1]+"\n")
            # get the dna sequence of the contig
            genome_seq  = dico_genome[seq_genomic[0][0]]
            raw_seq =genome_seq.seq
            # get the positions sequence of the protein
            seq_dna = seq_genomic[1]
            raw_seq_extract = ""
            debug = ""
            for i in seq_dna:
                debug += raw_seq[i]
            bio_debug = Seq(debug)
            if  cds_strand  == "+" :
                trans_debug = bio_debug.translate()
            else :
                trans_debug = bio_debug.reverse_complement().translate()
            shift_s = 19
            st = match.span()[0]
            #st = st + length_diff 
            en = match.span()[1]
            flog.write("Protein start end = "+str(st) + ","+str(en)+"\n")
            frame = 0
            if cds_strand  == "-" :
                flog.write("Protein is on reverse strand.\n")
                flog.write("Modified = "+str(modified)+"\n")
                flog.write("Frame first cds = "+str(frame_first_cds)+"\n")
                flog.write("Frame last cds = "+str(frame_last_cds)+"\n")
                st = len(full_sequence) - match.span()[1] 
                if len(full_sequence) > len(sequence) :
                    flog.write("Protein is not full\n")
                    st = st - 1 
                    st = st + frame_last_cds
                    if modified:
                        st = st - 3
                else :
                    st = st - length_diff 

                en = st + match.span()[1] - match.span()[0]
                #frame = - 3


                # if partial_end == True :
                #     flog.write("Protein is partial at the end.\n")
                #     flog.write("First cds frame "+str(frame_first_cds)+"\n")
                #     flog.write("Last cds frame "+str(frame_last_cds)+"\n")
                #     en = en - 2
                #     st = st - 2

                shift_s = en - 19 - st
                #shift_s -= 1
                flog.write("Protein start end for reverse strand = "+str(st) + ","+str(en)+"\n")
                
            flog.write("Last protein position in dna is " + str(en*3+1 + 2)+"\n")
            flog.write("Length of dna is " + str(len(seq_dna))+"\n") 
            if length_diff_codon == -2 :
                frame = length_diff_codon + 3 
            # build the dna sequence of the matching part of the protein
            if modified == False :
                for pos_prot in range(st, en):
                    raw_seq_extract += raw_seq[seq_dna[(pos_prot)*3+0 + frame]]
                    raw_seq_extract += raw_seq[seq_dna[(pos_prot)*3+1 + frame]]
                    raw_seq_extract += raw_seq[seq_dna[(pos_prot)*3+2 + frame]]
            else :
                for pos_prot in range(st, st + shift_s):
                    raw_seq_extract += raw_seq[seq_dna[(pos_prot)*3+0 + frame]]
                    raw_seq_extract += raw_seq[seq_dna[(pos_prot)*3+1 + frame]]
                    raw_seq_extract += raw_seq[seq_dna[(pos_prot)*3+2 + frame]]
                #for pos_prot in range(st + shift_s + 1, en):
                raw_seq_extract += "NNN"

                for pos_prot in range(st + shift_s , en):
                    raw_seq_extract += raw_seq[seq_dna[(pos_prot)*3+0 + frame]]
                    raw_seq_extract += raw_seq[seq_dna[(pos_prot)*3+1 + frame]]
                    raw_seq_extract += raw_seq[seq_dna[(pos_prot)*3+2 + frame]]

            flog.write("CDS STRAND :"+cds_strand+"\n")
            flog.write("Prot:   ")
            for aa in zf:
                flog.write(aa+"  ")
            flog.write("\n")
            flog.write("DNA:    ")
            flog.write(raw_seq_extract+"\n")
            bioseq_dna = Seq(raw_seq_extract)
            if  cds_strand  == "+" :
                bioseq_prot = bioseq_dna.translate()
                compl = bioseq_dna
            else :
                bioseq_prot = bioseq_dna.reverse_complement().translate()
                compl = bioseq_dna.reverse_complement()
            flog.write("Trans.: ")
            for aa in bioseq_prot:
                flog.write(aa+"  ")
            flog.write("\n")
            ratio = 1
            if str(bioseq_prot) != str(zf) :
                sys.stderr.write("\n\nWarning: translated match and protein match are different.\n")
                sys.stderr.write("Protein:\n")
                sys.stderr.write(str(zf))
                sys.stderr.write("\n")
                sys.stderr.write("Translation:\n")
                sys.stderr.write(str(bioseq_prot))
                sys.stderr.write("\n")
                ratio = diff_seq(str(zf),str(bioseq_prot))
                sys.stderr.write("Match ratio  :"+str(ratio))
                sys.stderr.write("\n")

                flog.write("\n\n********\nWarning: translated match and protein match are different.\n")
                flog.write("Protein:\n")
                flog.write(str(zf))
                flog.write("\n")
                flog.write("Translation:\n")
                flog.write(str(bioseq_prot))
                flog.write("\n")
                flog.write("Match ratio  :"+str(ratio))
                flog.write("\n********\n\n")
                if flag_identical :
                    flog.write("\n********\n\n")
                    flog.write("Error: the translated protein sequence was correct, something is wrong with the patterns.\n")
                    flog.write("\n********\n\n")
            if ratio < 0.9 :
                    flog.write("Error: translated match  and protein match  are too different\n")
                    #flog.write("Match is flaged as erroneous\n")
                    print("Error: translated match and protein match are too different")
                    #sys.exit("Error: translated match and protein  are too different")
                    f.write(seq_record.id+";"+seq_genomic[0][0]+";"+seq_genomic[0][1]+";"+status+";"+str(nb_zf_full)+";"+pattern+";"+str(pattern_nb)+";"+str(match_nb)+";"+str(tandem)+";"+str(match_tandem_nb)+";"+"UNCORRECTLY TRANSLATED MATCH"+";"+str(match.span()[0])+";"+str(match.span()[1])+";"+str(zf_length)+";"+zf+";"+str(match.group())+";"+seq_genomic[0][0]+";"+seq_genomic[0][1]+";"+raw_seq_extract+";"+str(compl)+";"+str(len(raw_seq_extract))+"\n")
                    #print("Match is flaged as erroneous")
                    #flag_match_ok = False
                    #continue
                    
            else:
                flog.write("\n\nCheck OK: Translated match and protein match are identical.\n\n")
                zfname = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcefghijklmnopqrsruvwxyz"[tandem]+str(match_tandem_nb)
                if modified :
                    zfname += "_28"
                else :
                    zfname += "_29"
                f.write(seq_record.id+";"+seq_genomic[0][0]+";"+seq_genomic[0][1]+";"+status+";"+str(nb_zf_full)+";"+pattern+";"+str(pattern_nb)+";"+str(match_nb)+";"+str(tandem)+";"+str(match_tandem_nb)+";"+zfname+";"+str(match.span()[0])+";"+str(match.span()[1])+";"+str(zf_length)+";"+zf+";"+str(match.group())+";"+seq_genomic[0][0]+";"+seq_genomic[0][1]+";"+raw_seq_extract+";"+str(compl)+";"+str(len(raw_seq_extract))+"\n")
            #f.write(seq_record.id+";"+seq_genomic[0][0]+";"+seq_genomic[0][1]+";"+status+";"+str(nb_zf_full)+";"+pattern+";"+str(pattern_nb)+";"+str(match_nb)+";"+str(tandem)+";"+str(match_tandem_nb)+";"+zfname+";"+str(match.span()[0])+";"+str(match.span()[1])+";"+str(zf_length)+";"+zf+";"+str(match.group())+";"+seq_genomic[0][0]+";"+seq_genomic[0][1]+";"+raw_seq_extract+";"+str(compl)+";"+str(len(raw_seq_extract))+"\n")
            #f.write(seq_record.id+";"+contig_mrna[0]+";"+contig_mrna[1]+";"+status+";"+pattern+";"+str(pattern_nb)+";"+str(match_nb)+";"+str(tandem)+";"+str(match_tandem_nb)+";"+zfname+";"+str(match.span()[0])+";"+str(match.span()[1])+";"+str(zf_length)+";"+zf+";"+str(match.group())+";"+seq_genomic[0][0]+";"+seq_genomic[0][1]+";"+raw_seq_extract+";"+str(compl)+";"+str(len(raw_seq_extract))+"\n")
        match_nb += 1
        match_tandem_nb += 1
    pattern_nb += 1
