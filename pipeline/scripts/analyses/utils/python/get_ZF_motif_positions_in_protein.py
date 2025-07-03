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
                        dico_cds_pos[(contig,parent)].append([start,end,strand])
                    else :
                        dico_cds_pos[(contig,parent)] = []
                        dico_cds_pos[(contig,parent)].append([start,end,strand])  


# output file
f = open(args.output, "w")
# log file 
flog = open(args.log, "w") 
f.write("SeqID;Pattern;Pattern num;ZF num;Start in prot;End in prot;Length;uniformised ZF string;original SF string;Contig;mrna;dna sequence;dna sequence reading strand;dna sequence length\n")
# file of warning                                          
fwarn = open(args.warnings, "w") 

# load dna sequence
print("load dna sequence...")
dico_genome = SeqIO.index(args.dna, 'fasta')
print("Ok.")

dico_exons_pos = {}     # Dico ["nom Contig","nom mrna"] => liste [start exon, end exon, strand]
dico_cds_pos = {}		# Dico ["nom Contig","nom mrna"] => liste [start cds, end cds, strand]
dico_prot = {}  		# Dico "nom protein" => liste ["nom Contig","nom mrna"]
                                                          
processExonsGff(args.gff)

pattern1 = r'..C..C.{12}H...H.{5}'
pattern2 = r'..C..C.{12}H....H.{5}' 



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
    flog.write(sequence)
    flog.write("\n")    
    dico_sequence[seq_record.id] = []
    if seq_record.id in dico_prot:
        contig_mrna_dico  = {}
        # Get the uniques couples (contig, mrna) associated to the protein
        for contig_mrna in dico_prot[seq_record.id]:
            if not (contig_mrna[0],contig_mrna[1]) in contig_mrna_dico:
                contig_mrna_dico[(contig_mrna[0],contig_mrna[1])] = "ok"
        unique_contig_mrna = contig_mrna_dico.keys()
        num_cds = 1
        for contig_mrna in unique_contig_mrna:
            sequence_pos = []
            print("SEQUENCE "+seq_record.id+" / transcrit "+str(num_cds))
            print(contig_mrna)
            if (contig_mrna[0],contig_mrna[1]) in dico_exons_pos:
                flog.write("Transcrit "+str(num_cds) + ": " +contig_mrna[0]+" "+contig_mrna[1]+"\n")
                # get the cds associated to the couple (contig, mrna)
                cds_feat = dico_cds_pos[(contig_mrna[0],contig_mrna[1])]
                flog.write("CDS   : ")
                cds_strand = cds_feat[0][2]
                flog.write("Strand = "+cds_strand+"\n") 
                for cds in cds_feat:
                    start = int(cds[0])
                    end = int(cds[1])
                    flog.write("["+str(start)+"-"+str(end)+"]")
                    if cds[2] != cds_strand :
                        sys.exit("Error: different strands")
                flog.write("\n")
                if  cds_strand  == "+" :
                    first_cds = cds_feat[0]
                    last_cds = cds_feat[len(cds_feat)-1]
                    start_prot = int(first_cds[0])
                    end_prot = int(last_cds[1])
                elif cds_strand  == "-" :
                    last_cds = cds_feat[0]
                    first_cds = cds_feat[len(cds_feat)-1]
                    start_prot = int(first_cds[0])
                    end_prot = int(last_cds[1])
                else :
                     sys.exit("Unknown strand")
                flog.write("Protein range = "+str(start_prot)+"-"+str(end_prot)+"\n") 
                flog.write("exons : ") 
                # get the exons associated to the couple (contig, mrna)
                exon_features = dico_exons_pos[(contig_mrna[0],contig_mrna[1])]
                if cds_strand == "+" :
                    for exon_feat in exon_features:
                        start = int(exon_feat[0])
                        end = int(exon_feat[1])
                        flog.write("["+str(start)+"-"+str(end)+"]")
                        # get the postion of exon from the cds start to the cds end
                        for pos in range(start, end + 1):
                            if pos >= start_prot - 1 and pos < end_prot :
                                sequence_pos.append(pos)
                else :
                    reverse = list(reversed(exon_features))
                    for exon_feat in exon_features:
                        start = int(exon_feat[1])
                        end = int(exon_feat[0])
                        flog.write("["+str(start)+"-"+str(end)+"]")
                        # get the postion of exon from the cds start to the cds end
                        for pos in range(start -1 , end -2 , -1):
                            if pos >= start_prot - 1  and pos < end_prot :
                                sequence_pos.append(pos)  
                    reverse_seq=list(reversed(sequence_pos))
                    sequence_pos = reverse_seq                                  
                flog.write("\n")
            else:
                print("Pas d'exons")
                sys.exit("Pas d'exons.")

            flog.write("DNA SEQUENCE LENGTH :"+str(len(sequence_pos))+"\n")   
            # get the dna sequence of the contig 
            genome_seq  = dico_genome[contig_mrna[0]]
            raw_seq = genome_seq.seq
            # build the coding sequence of the protein (for information)
            dna_seq = ""
            for i in sequence_pos:
                dna_seq += raw_seq[i]
            flog.write("DNA SEQUENCE:\n")
            flog.write(dna_seq)
            flog.write("\n")   
            bio_dna_seq = Seq(dna_seq)
            if  cds_strand  == "+" :
                trans_dna_seq = bio_dna_seq.translate()
            else :
                trans_dna_seq = bio_dna_seq.reverse_complement().translate()
            flog.write("TRANSLATED DNA SEQUENCE:\n")
            flog.write(str(trans_dna_seq))
            flog.write("\n") 
            trans_dna_seq_nostop = str(trans_dna_seq)[:-1]
            if trans_dna_seq_nostop != sequence :
            
                flog.write("\n\n********\nWarning: translated sequence and protein sequence are different.\n")
                flog.write(str(trans_dna_seq_nostop))
                flog.write("\n")
                flog.write(str(len(str(trans_dna_seq_nostop))))
                flog.write("\n")  
                flog.write(str(sequence))
                flog.write("\n")   
                flog.write(str(len(str(sequence))))    
                flog.write("\n")        
                ratio = diff_seq(trans_dna_seq_nostop,sequence)  
                #ratio = SequenceMatcher(None,str(trans_dna_seq_nostop), str(sequence)).ratio()
                flog.write("Match ratio  :"+str(ratio))
                flog.write("\n********\n\n")
                if ratio < 0.9 :                      
                        flog.write("Error: translated sequence and protein sequence are too different")
                        sys.exit("Error: translated sequence and protein sequence are too different")                        
            else :
                flog.write("Sequences are identical")
 			               
            
            dico_sequence[seq_record.id].append((contig_mrna,sequence_pos))
            num_cds +=1        
    else:
        print("Sequence inconnue.") 
        sys.exit("Sequence inconnue.")

    # get the protein sequence
    #sequence = str(seq_record.seq)

    # search patterns
    pattern_nb = 1
    for pattern in [pattern1,pattern2]:
        flog.write("\nPattern "+pattern+":\n")  
        matches = re.finditer(pattern, sequence)
        match_nb = 1
        for match in matches:
            flog.write("\nMatch "+str(match_nb)+" "+str(match.span()[0]) + "-"+str(match.span()[1]))  
            # get the matching part of the sequence
            zf_length = match.span()[1] - match.span()[0]
            zf = sequence[match.span()[0] : match.span()[1]]
            flog.write("("+str(zf_length)+")\n")
            flog.write("Original sequence : "+str(match.group())+"\n")
            modified = False
            # modify the matching part if needed
            if zf_length > 28 :
                modified = True
                zf = sequence[match.span()[0] : match.span()[0]+19] + sequence[match.span()[0] + 20 : match.span()[1]]
            flog.write("Modified sequence : "+zf+"\n")   
            print(seq_record.id+";"+str(match_nb)+";"+str(match.span()[0])+";"+str(match.span()[1])+";"+str(zf_length)+";"+zf+";"+repr(match.group()))
            for seq_genomic in dico_sequence[seq_record.id]:
                flog.write(seq_genomic[0][0]+" "+ seq_genomic[0][1]+":\n")   
                # get the dna sequence of the contig
                genome_seq  = dico_genome[seq_genomic[0][0]]
                print(genome_seq)
                raw_seq =genome_seq.seq
                # get the positions sequence of the protein 
                seq_dna = seq_genomic[1]
                raw_seq_extract = ""
                
                
                debug = ""
                for i in seq_dna:
                    debug += raw_seq[i]
                flog.write("DEBUG DNA SEQUENCE:\n")
                flog.write(debug)
                flog.write("\n")   
                bio_debug = Seq(debug)
                flog.write("CDS STRAND :"+cds_strand+"\n")
                if  cds_strand  == "+" :
                    trans_debug = bio_debug.translate()
                else :
                    trans_debug = bio_debug.reverse_complement().translate()
                flog.write("DEBUG TRANSLATED DNA SEQUENCE:\n")
                flog.write(str(trans_debug))
                flog.write("\n")   
                shift_s = 19
                st = match.span()[0]
                en = match.span()[1]
                flog.write("debug start end = "+str(st) + ","+str(en)+"\n")  
                if cds_strand  == "-" :
                    st = len(sequence) - match.span()[1] + 1 
                    en = st + match.span()[1] - match.span()[0]
                    # st + x = en - 19 => x = end - 19 -st  
                    shift_s = en - 19 - st
                flog.write("debug start end = "+str(st) + ","+str(en)+"\n") 
                
                # build the dna sequence of the matching part of the protein
                if modified == False :
                    for pos_prot in range(st, en):
                        raw_seq_extract += raw_seq[seq_dna[(pos_prot)*3+0]]
                        raw_seq_extract += raw_seq[seq_dna[(pos_prot)*3+1]]
                        raw_seq_extract += raw_seq[seq_dna[(pos_prot)*3+2]]
                else :
                    for pos_prot in range(st, st + shift_s):
                        raw_seq_extract += raw_seq[seq_dna[(pos_prot)*3+0]]
                        raw_seq_extract += raw_seq[seq_dna[(pos_prot)*3+1]]
                        raw_seq_extract += raw_seq[seq_dna[(pos_prot)*3+2]]
                    for pos_prot in range(st + shift_s + 1, en):
                        raw_seq_extract += raw_seq[seq_dna[(pos_prot)*3+0]]
                        raw_seq_extract += raw_seq[seq_dna[(pos_prot)*3+1]]
                        raw_seq_extract += raw_seq[seq_dna[(pos_prot)*3+2]]
                                                
                for aa in zf:
                    flog.write(aa+"  ")    
                flog.write("\n")    
                flog.write(raw_seq_extract+"\n")   
                bioseq_dna = Seq(raw_seq_extract)

                flog.write("CDS STRAND :"+cds_strand+"\n")
                if  cds_strand  == "+" :
                    bioseq_prot = bioseq_dna.translate()
                    compl = bioseq_dna
                else :
                    bioseq_prot = bioseq_dna.reverse_complement().translate()  
                    compl = bioseq_dna.reverse_complement()              
                #bioseq_prot = compl.translate()
                flog.write("Translation:\n")
                for aa in bioseq_prot:
                    flog.write(aa+"  ")    
                flog.write("\n") 
                if str(bioseq_prot) != str(zf) :
                    sys.stderr.write("\n\nWarning: translated sequence and protein sequence are different.\n")
                    sys.stderr.write("Protein:\n")
                    sys.stderr.write(str(zf))
                    sys.stderr.write("\n")
                    sys.stderr.write("Translation:\n")
                    sys.stderr.write(str(bioseq_prot))
                    sys.stderr.write("\n")
                    ratio = diff_seq(str(zf),str(bioseq_prot)) 
                    #ratio = SequenceMatcher(None, str(zf), str(bioseq_prot)).ratio()
                    sys.stderr.write("Match ratio  :"+str(ratio))
                    sys.stderr.write("\n")
                    
                    flog.write("\n\n********\nWarning: translated sequence and protein sequence are different.\n")
                    flog.write("Protein:\n")
                    flog.write(str(zf))
                    flog.write("\n")
                    flog.write("Translation:\n")
                    flog.write(str(bioseq_prot))
                    flog.write("\n")
                    flog.write("Match ratio  :"+str(ratio))
                    flog.write("\n********\n\n")
                    if ratio < 0.9 :
                        flog.write("Error: translated sequence and protein sequence are too different")                                         
                        sys.exit("Error: translated sequence and protein sequence are too different")
                else:
                    print("\n\nCheck OK: Translated sequence and protein sequence are identical.\n\n")
                f.write(seq_record.id+";"+pattern+";"+str(pattern_nb)+";"+str(match_nb)+";"+str(match.span()[0])+";"+str(match.span()[1])+";"+str(zf_length)+";"+zf+";"+str(match.group())+";"+seq_genomic[0][0]+";"+seq_genomic[0][1]+";"+raw_seq_extract+";"+str(compl)+";"+str(len(raw_seq_extract))+"\n")
            
            match_nb += 1
        pattern_nb += 1
            
