import argparse
import os
import sys
import pandas as pd
import numpy as np
import subprocess
from Bio import SeqIO
#from Bio import Cluster
from Bio.Seq import Seq
#from Bio.Cluster import Cluster
from Bio.SeqRecord import SeqRecord
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import AlignIO
#from sklearn import cluster
#from sklearn.cluster import AgglomerativeClustering

parser = argparse.ArgumentParser()

parser.add_argument('-i', '--input', type=str, required=True, help='zf results')
parser.add_argument('-o', '--output_dir', type=str, required=True, help='dir of fasta files')

args = parser.parse_args()


df = pd.read_csv(args.input, sep=';')
print(df.SeqID)
seqids = np.unique(df.SeqID)
print(seqids)
#calculator = DistanceCalculator('identity')
# list of output file
fl = open(args.output_dir+"/list_of_files.txt", "w")
for seqid in seqids:
    print(seqid)
    extract_seqid = df[df["SeqID"] == seqid]
    #patterns = np.unique(extract_seqid["Pattern num"])
    #print(patterns)
    contigs = np.unique(extract_seqid["Contig"])
    for contig in contigs :
        extract_contig = extract_seqid[extract_seqid["Contig"] == contig]
        #print(contig)
        #print(extract_contig)
        sequences = extract_contig["dna sequence reading strand"]
        print(sequences)
        zf_names = []
        for name in extract_contig["ZF name"]:
            zf_names.append(name)
            #zf_names = extract_contig["ZF name"]
        print(zf_names)
        print("SEQUENCE "+seqid+" CONTIG "+ contig+":\n" )
        #print(sequences)
        file_name = seqid+"_"+contig
        description = "Sequence :"+seqid+"; Contig:"+contig
        my_records = []
        zf_number = 1
        for sequence in sequences:
            dna_seq = Seq(sequence)

            record = SeqRecord(dna_seq, id=zf_names[zf_number -1],description="Zinc finger ;"+description,)
                #print(record)
            my_records.append(record)
            zf_number +=1
            sequences2 = sequences
            zf_number2 = 1
            for sequence2 in sequences2:
                print("DEBUG ",sequence,sequence2)
                dna_seq2 = Seq(sequence2)
                record2 = SeqRecord(dna_seq2, id=zf_names[zf_number2 -1],description="Zinc finger ;"+description,)
                zf_number2 += 1
                print("DEBUG ",record,record2)
        count = SeqIO.write(my_records, args.output_dir+"/"+file_name+".fasta", "fasta")
        fl.write(file_name+".fasta\n")
fl.close()
