import argparse
import os
import sys
import pandas as pd
import numpy as np
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import AlignIO

parser = argparse.ArgumentParser()

parser.add_argument('-i', '--input', type=str, required=True, help='zf results')
parser.add_argument('-o', '--output_dir', type=str, required=True, help='dir of fasta files')

args = parser.parse_args()

# Calculate the distance in position 3
def distance_pos_3(s1:str,s2:str):
    distance = 0
    l = len(s1)
    if l != len(s2):
      sys.exit("Error in distance_pos_3: sequences have different length.")
    for i in range(0,l-2,3) :
        if s1[i+2] != s2[i+2] :
            distance +=1
    return distance

df = pd.read_csv(args.input, sep=';')
seqids = np.unique(df.SeqID)
# list of output file
fl = open(args.output_dir+"/list_of_files.txt", "w")
for seqid in seqids:
    print(seqid)
    extract_seqid = df[df["SeqID"] == seqid]
    contigs = np.unique(extract_seqid["Contig"])
    for contig in contigs :
        extract_contig = extract_seqid[extract_seqid["Contig"] == contig]
        sequences = extract_contig["dna sequence reading strand"]
        sequences2 = sequences
        zf_names = []
        for name in extract_contig["ZF name"]:
            zf_names.append(name)
        file_name = seqid+"_"+contig
        fsilix=open(args.output_dir+"/"+file_name+".silix", "w")
        description = "Sequence :"+seqid+"; Contig:"+contig
        my_records = []
        zf_number = 0
        for sequence in sequences:
            dna_seq = Seq(sequence)
            record = SeqRecord(dna_seq, id=zf_names[zf_number],description="Zinc finger ;"+description,)
            my_records.append(record)
            zf_number2 = 0
            for sequence2 in sequences2:
                distance = distance_pos_3(sequence.upper(),sequence2.upper())
                if zf_number2 > zf_number:
                    if distance <= 2:
                        print(zf_names[zf_number]+" "+zf_names[zf_number2]+" "+ str(distance))
                        fsilix.write(zf_names[zf_number]+" "+zf_names[zf_number2]+"\n")
                zf_number2 += 1            
            zf_number +=1
        count = SeqIO.write(my_records, args.output_dir+"/"+file_name+".fasta", "fasta")
        fl.write(file_name+".fasta\n")
        fsilix.close()
        fclust=open(args.output_dir+"/"+file_name+".clust", "w")
        subprocess.run(["silixx", str(len(zf_names)),args.output_dir+"/"+file_name+".silix"],stdout=fclust) 
        fclust.close()
        families = {}
        with open(args.output_dir+"/"+file_name+".clust", 'r') as reader:
            for line in reader:
                buf = line.rstrip().split('\t')
                print(buf)
                if len(buf) < 2:
                    continue
                fam = buf[0]
                seq = buf[1]
                if fam in families:
                    families[fam].append(seq)
                else:
                    families[fam] = []
                    families[fam].append(seq)
            print(families)  
            clusters = list(families.keys())
            cluster_max = clusters[0]
            nb_seq_max = len(families[cluster_max])
            for cluster in clusters:
                seqs = families[cluster]
                nb_seq = len(families[cluster])
                print(cluster +": "+str(nb_seq))
                if nb_seq > nb_seq_max :
                    cluster_max = cluster
                    nb_seq_max = nb_seq
            print("Nb of clusters : "+str(len(clusters)))
            print("Cluster max : "+cluster_max)
            print("Size cluster max : "+ str(nb_seq_max)) 
            fclustsummary=open(args.output_dir+"/"+file_name+".clust_summary", "w") 
            fclustsummary.write("Nb of zf : "+str(len(sequences))+"\n")
            fclustsummary.write("Nb of clusters : "+str(len(clusters))+"\n")
            fclustsummary.write("Cluster max : "+cluster_max+"\n")
            fclustsummary.write("Size cluster max : "+ str(nb_seq_max)+"\n")
            seqs = families[cluster]
            fclustsummary.write("Cluster max contents:\n")
            list_array = []
            for seq in seqs:
                fclustsummary.write(seq+"\n")
                arr = seq[0]
                print(arr)
                if not (arr in list_array) :
                    list_array.append(arr)
            fclustsummary.write("Nb of arrays: "+str(len(list_array))+"\n")       
fl.close()
