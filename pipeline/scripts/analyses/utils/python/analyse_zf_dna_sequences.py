import argparse
import os
import sys
import pandas as pd
import numpy as np
import subprocess
from Bio import SeqIO
from Bio import Cluster
from Bio.Seq import Seq
#from Bio.Cluster import Cluster
from Bio.SeqRecord import SeqRecord
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import AlignIO
from sklearn import cluster
from sklearn.cluster import AgglomerativeClustering

parser = argparse.ArgumentParser()

parser.add_argument('-i', '--input', type=str, required=True, help='zf results')

args = parser.parse_args()

current_working_directory = os.getcwd()

df = pd.read_csv(args.input, sep=';')
print(df.SeqID)
seqids = np.unique(df.SeqID)
print(seqids)
calculator = DistanceCalculator('identity')
for seqid in seqids:
    print(seqid)
    extract_seqid = df[df["SeqID"] == seqid]
    patterns = np.unique(extract_seqid["Pattern num"])
    print(patterns)
    for pattern in patterns : 
        extract_pattern = extract_seqid[extract_seqid["Pattern num"] == pattern]
        #print(str(pattern))
        #print(extract_pattern)
        contigs = np.unique(extract_pattern["Contig"])
        for contig in contigs : 
            extract_contig = extract_pattern[extract_pattern["Contig"] == contig]
            #print(contig)
            #print(extract_contig)
            sequences = extract_contig["dna sequence reading strand"]
            print("SEQUENCE "+seqid+" PATTERN "+str(pattern)+" CONTIG "+ contig+":\n" )
            #print(sequences)
            file_name = seqid+"_p"+str(pattern)+"_"+contig
            description = "Sequence :"+seqid+"; Pattern :"+str(pattern)+"; Contig:"+contig
            my_records = []
            zf_number = 1
            for sequence in sequences:
                dna_seq = Seq(sequence)
                
                record = SeqRecord(dna_seq, id="ZF_"+str(zf_number),description="Zinc finger "+str(zf_number)+";"+description,)
                #print(record)
                my_records.append(record)
                zf_number +=1
            #print(my_records)
            count = SeqIO.write(my_records, file_name+".fasta", "fasta")
            aln = AlignIO.read(open(file_name+".fasta"), 'fasta')
            dm = calculator.get_distance(aln)
            # output file
            fdist = open(file_name+".dist", "w")
            fsilix = open(file_name+".silix", "w")
            fhist = open(file_name+".csv", "w")
            fhm = open(file_name+".hm", "w")
            fhist.write("couple;dist\n")
            fdist.write(str(dm))
            names = dm.names
            min = 1
            max = 0
            somme = 0
            
            for i in range(0,len(names)):
                for j in range(0,i):
                    if dm[i][j] > max:
                        max = dm[i][j]
                    if dm[i][j] < min:
                        min = dm[i][j]  
                    somme +=  dm[i][j] 
            mean = somme /   len(names)               
            for i in range(0,len(names)):
                for j in range(0,i):
                    print("COUPLE "+names[i]+" "+names[j])
                    print(dm[i][j])                       
                    fhist.write(names[i]+"_"+names[j]+";"+str(dm[i][j])+"\n")
                    if dm[i][j] < mean / 10 :
                        fsilix.write(names[i]+" "+names[j]+"\n")
            
            for i in range(0,len(names)):
                fhm.write("\t"+names[i])
            fhm.write("\n")
            for i in range(0,len(names)):
                fhm.write(names[i])    
                for j in range(0,len(names)):
                    fhm.write("\t"+str(dm[i][j]))
                fhm.write("\n")       
            for name in names:
                print("NAME "+name)
                print(dm[0])
            for i in dm:
                print("DEBUG",i)
                print(i)
            print(dm.names)
#file_name+".cluster"
            fhm.close()
            #model = AgglomerativeClustering(n_clusters=None,linkage='complete',distance_threshold=0.05).fit(dm)
            if len(names) > 1 :
                nbclust = 2
                if len(names) > 4 :
                    nbclust = 3
                if len(names) > 5 :
                    nbclust = 4    
                print("./cluster.R "+file_name+".hm 4 ./"+file_name+".clusters")
                subprocess.run (["./cluster.R", current_working_directory+"/"+file_name+".hm", str(nbclust), current_working_directory+"/"+file_name+".clusters"])
            #print(model.labels_)

