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

# summay of zf in some positions

def zfsum(zf:[], positions:[]):
    sum = 0
    for pos in positions:
        sum = sum + zf[pos]
    total = 0
    for val in zf:
        total = total + val
    sum = sum / total
    return sum
 

# Calculate the zfd
def zfd(seqs:[]):
    zfd = []
    l = len(seqs[0])
    nbseq =  len(seqs)
    for i in range(0,l) :
        site = []
        print(i)
        for seq in seqs:
            site.append(seq[i])
        print("Site = ",site)   
        df = pd.DataFrame({'aa':site})
        unic_site = list(np.unique(df.aa))
        print("Unic = ",unic_site) 
        div = 0
        for aa in unic_site:
            nb_aa = site.count(aa)
            print(aa+" => "+str(nb_aa))
            div  = div + (nb_aa / nbseq) ** 2
        div = 1 - div
        print("Site "+str(i) +" => "+str(div))
        zfd.append(div)
    print(zfd)
    return zfd


# Calculate the zfd on codon
def zfd_codon(seqs:[]):
    zfd = []
    l = len(seqs[0])
    nbseq =  len(seqs)
    for i in range(0,l,3) :
        site = []
        print(i)
        for seq in seqs:
            codon = seq[i]+seq[i+1]+seq[i+2]
            codon = codon.upper()
            site.append(codon)
        print("Site = ",site)  
        df = pd.DataFrame({'codon':site})
        unic_site = list(np.unique(df.codon))
        print("Unic = ",unic_site) 
        div = 0
        for codon in unic_site:
            nb_codon = site.count(codon)
            print(codon+" => "+str(nb_codon))
            div  = div + (nb_codon / nbseq) ** 2
        div = 1 - div
        print("SIte "+str(i) +" => "+str(div))
        zfd.append(div)
    print(zfd)
    return zfd
    
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
    
def write_divindex(zfd_data, zfd_codon_data, zfd_name):    
    fclustsummary.write('{:45}'.format("# Position : "))
    ipos = 1 
    for div in zfd_data:
        fclustsummary.write('%5d' % (ipos))
        ipos += 1       
    fclustsummary.write("\n")      
          
    fclustsummary.write('{:45}'.format("# " + zfd_name + " div. index aa : "))
    for div in zfd_data:
        fclustsummary.write('%5.2f' % (div))
    fclustsummary.write("\n")
    
    fclustsummary.write('{:45}'.format("# " + zfd_name + " div. index  codons : "))
    for div in zfd_codon_data:
        fclustsummary.write('%5.2f' % (div))
    fclustsummary.write("\n")    

def write_divindex_pos(zfd_data,zfd_name):
            zfd_pos = zfsum(zfd_data,positions_contact )
            pos_string = ' '.join(str(item+1) for item in positions_contact)
            #fclustsummary.write(zfd_name + " d. i. for " + pos_string)
            fclustsummary.write('{:45}'.format("# " + zfd_name + " sum/tot " + pos_string+ " : "))
            fclustsummary.write('%6.3f' % (zfd_pos))
            fclustsummary.write("\n")

def calcul_synonym_divindex(zfd_data, zfd_codon_data, zfd_name):    
    zdf_syno = []
    zdf_ratio = []
    total_syno = 0.0
    total_ratio = 0.0
    for i in range(0,len(zfd_data)):
        syno = zfd_codon_data[i] - zfd_data[i]
        if syno < 0:
            if (-syno) < 0.000000001:
                syno = -syno
            else:
                sys.exit("Error in calcul_synonym_divindex, nefatve value")
        if syno < 0.000000001:
            syno = 0
        
        ratio = 0
        if syno > 0:
            ratio = zfd_data[i] / syno
            total_ratio += ratio
        zdf_ratio.append(ratio)   
        total_syno += syno        
        zdf_syno.append(syno)
    mean_syno = total_syno / len(zfd_data)   
    mean_ratio = total_ratio / len(zfd_data)   
    fclustsummary.write('{:45}'.format("# " + zfd_name + " syno (codon - aa) : "))
    for div in zdf_syno:
        fclustsummary.write('%5.2f' % (div))
    fclustsummary.write("\n") 
    fclustsummary.write('{:45}'.format("# "+ zfd_name + " ratio (aa / syno) : "))
    for div in zdf_ratio:
        fclustsummary.write('%5.2f' % (div))
    fclustsummary.write("\n")       
    fclustsummary.write('{:45}'.format("# Mean syno : "))     
    fclustsummary.write('%5.2f' % (mean_syno))       
    fclustsummary.write("\n") 
    fclustsummary.write('{:45}'.format("# Mean ratio : "))     
    fclustsummary.write('%5.2f' % (mean_ratio))       
    fclustsummary.write("\n")
    
    for pos in positions_contact:
        fclustsummary.write('{:45}'.format("# " + zfd_name + " position "+str(pos+1)))
        fclustsummary.write(" Ratio : ")     
        fclustsummary.write('%5.2f' % (zdf_ratio[pos])) 
        fclustsummary.write(" Mean ratio : ")    
        fclustsummary.write('%5.2f' % (mean_ratio)) 
        fclustsummary.write("\n")
    fclustsummary.write("\n")
       
positions_contact =  [11,13,14,17] #(dans R 12, 14, 15, 18, -1 pour l'index
positions_to_exclide = [19]  #(dans R 20 -1 pour l'index
 
df = pd.read_csv(args.input, sep=';')
seqids = np.unique(df.SeqID)
# list of output file
fl = open(args.output_dir+"/list_of_files.txt", "w")
for seqid in seqids:
    print(seqid)
    extract_seqid = df[df["SeqID"] == seqid]
    contigs = np.unique(extract_seqid["Contig"])
    for contig in contigs :
    
        # prefixe des noms de fichiers
        file_name = seqid+"_"+contig
        
        extract_contig = extract_seqid[extract_seqid["Contig"] == contig]
        
        # liste des noms
        zf_names = []
        for name in extract_contig["ZF name"]:
            zf_names.append(name)
        
        # liste des proteines    
        proteins = extract_contig["uniformised ZF string"]
        
        # calcul du zfd sur toutes les sequences de proteines
        zfd_all = zfd(list(proteins))
        
        # creation d'un dico des proteines
        dico_protein = {}
        zf_number = 0
        for protein in proteins:
            protein_seq = Seq(protein)
            record = SeqRecord(protein_seq, id=zf_names[zf_number])
            dico_protein[zf_names[zf_number]] = record
            zf_number +=1
            
        # sequences dna
        sequences = extract_contig["dna sequence reading strand"]   
        # calcul du zfd sur toutes les sequences nucleiques
        zfd_codon_all = zfd_codon(list(sequences))           
        sequences2 = sequences # pour la boucle
        
        # creation d'un dico des noms de sequence dna, d'un dico des array,  
        # du fichier d'entrre pour silixx.
        fsilix=open(args.output_dir+"/"+file_name+".silix", "w")
        description = "Sequence :"+seqid+"; Contig:"+contig
        my_records = []
        zf_number = 0
        dico_sequence = {}
        list_array_total = {}
        nb_29 = 0
        nb_28 = 0
        for sequence in sequences:
            dna_seq = Seq(sequence)
            record = SeqRecord(dna_seq, id=zf_names[zf_number],description="Zinc finger ;"+description,)
            my_records.append(record)
            dico_sequence[zf_names[zf_number]] = record
            arr = zf_names[zf_number][0]
            if arr in list_array_total :
                list_array_total[arr].append(zf_names[zf_number])
            else :
                list_array_total[arr] = []
                list_array_total[arr].append(zf_names[zf_number])
            suff = zf_names[zf_number][-2:]
            if suff == "28":
                nb_28 +=1
            elif suff == "29":
                nb_29 +=1    
            else :
                sys.exit("Unknown suffix for "+zf_names[zf_number])
            zf_number2 = 0
            for sequence2 in sequences2:
                distance = distance_pos_3(sequence.upper(),sequence2.upper())
                if zf_number2 > zf_number:
                    if distance <= 2:
                        fsilix.write(zf_names[zf_number]+" "+zf_names[zf_number2]+"\n")
                zf_number2 += 1            
            zf_number +=1
        array_max = list(list_array_total.keys())[0]
        size_array_max = 0
        for array in list(list_array_total.keys()):
            if len(list_array_total[array]) > size_array_max:
                array_max = array
                size_array_max = len(list_array_total[array])
                
        # Calcul de la diversite sur le plus grand array
        print("For array max "+array_max)  
        proteins_array_max =  []
        dna_array_max =  []
        for name in  list_array_total[array_max]:
            print(name)
            protein_seq = dico_protein[name].seq
            print(protein_seq) 
            proteins_array_max.append(protein_seq)
            dna_seq = dico_sequence[name].seq
            dna_array_max.append(dna_seq)
        zfd_array = zfd(proteins_array_max)    
        zfd_codon_array = zfd_codon(dna_array_max)  
                  
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
                if len(buf) < 2:
                    continue
                fam = buf[0]
                seq = buf[1]
                if fam in families:
                    families[fam].append(seq)
                else:
                    families[fam] = []
                    families[fam].append(seq)
            clusters = list(families.keys())

            fclustsummary=open(args.output_dir+"/"+file_name+".clust_summary", "w") 
            fclustsummary.write("# Nb of zf : "+str(len(sequences))+"\n")   
            fclustsummary.write("\n")
            fclustsummary.write("#################################################\n")
            fclustsummary.write("# All ZF\n")
            fclustsummary.write("#################################################\n")
            write_divindex(zfd_all, zfd_codon_all, "Global")

            calcul_synonym_divindex(zfd_all, zfd_codon_all, "Global")
            
            write_divindex_pos(zfd_all, "Global")    
             
            write_divindex_pos(zfd_codon_all, "Global codons")       
      
                                
            fclustsummary.write("# Nb of zf 28 : "+str(nb_28)+"\n")
            fclustsummary.write("# Nb of zf 29 : "+str(nb_29)+"\n")                
            fclustsummary.write("# Nb of arrays : "+str(len(list(list_array_total.keys())))+"\n")
            fclustsummary.write("# Longest array : " + array_max + "\n")
            fclustsummary.write("# Size longest array : " + str(size_array_max)+ "\n") 
            
            fclustsummary.write("\n")
            fclustsummary.write("#################################################\n")
            fclustsummary.write("# Longest ZF array\n")
            fclustsummary.write("#################################################\n")
            
            write_divindex(zfd_array, zfd_codon_array, "Longest array")               
            
            calcul_synonym_divindex(zfd_array, zfd_codon_array, "Longest array")            
            
            write_divindex_pos(zfd_array, "Longest array")        

            write_divindex_pos(zfd_codon_array, "Longest array codons")  
            
                          
            if len(clusters) == 0: 
                fclustsummary.write("# Nb of clusters : 0\n")             
                fclustsummary.close()
            else:    
                cluster_max = clusters[0]
                nb_seq_max = len(families[cluster_max])
                new_records = []
                dico_sequence_bck  = dico_sequence.copy() # je garde une copie, car je vide le dico
                print("debug avant")
                print(dico_sequence_bck)
                for cluster in clusters:
                    seqs = families[cluster]
                    for seq in seqs:
                        #record = dico_sequence[seq]
                        record = dico_sequence.pop(seq)
                        dna_seq = record.seq
                        new_record = SeqRecord(dna_seq, id=seq+"_C"+cluster,description="Zinc finger ;"+description,)
                        new_records.append(new_record)
                    nb_seq = len(families[cluster])
                    print(cluster +": "+str(nb_seq))
                    if nb_seq > nb_seq_max :
                        cluster_max = cluster
                        nb_seq_max = nb_seq
                orphans = list(dico_sequence.keys())
                for orphan in orphans:
                    record = dico_sequence.pop(orphan)
                    dna_seq = record.seq
                    new_record = SeqRecord(dna_seq, id=orphan+"_singleton",description="Zinc finger ;"+description,)
                    new_records.append(new_record)             
                print("Nb of clusters : "+str(len(clusters)))
                print("Identifier cluster max : C"+cluster_max)
                print("Size cluster max : "+ str(nb_seq_max)) 
                print("Nb of singletons : "+ str(len(orphans))) 
                #print(dico_dna_clust)
                print(dico_sequence_bck)
                fclustsummary.write("# Nb of clusters : "+str(len(clusters))+"\n")
                fclustsummary.write("# Identifier cluster max : C_"+cluster_max+"\n")
                fclustsummary.write("# Size cluster max : "+ str(nb_seq_max)+"\n")
                fclustsummary.write("# Nb of singletons : "+ str(len(orphans))+"\n")
                seqs = families[cluster_max]
                fclustsummary.write("# Cluster max contents: ")
                list_array = []
                cluster_max_protein_sequence = []
                cluster_max_dna_sequence = []
                for seq in seqs:
                    cluster_max_protein_sequence.append(dico_protein[seq].seq)
                    cluster_max_dna_sequence.append(dico_sequence_bck[seq].seq)
                    fclustsummary.write(seq+"_C"+cluster_max+" ")
                    arr = seq[0]
                    if not (arr in list_array) :
                        list_array.append(arr)
                fclustsummary.write("\n")
                zfd_cluster = zfd(cluster_max_protein_sequence) 
                zfd_codon_cluster = zfd_codon(cluster_max_dna_sequence)  
                
                fclustsummary.write("\n")
                fclustsummary.write("#################################################\n")
                fclustsummary.write("# Longest ZF cluster (based on similarity at 3rd codon positions)\n")             
                fclustsummary.write("#################################################\n")
                                   
                write_divindex(zfd_cluster, zfd_codon_cluster, "Longest cluster")

                calcul_synonym_divindex(zfd_cluster, zfd_codon_cluster, "Longest cluster") 
                            
                write_divindex_pos(zfd_cluster, "Longest cluster")
                
                write_divindex_pos(zfd_codon_cluster, "Longest cluster codons")                    
                
                fclustsummary.write("# Nb of arrays in cluster max: "+str(len(list_array))+"\n") 
                fclustsummary.close() 
                count = SeqIO.write(new_records, args.output_dir+"/"+file_name+"_cluster.fasta", "fasta")    
fl.close()
