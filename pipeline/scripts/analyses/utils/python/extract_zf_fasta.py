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
    for i in range(0,len(zf)):
        if i not in positions_to_exclude:
            total = total + zf[i]
    if total == 0 :
        if sum > 0 :
              sys.exit("Error in zfsum")    
    else :
        sum = sum / total
    return sum
 
def zfmean(zf:[], positions:[]):
    sum = 0
    for pos in positions:
        sum = sum + zf[pos]
    mean = sum / len(position)
    return sum
    
# Calculate the zfd
def zfd(seqs:[]):
    zfd = []
    l = len(seqs[0])
    nbseq =  len(seqs)
    for i in range(0,l) :
        site = []
        for seq in seqs:
            site.append(seq[i]) 
        df = pd.DataFrame({'aa':site})
        unic_site = list(np.unique(df.aa))
        div = 0
        for aa in unic_site:
            nb_aa = site.count(aa)
            div  = div + (nb_aa / nbseq) ** 2
        div = 1 - div
        div = round(div, 5)
        zfd.append(div)
    return zfd


# Calculate the zfd on codon
def zfd_codon(seqs:[]):
    zfd = []
    l = len(seqs[0])
    nbseq =  len(seqs)
    for i in range(0,l,3) :
        site = []
        for seq in seqs:
            codon = seq[i]+seq[i+1]+seq[i+2]
            codon = codon.upper()
            site.append(codon)
        df = pd.DataFrame({'codon':site})
        unic_site = list(np.unique(df.codon))
        div = 0
        for codon in unic_site:
            nb_codon = site.count(codon)
            div  = div + (nb_codon / nbseq) ** 2
        div = 1 - div
        div = round(div, 5)
        zfd.append(div)
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
    fclustsummary.write('{:50}'.format("# Position : "))
    ipos = 1 
    for div in zfd_data:
        fclustsummary.write('%6d' % (ipos))
        ipos += 1       
    fclustsummary.write("\n")      
          
    fclustsummary.write('{:50}'.format("# " + zfd_name + " div. index aa : "))
    for div in zfd_data:
        fclustsummary.write('%6.2f' % (div))
    fclustsummary.write("\n")
    
    fclustsummary.write('{:50}'.format("# " + zfd_name + " div. index  codons : "))
    for div in zfd_codon_data:
        fclustsummary.write('%6.2f' % (div))
    fclustsummary.write("\n")    
    
    fclustsummary.write('{:50}'.format("# Legend : "))
    ipos = 0 
    legend_loc = legend.copy()
    legend_cons_loc = legend_cons.copy()    
    for div in zfd_data:
        leg = "."
        if ipos in positions_contact:
           leg = legend_loc.pop()
        if ipos in conserved_positions:
           leg = legend_cons_loc.pop()
        if ipos in positions_to_exclude:
           leg = "X"
        fclustsummary.write('{:>6}'.format(leg))
        ipos += 1       
    fclustsummary.write("\n") 

    
def write_divindex_pos(zfd_data,zfd_name):
            zfd_pos = zfsum(zfd_data,positions_contact )
            pos_string = ' '.join(str(item+1) for item in positions_contact)
            pos_string = ' '.join(str(item) for item in reversed(legend))
            #fclustsummary.write(zfd_name + " d. i. for " + pos_string)
            fclustsummary.write('{:50}'.format("# " + zfd_name + " ZFD " + pos_string+ " : "))
            fclustsummary.write('%6.3f' % (zfd_pos))
            fclustsummary.write("\n")
            df_csv.loc[df_csv['ZF_Dataset'] == zfd_name, "ZFD"] = round(zfd_pos, 2)

def calcul_synonym_divindex(zfd_data, zfd_codon_data, zfd_name):    
    zdf_syno = []
    zdf_ratio = []
    total_syno = 0.0
    total_ratio = 0.0
    total_ratio_no_pos = 0.0 # on exclut les position d'interet
    total_ratio_pos = 0.0 # on garde les position d'interet
    for i in range(0,len(zfd_data)):
        syno = zfd_codon_data[i] - zfd_data[i]
        if syno < 0:
            print(zfd_codon_data[i])
            print(zfd_data[i])
            sys.exit("Error in calcul_synonym_divindex, negative value")
        
        zdf_syno.append(syno)   
        if i not in positions_to_exclude : 
            total_syno += syno        
            
    mean_syno = total_syno / (len(zfd_data)   - len(positions_to_exclude)) 

    if mean_syno == 0 :
        mean_syno = 0.0000001

    for i in range(0,len(zfd_data)):        
        ratio = zfd_data[i] / mean_syno
        if ratio > 99 :
            ratio = 99.0
        zdf_ratio.append(ratio) 
        if i not in positions_to_exclude :       
            total_ratio += ratio
            if i not in positions_contact :
                total_ratio_no_pos += ratio
            else :
                total_ratio_pos += ratio

    mean_ratio = total_ratio / (len(zfd_data)   - len(positions_to_exclude))    
    mean_ratio_no_pos = total_ratio_no_pos / (len(zfd_data)   - len(positions_to_exclude) -len(positions_contact) )    
    mean_ratio_pos = total_ratio_pos / len(positions_contact)
    fclustsummary.write('{:50}'.format("# " + zfd_name + " syno (codon - aa) : "))
    for div in zdf_syno:
        fclustsummary.write('%6.2f' % (div))
    fclustsummary.write("\n") 
    fclustsummary.write('{:50}'.format("# "+ zfd_name + " ratio (aa / mean syno) : "))
    for div in zdf_ratio:
        fclustsummary.write('%6.2f' % (div))
    fclustsummary.write("\n")       
    fclustsummary.write('{:50}'.format("# Mean syno : "))     
    fclustsummary.write('%5.2f' % (mean_syno))       
    fclustsummary.write("\n") 
    fclustsummary.write('{:50}'.format("# Mean ratio : "))     
    fclustsummary.write('%5.2f' % (mean_ratio))       
    fclustsummary.write("\n")
    
    
    df_csv.loc[df_csv['ZF_Dataset'] == zfd_name, "Mean_divS"] = round(mean_syno, 2)
    df_csv.loc[df_csv['ZF_Dataset'] == zfd_name, "Mean_ratio_AS"] = round(mean_ratio, 2)
    jpos = 3
    for pos in positions_contact:
        fclustsummary.write('{:50}'.format("# " + zfd_name + " position "+str(pos+1)))
        fclustsummary.write(" Ratio : ")     
        fclustsummary.write('%5.2f' % (zdf_ratio[pos])) 
        fclustsummary.write(" Mean ratio excluding contact positions :")    
        fclustsummary.write('%5.2f' % (mean_ratio_no_pos)) 
        fclustsummary.write("\n")
        as_name = "Ratio_AS"+legend[jpos]
        df_csv.loc[df_csv['ZF_Dataset'] == zfd_name, as_name] = round(zdf_ratio[pos], 2)
        jpos = jpos -1
    df_csv.loc[df_csv['ZF_Dataset'] == zfd_name, "Mean_ratio_AS_non-1236"] = round(mean_ratio_no_pos, 2)
    df_csv.loc[df_csv['ZF_Dataset'] == zfd_name, "Mean_ratio_AS-1236"] = round(mean_ratio_pos, 2)    
    #fclustsummary.write("\n")
    
       
positions_contact =  [11,13,14,17] #(dans R 12, 14, 15, 18, -1 pour l'index
legend = ["-1", "2", "3", "6"]
legend.reverse()
positions_to_exclude = [19]  #(dans R 20 -1 pour l'index
conserved_positions = [2,5,18]
legend_cons = ["C", "C", "H"]
legend_cons.reverse()
min_ZF_nb = 4

df = pd.read_csv(args.input, sep=';')
seqids = np.unique(df.SeqID)
# list of output file
fl = open(args.output_dir+"/list_of_files.txt", "w")
for seqid in seqids:
    print(seqid)
    extract_seqid = df[df["SeqID"] == seqid]
    contigs = np.unique(extract_seqid["Contig"])
    for contig in contigs :
    
        # Creation du df pour le csv

        df_csv = pd.DataFrame(columns=["SeqID","Status","ZF_Dataset","ID_set", "Nb_ZF","Nb_ZF_28","Nb_ZF_29","Nb_Arrays","Nb_Clusters","Nb_singletons","Mean_divS", "Mean_ratio_AS","Ratio_AS-1","Ratio_AS2","Ratio_AS3","Ratio_AS6","Mean_ratio_AS-1236" ,"Mean_ratio_AS_non-1236" ,"ZFD"])
        new_row = pd.DataFrame({"ZF_Dataset" : "All_ZFs" }, index=[0])
        df_csv = pd.concat([df_csv.loc[:],new_row]).reset_index(drop=True)

        new_row = pd.DataFrame({"ZF_Dataset" : "Longest_ZF_array" }, index=[0])
        df_csv = pd.concat([df_csv.loc[:],new_row]).reset_index(drop=True)
        
        new_row = pd.DataFrame({"ZF_Dataset" : "Largest_ZF_cluster" }, index=[0])
        df_csv = pd.concat([df_csv.loc[:],new_row]).reset_index(drop=True)
          
          
        df_csv["SeqID"] = seqid    
        # prefixe des noms de fichiers
        file_name = seqid+"_"+contig
        
        extract_contig_all = extract_seqid[extract_seqid["Contig"] == contig]
        
        extract_contig = extract_contig_all[extract_contig_all["Status"] == "Ok"]
        
        # liste des noms
        zf_names = []
        for name in extract_contig["ZF name"]:
            zf_names.append(name)
        
        # liste des proteines    
        proteins = extract_contig["uniformised ZF string"]
        
        # calcul du zfd sur toutes les sequences de proteines
        if len(list(proteins)) > 0 :
            df_csv.loc[df_csv['SeqID'] == seqid, "Status"] = 0
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
            if len(list(sequences)) > 0 :
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
                            
                df_csv.loc[df_csv['ZF_Dataset'] == "All_ZFs", "Nb_ZF_28"] = nb_28
                df_csv.loc[df_csv['ZF_Dataset'] == "All_ZFs", "Nb_ZF_29"] = nb_29
            
            # Calcul de la diversite sur le plus grand array
            proteins_array_max =  []
            dna_array_max =  []
            nb_28 = 0
            nb_29 = 0
            for name in  list_array_total[array_max]:
                protein_seq = dico_protein[name].seq
                proteins_array_max.append(protein_seq)
                dna_seq = dico_sequence[name].seq
                dna_array_max.append(dna_seq)
                suff = name[-2:]
                if suff == "28":
                    nb_28 +=1
                elif suff == "29":
                    nb_29 +=1  
            df_csv.loc[df_csv['ZF_Dataset'] == "Longest_ZF_array", "Nb_ZF_28"] = nb_28
            df_csv.loc[df_csv['ZF_Dataset'] == "Longest_ZF_array", "Nb_ZF_29"] = nb_29                
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
                fclustsummary.write("#\n")
                fclustsummary.write("#################################################\n")
                fclustsummary.write("# All ZF\n")
                fclustsummary.write("#################################################\n")
                
                if zf_number >= min_ZF_nb:
                
                    write_divindex(zfd_all, zfd_codon_all, "All_ZFs")

                    calcul_synonym_divindex(zfd_all, zfd_codon_all, "All_ZFs")
                
                    write_divindex_pos(zfd_all, "All_ZFs")    
                
                    write_divindex_pos(zfd_codon_all, "All_ZFs codons")       
        
                else :
                
                    fclustsummary.write("# Less than 4 zf:  no stats.\n")  
                                            
                fclustsummary.write("# Nb of zf 28 : "+str(nb_28)+"\n")
                fclustsummary.write("# Nb of zf 29 : "+str(nb_29)+"\n")                
                fclustsummary.write("# Nb of arrays : "+str(len(list(list_array_total.keys())))+"\n")
                fclustsummary.write("# Longest array : " + array_max + "\n")
                fclustsummary.write("# Size longest array : " + str(size_array_max)+ "\n") 
                
                df_csv.loc[df_csv['ZF_Dataset'] == "All_ZFs", "Nb_ZF"] = zf_number
                df_csv.loc[df_csv['ZF_Dataset'] == "Longest_ZF_array", "ID_set"] = array_max
                df_csv.loc[df_csv['ZF_Dataset'] == "All_ZFs", "Nb_Arrays"] = len(list(list_array_total.keys()))
                
                fclustsummary.write("#\n")
                fclustsummary.write("#################################################\n")
                fclustsummary.write("# Longest ZF array\n")
                fclustsummary.write("#################################################\n")
                
                df_csv.loc[df_csv['ZF_Dataset'] == "Longest_ZF_array", "Nb_ZF"] = size_array_max
                if size_array_max >= min_ZF_nb :
                    write_divindex(zfd_array, zfd_codon_array, "Longest_ZF_array")               
                
                    calcul_synonym_divindex(zfd_array, zfd_codon_array, "Longest_ZF_array")            
                
                    write_divindex_pos(zfd_array, "Longest_ZF_array")        

                    write_divindex_pos(zfd_codon_array, "Longest_ZF_array codons")  
                
                else :
                    fclustsummary.write("# Less than 4 zf:  no stats.\n")              
                
                df_csv.loc[df_csv['ZF_Dataset'] == "All_ZFs", "Nb_Clusters"] = len(clusters)
                if len(clusters) == 0: 
                    fclustsummary.write("# Nb of clusters : 0\n") 
                    cluster_max = "None"            
                    #fclustsummary.close()
                else:    
                    cluster_max = clusters[0]
                    nb_seq_max = len(families[cluster_max])
                    new_records = []
                    dico_sequence_bck  = dico_sequence.copy() # je garde une copie, car je vide le dico
                    for cluster in clusters:
                        seqs = families[cluster]
                        for seq in seqs:
                            #record = dico_sequence[seq]
                            record = dico_sequence.pop(seq)
                            dna_seq = record.seq
                            new_record = SeqRecord(dna_seq, id=seq+"_C"+cluster,description="Zinc finger ;"+description,)
                            new_records.append(new_record)
                        nb_seq = len(families[cluster])
                        if nb_seq > nb_seq_max :
                            cluster_max = cluster
                            nb_seq_max = nb_seq
                    orphans = list(dico_sequence.keys())
                    for orphan in orphans:
                        record = dico_sequence.pop(orphan)
                        dna_seq = record.seq
                        new_record = SeqRecord(dna_seq, id=orphan+"_singleton",description="Zinc finger ;"+description,)
                        new_records.append(new_record)             
    
                    fclustsummary.write("# Nb of clusters : "+str(len(clusters))+"\n")
                    fclustsummary.write("# Identifier largest cluster : C_"+cluster_max+"\n")
                    fclustsummary.write("# Size largest cluster : "+ str(nb_seq_max)+"\n")
                    fclustsummary.write("# Nb of singletons : "+ str(len(orphans))+"\n")
                    seqs = families[cluster_max]
                    fclustsummary.write("# Cluster max contents: ")
                    
                    df_csv.loc[df_csv['ZF_Dataset'] == "Largest_ZF_cluster", "ID_set"] = "C_"+cluster_max
                    df_csv.loc[df_csv['ZF_Dataset'] == "Largest_ZF_cluster", "Nb_ZF"] = nb_seq_max
                    df_csv.loc[df_csv['ZF_Dataset'] == "All_ZFs", "Nb_singletons"] = len(orphans)
                    
                    
                    
                    list_array = []
                    cluster_max_protein_sequence = []
                    cluster_max_dna_sequence = []
                    nb_28 = 0
                    nb_29 = 0
                    for seq in seqs:
                        suff = seq[-2:]
                        if suff == "28":
                            nb_28 +=1
                        elif suff == "29":
                            nb_29 +=1  
                        cluster_max_protein_sequence.append(dico_protein[seq].seq)
                        cluster_max_dna_sequence.append(dico_sequence_bck[seq].seq)
                        fclustsummary.write(seq+"_C"+cluster_max+" ")
                        arr = seq[0]
                        if not (arr in list_array) :
                            list_array.append(arr)
                            
                    df_csv.loc[df_csv['ZF_Dataset'] == "Largest_ZF_cluster", "Nb_ZF_28"] = nb_28
                    df_csv.loc[df_csv['ZF_Dataset'] == "Largest_ZF_cluster", "Nb_ZF_29"] = nb_29                          
                    fclustsummary.write("\n")
                    zfd_cluster = zfd(cluster_max_protein_sequence) 
                    zfd_codon_cluster = zfd_codon(cluster_max_dna_sequence)  
                    
                    fclustsummary.write("#\n")
                    fclustsummary.write("#################################################\n")
                    fclustsummary.write("# Largest ZF cluster (based on similarity at 3rd codon positions)\n")             
                    fclustsummary.write("#################################################\n")
                    
                    if nb_seq_max >= min_ZF_nb :              
                        write_divindex(zfd_cluster, zfd_codon_cluster, "Largest_ZF_cluster")

                        calcul_synonym_divindex(zfd_cluster, zfd_codon_cluster, "Largest_ZF_cluster") 
                                
                        write_divindex_pos(zfd_cluster, "Largest_ZF_cluster")
                    
                        write_divindex_pos(zfd_codon_cluster, "Largest_ZF_cluster codons")                    
                    
                    else :
                        fclustsummary.write("# Less than 4 zf:  no stats.\n")    
                    
                    fclustsummary.write("# Nb of arrays in largest cluster : "+str(len(list_array))+"\n")
                    count = SeqIO.write(new_records, args.output_dir+"/"+file_name+"_cluster.fasta", "fasta")     
                df_csv.to_csv(fclustsummary,sep=';' , index=False , na_rep="NA")
                fclustsummary.close()



        extract_contig_errors = extract_contig_all[extract_contig_all["Status"] != "Ok"]
        erroneous_protein_names = extract_contig_errors["SeqID"]
        for erroneous_protein_name in erroneous_protein_names :
            print("Problem occurs with " + erroneous_protein_name)
            fclustsummary=open(args.output_dir+"/"+file_name+".clust_summary", "w")    
            fclustsummary.write("#################################################\n")
            fclustsummary.write("# Problems\n")
            err_status = list(extract_contig_all[extract_contig_all["SeqID"] == erroneous_protein_names]["Status"])
            fclustsummary.write("# "+str(err_status) +"\n")
            fclustsummary.write("# Status 1 : No exon.\n")
            fclustsummary.write("# Status 2 : Protein from fasta input file is different from translated DNA equence.\n")
            fclustsummary.write("# Status 3 : Other.\n")             
            fclustsummary.write("#################################################\n")
            status = 3
            if str(err_status) == "['no exon']" :
                status = 1
            if str(err_status) == "['sequence transl. problem']" :
                status = 2

            df_csv.loc[df_csv['SeqID'] == seqid, "Status"] = status
            df_csv.to_csv(fclustsummary,sep=';' , index=False , na_rep="NA")
            fclustsummary.close() 
fl.close()
