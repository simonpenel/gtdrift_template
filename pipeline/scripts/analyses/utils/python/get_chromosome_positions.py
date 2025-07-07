import argparse
import time
import os
import pandas as pd
from Bio import SeqIO
 
parser = argparse.ArgumentParser()

parser.add_argument('-i', '--input', type=str, required=True, help='whole_summary.csv from protein analysis')

parser.add_argument('-p', '--protein', type=str, required=True, help='path to fasta protein file')

parser.add_argument('-g', '--gff', type=str,required=True, help='path to gff file, REQUIRED if type=prot')
parser.add_argument('-o', '--output', type=str, required=True, help='output file path')

args = parser.parse_args()

## Function to read gff files and get chromosome,
# start and end for locus.

def preProcessGff(gff:str):
    with open(gff, 'r') as reader:
        print("Preprocessing gff... (This may take several minutes)")        
        print("Reading gff...")
        for line in reader:
            if line.startswith('#'):
                continue
            if line.split('\t')[2] == 'gene' or line.split('\t')[2] == 'pseudogene':
            
                splitline = line.split('\t')
                chrom     = splitline[0]
                strand    = splitline[6]
                genepos   = [splitline[3], splitline[4]]
                gene_info   = splitline[8]
                pseudo = 0
                if line.split('\t')[2] == 'pseudogene':
                    pseudo = 1
                #identifiant = gene_info.split(';')[0].split("=")
                #if identifiant[0] == "ID" :
                #        gene_id = identifiant[1]
                #        dico_gene[gene_id] = [chrom,strand,genepos]   
                flag = 0                     
                xrefs = gene_info.split(';')[1].split(",")
                for xref in xrefs:
                    ref=xref.split(':')
                    if ref[0] == "Dbxref=GeneID" or ref[0] == "GeneID":
                        gene_id = ref[1]
                        dico_gene[gene_id] = [chrom,strand,genepos,pseudo]
                        flag = 1
                if flag == 0 :
                    xrefs = gene_info.split(';')
                    for xref in xrefs:
                        ref=xref.split('=')
                        if ref[0] == "locus_tag" :
                            gene_id = ref[1].rstrip()
                            dico_gene[gene_id] = [chrom,strand,genepos,pseudo]
                            flag = 1
                if flag == 0 :
                    print("debug2 "+gene_info)
                    
def processGff(gff:str):
    with open(gff, 'r') as reader:
        print("Processing gff... (This may take several minutes)")        
        print("Reading gff...")
        for line in reader:
            if line.startswith('#'):
                continue
            if line.split('\t')[2] == 'CDS':
                cds_gene_id = "none"
                prot_name = "none"
                splitline = line.split('\t')
                cds_info   = splitline[8]
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
                            
                if flag_gene == 0 :
                    print("Warning : No gene")    
                if flag_prot == 0 :
                    print("Warning : No protein")                      
                if cds_gene_id == "none":
                    print("WARNING: Unable to get gene id : ")
                    print(cds_info)
                if prot_name == "none":
                    print("WARNING: Unable to get protein name : ")
                    print(cds_info)
                if prot_name != "none" and cds_gene_id != "none":
                    if not prot_name in dico_prot:
                        if  cds_gene_id in dico_gene:
                            dico_prot[prot_name] = dico_gene[cds_gene_id]
                        else :
                            print("ERROR "+cds_gene_id + "(" + cds_info+ ")")    


record_dict = SeqIO.to_dict(SeqIO.parse(args.protein, "fasta"))
dico_gene = {}                                                              
preProcessGff(args.gff)
dico_prot = {}     
processGff(args.gff)
df = pd.DataFrame() 
with open(args.input, 'r') as reader:
    lines = reader.readlines()[1:]
    for line in lines:
        splitline = line.split(';')
        prot_name = splitline[1]
        if prot_name in dico_prot :
            chromo_info = dico_prot[prot_name]
        else :
            print("Warning: No chromosome information for " + prot_name)
            chromo = ["NA", "NA", [0,0], 0]
        chromo = chromo_info[0]
        strand = chromo_info[1]
        positions = chromo_info[2]
        pseudo = chromo_info[3]
        from_fasta = record_dict[prot_name]
        protein_length = len(from_fasta.seq)
        new_row = {"SeqID": prot_name, "Chromosome": chromo, "Chr Start": positions[0],"Chr End":  positions[1], "Strand" : strand, "Protein Length" : protein_length, "Pseudo" : pseudo}
        df = df._append(new_row, ignore_index=True)
        

df_summary = pd.read_csv(args.input, sep=';', header=0)
df_summary = df_summary.drop(["Unnamed: 0"],axis=1) 

if len(df.index) > 0 :
    newdf = df_summary.join(df.set_index('SeqID'),on='SeqID', how="outer",lsuffix='_caller', rsuffix='_other')
else:
    newdf = df_summary
    newdf['Chromosome'] = []
    newdf['Chr Start'] = []
    newdf['Chr End'] = []
    newdf['Strand'] = []
    newdf['Protein Length'] = []
    newdf['Pseudo'] = []
newdf.to_csv(args.output, sep=';')




