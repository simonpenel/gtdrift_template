import argparse
import time
import os
 
parser = argparse.ArgumentParser()

parser.add_argument('-i', '--input', type=str, required=True, help='blast output directory path')
parser.add_argument('-o', '--output', type=str, required=True, help='output file path')
parser.add_argument('-t', '--type', type=str, default='nucl', help='type of blast [nucl, prot], default: nucl')
parser.add_argument('-g', '--gff', type=str, help='path to gff file, REQUIRED if type=prot')

args = parser.parse_args()

## Function to read gff files and get chromosome,
# start and end for locus.

## WARNING: Reading the gff takes a while.
## This code can be HEAVILY optimized.
def preProcessGff(gff:str, out:str, names:list):
    with open(gff, 'r') as reader, open(out, 'w') as writer:
        print("Preprocessing gff... (This may take several minutes)")
        totstart = time.time()
        gene     = {}
        genelist = []
        donelist = []

        protparse = []
        
        print("Reading gff...")
        timestart = time.time()
        for line in reader:
            if line.startswith('#'):
                continue
            if line.split('\t')[2] == 'gene':
                splitline = line.split('\t')
                gene_id   = splitline[8].split(';')[0].split('=')[1].split('-')[1].strip()
                chrom     = splitline[0].strip()
                strand    = splitline[6].strip()
                genepos   = [splitline[3].strip(), splitline[4].strip()]
                gene_ref  = splitline[1].strip()
                if '%2C' in gene_ref:
                    gene_ref.replace('%2C', ',')
                if not gene_id in gene.keys():
                    gene[gene_id] = []
                    gene[gene_id].append([strand, genepos, chrom, gene_ref])
                else:
                    gene[gene_id].append([strand, genepos, chrom, gene_ref])
                #print("DEBUG3 ",gene_id,gene[gene_id])
                genelist.append(gene_id)
            elif line.split('\t')[2] == 'CDS':
                protparse.append(line)
        timeend = time.time()
        print(f'Done reading gff. Parsing time: {round(timeend-timestart, 2)} seconds / {round((timeend-timestart)/60, 2)} minutes\n')
        #for toto in genelist:
        #    print("DEBUG " +toto+" ",gene[toto])
        #print("DEBUG",protparse)    
        print("Getting all proteins...")
        timestart = time.time()
        protparse.sort(key=lambda x: (x.split('\t')[0], x.split('\t')[3]))
        intpos      = {}

        
        for p in protparse:
           
            splitline   = p.split('\t')
            parent_gene = splitline[8].split(';')[5].split('=')[1] 
            prot_id     = splitline[8].split(';')[0].split('=')[1].split('-')[1]
            if parent_gene in genelist and prot_id in names:
                print("DEBUG4 "+parent_gene +" "+prot_id+ " "+p)
                if not prot_id in intpos.keys():
                    intpos[prot_id] = []
            
                    intpos[prot_id].append(splitline[4]) 
                    int_start = splitline[3]
                    

                    i = 0
                    print("DEBUG5 "+prot_id+ " "+str(len(gene[parent_gene])))
                    #??
                    while i < len(gene[parent_gene]):
                        if (gene[parent_gene][i][1][0] <= splitline[3] <= gene[parent_gene][i][1][1]) and (splitline[0] == gene[parent_gene][i][2]):
                            break
                        i += 1

                    if not [prot_id, splitline[3], splitline[4]] in donelist and i < len(gene[parent_gene]):
                        writer.write(f'{gene[parent_gene][i][2]}\t{gene[parent_gene][i][1][0]}\t{gene[parent_gene][i][1][1]}')
                        writer.write(f'\t{gene[parent_gene][i][0]}\t{gene[parent_gene][i][3]}\tNA\t{prot_id}\t{int_start}\n')
                        donelist.append([prot_id, splitline[3], splitline[4]])
                        print(f'DEBUG6 {gene[parent_gene][i][2]}\t{gene[parent_gene][i][1][0]}\t{gene[parent_gene][i][1][1]}\t{gene[parent_gene][i][0]}\t{gene[parent_gene][i][3]}\tNA\t{prot_id}\t{int_start}')
                        

                else:
                    intpos[prot_id].append(splitline[3])
                    int_start = splitline[3]
                    
                    i = 0
                    while i < len(gene[parent_gene]):
                        if gene[parent_gene][i][1][0] <= splitline[3] <= gene[parent_gene][i][1][1]:
                            break
                        i += 1
                    
                    if not [prot_id, splitline[3], splitline[4]] in donelist and i < len(gene[parent_gene]):
                        writer.write(f'{gene[parent_gene][i][2]}\t{gene[parent_gene][i][1][0]}\t{gene[parent_gene][i][1][1]}')
                        writer.write(f'\t{gene[parent_gene][i][0]}\t{gene[parent_gene][i][3]}\t{intpos[prot_id][0]}-{intpos[prot_id][1]}\t{prot_id}\t{int_start}\n')
                        donelist.append([prot_id, splitline[3], splitline[4]])

                    intpos[prot_id] = []
                    intpos[prot_id].append(splitline[4])
    timeend = time.time()
    print(f"Proteins done. Parsing time: {round(timeend-timestart, 2)} seconds / {round((timeend-timestart)/60, 2)} minutes\n\n")
    print(f"Preprocessing done. Total time: {round(timeend-totstart, 2)} seconds / {round((timeend-totstart)/60, 2)} minutes")

def getPosGff(hits:list, gff:str):
        preprocess = args.output+'.temporary'
        buffer   = {}
        output   = []
        names    = set()
        strnames = set()
        print(f'\n{len(hits)} hits to parse (this may take a while)')
        unique = [list(x) for x in set(tuple(x) for x in hits)]
        unique.sort(key=lambda x: (x[0], x[1]))
        for hit in unique:
            print("GETPOS  DEBUG ", hit)
        print("Parsing through names...")
        for elt in unique:
                names.add(elt[0])
                strnames.add('Name='+elt[0])
        print("Done parsing\n")

        preProcessGff(gff, preprocess, names)

        print("Getting all hits...")
        timestart = time.time()
        prev_buffer = {}
        for obj in unique:
            print("GETPOS DEBUG2 ", obj)
            exon = obj[4]
            name = obj[0]
            info = []

            index_file = open(preprocess, 'r')
            for hit in index_file:
                if hit.split('\t')[6].strip() == name:
                    info.append(hit.strip().split('\t'))
                    print("GETPOS  DEBUG3 ", hit)
            index_file.close()

            if info != []:
                chr    = info[0][0]
                strand = info[0][3]
                ref    = info[0][4]
                prot   = info[0][6]
                intron = []
                for line in info:
                    intron.append(line[5])
                if strand == '+':
                    start  = info[0][1]
                    end    = info[0][2]
                elif strand == '-':
                    start  = info[0][2]
                    end    = info[0][1]
                
                start_int = info[0][7]
                    
                if not prot in buffer.keys():

                    if buffer != {}:
                        out_elt = prev_buffer
                        #print("DEBUG8 ",out_elt)
                        output.append(out_elt)

                    buffer[prot] = [chr, start, end, strand, exon, prot, intron, ref, start_int]
                    prev_buffer  = buffer[prot]
                else:
                    prev_buffer[2]  = end
                    if len(prev_buffer[6]) < len(intron):
                        prev_buffer[6] = (intron)
                    if int(start_int) < int(prev_buffer[8]):
                        prev_buffer[8] = start_int
        out_elt = prev_buffer
        #print("DEBUG8 ",out_elt)
        output.append(out_elt)
        timeend = time.time()
        print(f'Done. Total time: {round(timeend-timestart, 2)} seconds / {round((timeend-timestart)/60, 2)} minutes')
        #os.remove(preprocess)
        return output


## Function to get a list of all hits from file.
# The file has to be a blast file.

def getHits(path2file:str) -> list:
# We make a list of all hits split by tabulation
    hits = []
    with open(path2file) as f:
        print(f"Opening file {path2file}")
        l = f.readline()

        while l:
# In the blast file, hits are the only lines
# to not start with a pound character ("#")
# Therefore we extract all lines not starting with one
            if l[0] != '#':

                hit  = l.split('\t')
                if "|" in hit[2]:
                    name = hit[2].split('|')[1].strip()
                else:
                    name = hit[2]

# Must make sure to use gff annotation for prot

                current = [name,hit[11],hit[12]]
                
# Our hits will have chr, start, end and strand

                if hit[11] < hit[12]:
                    current.append("+")
                else:
                    current.append("-")
                hits.append(current)
            l = f.readline()
    return hits


## Extra function to add exon number to hits.
# Exon number is passed as the num parameter.
# While a bit redundant, it is more practical.
def exonRead(path2file:str, num:int) -> list:
    hits = getHits(path2file)
    for hit in hits:
        hit.append(num)
    return hits


## Hits from different strands must be handled separately.
# Therefore, we make a function to separate hits based on strand.
def strandSeparate(hits:list) -> list:
    plus = []
    minus = []
    for hit in hits:
        if len(hit) > 0 :
            print(hit)
            if hit[3] == "+":
                plus.append(hit)
            elif hit[3] == "-":
                minus.append(hit)
    
# Sort all hits by chromosome and start position.
    plus.sort(key=lambda x: (x[0], x[1]))
    minus.sort(key=lambda x: (x[0], -int(x[1])))
    return [plus, minus]


## Function to define each potential locus.
def exonCount(hits:list):
    plus = hits[0]
    minus = hits[1]
    i = 0
# Starting parameters for our loci.
# We arbitrarily start with + strand.
    chr   = plus[0][0]
    start = plus[0][1]
    end   = plus[0][2]
    counted = [["locus"+str(i)+'_'+args.type, chr, start, end, "+",
               0,0,0,0,0,0,0,0,0,0]]
    counted[0][int(plus[0][4])+3] += 1
    if args.type == 'prot':
        counted[0].append(plus[0][5])
        #print("debug  test exon ",plus[0][6])
        #if plus[0][6] != ['NA']:
        #    counted[0].append(plus[0][6])
        counted[0].append(plus[0][6])    
        counted[0].append(plus[0][7])
        counted[0].append(plus[0][8])

# The 10 zeros in the list are for each exon counts,
# from exon2 to exon11.

    for hit in plus:
# We consider that the maximum distance between first and
# last exon should be 100k bases.
# We also check for same chromosome.
        print("DEBUG 9 ",hit)
        if int(hit[1]) < (int(start) + 100000) and int(hit[2]) > int(start) and hit[0] == chr:
            if args.type == 'nucl' or (args.type == 'prot' and hit[5] == counted[i][-4]):
# Increment correct exon.
                counted[i][int(hit[4])+3] += 1
                counted[i][3] = hit[2]
                if args.type == 'prot':
                    #print("debug  test exon ",plus[0][6])
                    if int(counted[i][-1]) > int(hit[8]):
                        counted[i][-1] = hit[8]
                    #if hit[6] != ['NA'] and not hit[6] in counted[i][-3]:
                    if  not hit[6] in counted[i][-3]:                    
                        counted[i][-3] += hit[6]
                        counted[i][-3] = list(set(counted[i][-3]))
        else:
# If out of range/chromosome, redefine starting parameters
# for next locus.
            start = hit[1]
            end = hit[2]
            chr = hit[0]
            i += 1
            counted.append(["locus"+str(i)+'_'+args.type, chr, start, end, "+",
                            0,0,0,0,0,0,0,0,0,0])
            counted[i][int(hit[4])+3] += 1
            if args.type == 'prot':
                counted[i].append(hit[5])
#                if hit[6] != ['NA']:
#                    counted[i].append(hit[6])
                counted[i].append(hit[6])    
                counted[i].append(hit[7])
                counted[i].append(hit[8])
# Repeat for - strand.

    i+=1
    chr = minus[0][0]
    start = minus[0][1]
    end = minus[0][2]
    print("DEBUG BORDEL ",minus[0])
    counted.append(["locus"+str(i)+'_'+args.type, chr, start, end, "-",
               0,0,0,0,0,0,0,0,0,0])
    counted[i][int(minus[0][4])+3] += 1
    if args.type == 'prot':
        counted[i].append(minus[0][5])
#        if minus[0][6] != ['NA']:
#           counted[i].append(minus[0][6])
        counted[i].append(minus[0][6])           
        counted[i].append(minus[0][7])
        counted[i].append(minus[0][8])
        print("DEBUG BORDEL "+str(i)+" ",counted[i])

	#EXPLORE gene:hits blast:brin:matching exon:protif:introns
    for hit in minus:
        print("DEBUG 9 ",hit)
        print("DEBUG 9 TEST "+str(hit[1]) +" > "+ str(int(start) - 100000) +" and "+str(hit[2]) + " < "+  str(start) + " and " + hit[0] +  " = "+chr)
        if int(hit[1]) > (int(start) - 100000) and int(hit[2]) < int(start) and hit[0] == chr:
            if args.type == 'prot' :
                print("DEBUG 9 COMPARE "+str(i)+" : ",hit," WITH   ",counted[i] )
                print("DEBUG 9 SELECTED TEST ",hit," ", hit[5]," = ",counted[i][-4] )
            if args.type == 'nucl' or (args.type == 'prot' and hit[5] == counted[i][-4]):
                print("DEBUG 9 SELECTED ",hit)
                print("DEBUG BORDEL COUNTED AVANT"+str(i)+" ",counted[i])
                counted[i][int(hit[4])+3] += 1
                counted[i][3] = hit[2]
                if args.type == 'prot':
                    if int(counted[i][-1]) > int(hit[8]):
                        counted[i][-1] = hit[8]
 #                   if hit[6] != ['NA'] and not hit[6] in counted[i][-3]:
                    if  not hit[6] in counted[i][-3]:                    
                        counted[i][-3] += hit[6]
                        counted[i][-3] = list(set(counted[i][-3]))
                    print("DEBUG BORDEL COUNTED APRES "+str(i)+" ",counted[i])    

        else:
            print("DEBUG 9 UNSELECTED ",hit)
            start = hit[1]
            end = hit[2]
            chr = hit[0]
            i += 1
            counted.append(["locus"+str(i)+'_'+args.type, chr, start, end, "-",
               0,0,0,0,0,0,0,0,0,0])
            counted[i][int(hit[4])+3] += 1
            if args.type == 'prot':
                print("DEBUG 12 SELECTED ",hit)
                counted[i].append(hit[5])
                #print("debug  test exon 4 ",hit[6])
                #if hit[6] != ['NA']:
                #    counted[i].append(hit[6])
                counted[i].append(hit[6])   
                counted[i].append(hit[7])
                counted[i].append(hit[8])
    for deb in counted:
        print("DEBUG RETURN ",deb)            
    return counted

####################################
### Proceed to execution

i = 2
hitlist = []

if args.type == 'nucl':
    entry_type = 'tblastn'
elif args.type == 'prot':
    entry_type = 'blastp'
else:
    raise SyntaxError
if args.type == 'prot' and not args.gff:
    raise SyntaxError
# Perform from exon 2 to 11.
while i < 12:
    hitlist.append(exonRead(args.input+"/PRDM9_CDS_exon"+str(i)+"."+entry_type+".fmt7", i))
    
    i += 1

hits = []
# Concatenate all hits in a same list of lists.
for hit in hitlist:
    hits += hit

if args.type == 'prot' and hits != []:
    print(f"Opening gff file {args.gff}")
    hits = getPosGff(hits, args.gff)

# Sort, separate, count exons.
#print("Number of hits : "+str(len(hits)))
hits = strandSeparate(hits)
for hit in hits[0]:
    print(hit)
for hit in hits[1]:
    print(hit)
#print("Number of hits : "+str(len(hits)))

if len(hits[0]) > 0 and len(hits[1]) > 0 :
    print("Avant "+ str(len(hits)))
    hits = exonCount(hits)
    print("Apres "+ str(len(hits)))
for hit in hits:
    print(hit)    
# Write results in output file.
with open(args.output, "w") as final:
    if args.type == 'nucl':
        final.write("LocusID\tchr\tstart\tend\tstrand\texon2\texon3\texon4\texon5\texon6\texon7\texon8\texon9\texon10\texon11\n")
        for hit in hits:
    # Will first skip hit if no SET hit.
    # Will then skip if no other hit than SET.
            if hit[11] == 0 and hit[12] == 0 and hit[13] == 0:
                continue
            if hit[5] == 0 and hit[6] == 0 and hit[7] == 0 and hit[8] == 0 and hit[9] == 0 and hit[10] == 0 and hit[14] == 0:
                continue
            for elem in hit:
                final.write(str(elem)+"\t")
            final.write("\n")
    elif args.type == 'prot':
        final.write("LocusID\tchr\tstart\tend\tstrand\tprot_id\tintrons\tref\tlocus_start\n")
        for hit in hits:
            if len(hit) > 0 :
                print(hit)
                #print("debug len" + str( len(hit)))
                #if len(hit) == 19:
                no_exon = [hit[0], hit[1], hit[2], hit[3], hit[4], hit[-4], hit[-3], hit[-2], hit[-1]]
                #elif len(hit) == 18:
                #print("debug "+ str(hit[-4]))
                for elem in no_exon:
                    final.write(str(elem)+"\t")
            final.write("\n")
