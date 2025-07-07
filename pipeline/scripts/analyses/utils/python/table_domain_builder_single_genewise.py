import pandas as pd

accession_number = snakemake.params.accession

organisms_file = snakemake.input.organisms_file

domain_per_sequence_tabulated_file = snakemake.input.domain_per_sequence_tabulated

domain_per_domain_summary_file = snakemake.input.domain_per_domain_summary

output_file = snakemake.output[0]

def process_domain_tabulated(domain, domain_tabulated_file, accession_number=accession_number):
    '''
    Lit les fichiers résultat de hmm_search après mise en forme (1 fichier pour chaque domain protéique) et saisit les valeurs d'intérêt (E-value, Score) dans un data frame
    '''
    with open(domain_tabulated_file) as reader:    
        for line in reader.readlines():
            line_data = line.split('\t')
            seq_id = line_data[0]
            if seq_id in summarised_data['SeqID'].values:
                summarised_data.loc[summarised_data['SeqID'] == seq_id, f"{domain} Query"] = line_data[3]
                summarised_data.loc[summarised_data['SeqID'] == seq_id, f"{domain} E-value"] = line_data[6]
                summarised_data.loc[summarised_data['SeqID'] == seq_id, f"{domain} Score"] = line_data[7]

def process_domain_summary(domain, domain_summary_file, accession_number=accession_number):
    '''
    Récupère les informations importantes (nombre de domaines identifiés, position) dans les fichiers résultat de hmm_search --domtblout et les saisit dans un dataframe
    '''
    with open(domain_summary_file) as reader:
        nb_lines = 0 
        for line in reader.readlines()[1:]:
            nb_lines += 1
            line_data = line.split('\t')
            nb_domains = line_data[9]
            seq_id = line_data[0]
            pseudo = line_data[22].split(',')        
            shift = 0
            intron = 0
            if seq_id in summarised_data['SeqID'].values:
                summarised_data.loc[summarised_data['SeqID'] == seq_id, f"Nb {domain} domains"] = nb_domains
                summarised_data.loc[summarised_data['SeqID'] == seq_id, f"{domain} domain start"] = line_data[17]
                summarised_data.loc[summarised_data['SeqID'] == seq_id, f"{domain} domain end"]= line_data[18]
                summarised_data.loc[summarised_data['SeqID'] == seq_id, f"Chromosome"]= pseudo[0]
                summarised_data.loc[summarised_data['SeqID'] == seq_id, f"Chr Start"]= pseudo[4].split('-')[0]
                summarised_data.loc[summarised_data['SeqID'] == seq_id, f"Chr End"]= pseudo[4].split('-')[1]
                summarised_data.loc[summarised_data['SeqID'] == seq_id, f"Strand"]= pseudo[5]
                summarised_data.loc[summarised_data['SeqID'] == seq_id, f"Protein Length"]= pseudo[3]
                summarised_data.loc[summarised_data['SeqID'] == seq_id, f"Genewise index"]= pseudo[8]
                summarised_data.loc[summarised_data['SeqID'] == seq_id, f"ProtRefID"]= pseudo[6]
                summarised_data.loc[summarised_data['SeqID'] == seq_id, f"Pseudogene (Genewise)"] = pseudo[7]
                summarised_data.loc[summarised_data['SeqID'] == seq_id, f"Stop/Shift Positions"] = pseudo[1]
      
                posint = pseudo[2].split(';')
                if posint == ['NA']:
                    posint = []
                posshift = pseudo[1].split(';')
                if posshift == ['NA']:
                    posshift = []
                    
                #print(posint)
                if posint != [''] :
                    for i in range(len(posint)):
                        if int(summarised_data.loc[summarised_data['SeqID'] == seq_id, f"{domain} domain start"]) <= int(posint[i]) <= int(summarised_data.loc[summarised_data['SeqID'] == seq_id, f"{domain} domain end"]):
                            intron += 1
                    for i in range(len(posshift)):
                        if int(summarised_data.loc[summarised_data['SeqID'] == seq_id, f"{domain} domain start"]) <= int(posshift[i]) <= int(summarised_data.loc[summarised_data['SeqID'] == seq_id, f"{domain} domain end"]):
                            shift += 1

                summarised_data.loc[summarised_data['SeqID'] == seq_id, f"{domain} Intron"] = intron
                summarised_data.loc[summarised_data['SeqID'] == seq_id, f"{domain} Stop/Frameshift"] = shift

## Counting non-truncated domains
                if posshift != []:
                    if int(summarised_data.loc[summarised_data['SeqID'] == seq_id, f"{domain} domain end"].iloc[0]) <= int(posshift[0]):
                        if not summarised_data.loc[summarised_data['SeqID'] == seq_id, f"{domain} non-truncated"].isnull().iloc[0]:
                            summarised_data.loc[summarised_data['SeqID'] == seq_id, f"{domain} non-truncated"] += 1
                        else:
                            summarised_data.loc[summarised_data['SeqID'] == seq_id, f"{domain} non-truncated"] = 1
                else:
                    if not summarised_data.loc[summarised_data['SeqID'] == seq_id, f"{domain} non-truncated"].isnull().iloc[0]:
                        summarised_data.loc[summarised_data['SeqID'] == seq_id, f"{domain} non-truncated"] += 1
                    else:
                        summarised_data.loc[summarised_data['SeqID'] == seq_id, f"{domain} non-truncated"] = 1      
        if nb_lines == 0 :
            summarised_data.insert(loc=len(summarised_data.columns)-2, column=f"{domain} Stop/Frameshift",value="NA")
            summarised_data.insert(loc=len(summarised_data.columns)-2, column=f"{domain} Intron",value="NA")
            #noms_colonnes.append(f"{domain} Stop/Frameshift Positions")

def process_hmm_cov(domain, domain_summary_file, accession_number=accession_number):
    '''
    '''
    dico = {}
    with open(domain_summary_file) as reader:
        for line in reader.readlines()[1:]:
            line_data = line.split('\t')
            seq_id = line_data[0]
            hmm_id = line_data[3]
            prot_len = int(line_data[2])
            hmm_len = int(line_data[5])
            hmm_range = [int(line_data[15]),int(line_data[16])]
            prot_range = [int(line_data[19]),int(line_data[20])]
            if seq_id in dico :
                _val = dico[seq_id]
                hmm = _val[0]
                sequence = _val[1]
                for i in range(hmm_range[0],hmm_range[1]) :
                    hmm[i-1] += 1
                    if hmm[i-1] > 1 :
                        hmm[i-1] = 1
                for i in range(prot_range[0],prot_range[1]) :
                    sequence[i-1] += 1
                    if sequence[i-1] > 1 :
                        sequence[i-1] = 1
                dico[seq_id] = [hmm,sequence]

            else :
                hmm = [0] * hmm_len
                for i in range(hmm_range[0],hmm_range[1]) :
                    hmm[i-1] = 1
                sequence  = [0] * prot_len
                for i in range(prot_range[0],prot_range[1]) :
                    sequence[i-1] = 1
                dico[seq_id] = [hmm,sequence]

    for seq_id in  dico:
        _val = dico[seq_id]
        hmm = _val[0]
        couv_hmm = 0
        ii = 0 # (index)
        limit = [] # [debut, fin] de la couverture
        limits = []
        curr_val = hmm[ii]
        
        if curr_val == 1 :
            limit.append( ii + 1 ) #ajoute 1 car la 1ere position est 1
        for i in hmm:
            if i != curr_val:
                if i == 1 :
                    limit = []
                    limit.append( ii + 1 )
                    curr_val = i
                if i == 0 :
                    limit.append( ii  )
                    limits.append(limit)
                    curr_val = i   # on ne rajoute pas 1 ici car i+1 est un 0
            if i > 0 :
                couv_hmm += 1
            ii += 1
        score_hmm = int (100 * couv_hmm / len(hmm))/100

        prot = _val[1]
        couv_prot = 0
        for i in prot:
            if i > 0 :
                couv_prot += 1
        score_prot = int (100 * couv_prot / len(hmm))/100        

        if seq_id in summarised_data['SeqID'].values:
            summarised_data.loc[summarised_data['SeqID'] == seq_id, f"{domain} HMM cov."] = score_hmm 
            summarised_data.loc[summarised_data['SeqID'] == seq_id, f"{domain} HMM cov. pos."] = str(limits) 
            summarised_data.loc[summarised_data['SeqID'] == seq_id, f"{domain} Prot cov."] = score_prot 
        
        
def getTaxid(accession_number=accession_number,input_file=organisms_file):
    df = pd.read_csv(organisms_file, sep='\t', header=0)
    taxid = df.loc[df['Assembly Accession'] == accession_number, 'Taxid'].values[0]    
    summarised_data["Taxid"] = taxid

def getSpecies(accession_number=accession_number,input_file=organisms_file):
    df = pd.read_csv(organisms_file, sep='\t', header=0)
    species = df.loc[df['Assembly Accession'] == accession_number, 'Species Name'].values[0]    
    summarised_data["Species"] = species    

domain = domain_per_sequence_tabulated_file.split("/")[-1].split("_")[0]  

noms_colonnes = ['SeqID']
noms_colonnes.append('Assembly')
noms_colonnes.append(domain+' Query')
noms_colonnes.append(domain+' E-value')
noms_colonnes.append(domain+' Score')
noms_colonnes.append('Nb '+domain+' domains')
noms_colonnes.append(domain+' domain start')
noms_colonnes.append(domain+' domain end')

noms_colonnes.append(domain+' HMM cov.')
noms_colonnes.append(domain+' HMM cov. pos.')
noms_colonnes.append(domain+' Prot cov.')

noms_colonnes.append(domain+' non-truncated')

noms_colonnes.append('Chromosome')
noms_colonnes.append('Chr Start')
noms_colonnes.append('Chr End')
noms_colonnes.append('Strand')
noms_colonnes.append('Protein Length')
noms_colonnes.append('Genewise index')
noms_colonnes.append('ProtRefID')
noms_colonnes.append('Pseudogene (Genewise)')
noms_colonnes.append('Stop/Shift Positions')

data_list = []


print(".... Accession  " + accession_number )
print(".... Domain  " + domain )
# All candidates must have the domain
print(".... Processing file "+domain_per_sequence_tabulated_file)

with open(domain_per_sequence_tabulated_file) as reader:
    for line in reader:
        line_data = line.strip().split('\t')
        to_add = {'SeqID': line_data[0], 'Assembly':accession_number, domain+' Query': line_data[2], domain+' E-value': line_data[7], domain+' Score': line_data[8]}
        data_list.append(to_add)
    summarised_data = pd.DataFrame(data_list, columns=noms_colonnes)


process_domain_tabulated(domain,domain_per_domain_summary_file)
process_domain_summary(domain,domain_per_domain_summary_file)  
process_hmm_cov(domain,domain_per_domain_summary_file)   
 
getTaxid()                
getSpecies()
# On vire les sequences qui sont dans tbl mais pas dans domtbl ce qui induit une sequence 
# avec in Chromosome = NA
print(summarised_data[summarised_data['Chromosome'].isnull()]) 
missing = summarised_data[summarised_data['Chromosome'].isnull()]
missing_seqids = missing["SeqID"].values
print("Missing seq = "+missing_seqids)
summarised_data.drop(summarised_data[summarised_data['Chromosome'].isnull()].index, inplace = True)
summarised_data = summarised_data.fillna(0)    
print("Output file = "+output_file)        
summarised_data.drop(summarised_data.columns[summarised_data.columns.str.contains('unnamed', case=False)], axis=1, inplace=True)         

summarised_data.to_csv(output_file, sep=';')
#summarised_data.to_csv(output_file, sep=';', index = False)
