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
        for line in reader.readlines()[1:]:
            line_data = line.split('\t')
            nb_domains = line_data[9]
            seq_id = line_data[0]
            if seq_id in summarised_data['SeqID'].values:
                summarised_data.loc[summarised_data['SeqID'] == seq_id, f"Nb {domain} domains"] = nb_domains
                summarised_data.loc[summarised_data['SeqID'] == seq_id, f"{domain} domain start"] = line_data[17]
                summarised_data.loc[summarised_data['SeqID'] == seq_id, f"{domain} domain end"]= line_data[18]
      
      
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

#domain = domain_per_sequence_tabulated_file.split("/")[-1].split("_")[0]  
domain_split = domain_per_sequence_tabulated_file.split("/")[-1].split("_")
domain_split.pop()
domain = "_".join(domain_split)
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
summarised_data = summarised_data.fillna(0)    
print("Output file = "+output_file)                 
summarised_data.to_csv(output_file, sep=';')

