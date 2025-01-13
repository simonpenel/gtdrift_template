import pandas as pd

accession_number = snakemake.params.accession

organisms_file = snakemake.input.organisms_file

domain_per_sequence_tabulated_files = snakemake.input.domain_per_sequence_tabulated

domain_per_domain_summary_files = snakemake.input.domain_per_domain_summary

output_files = snakemake.output

# Build an array of domains
domains = [] 
for file in domain_per_sequence_tabulated_files:
    domain = file.split("/")[-1].split("_")[0]
    if  not domain in domains:
        domains.append(domain)

# Build an array of accessions        
accessions = []

# Build a dictionary accession -> [per_sequence files]
dico_per_seq = {}
for file in domain_per_sequence_tabulated_files:
    _buflist = file.split("/")[:-5]
    accession = _buflist.pop() 
    if not   accession in dico_per_seq:
        dico_per_seq[accession] = []
    dico_per_seq[accession].append(file)
    
# Build a dictionary accession -> [per_domain files]    
dico_per_dom = {}
for file in domain_per_domain_summary_files:
    _buflist = file.split("/")[:-5]
    accession = _buflist.pop() 
    if not   accession in dico_per_dom:
        dico_per_dom[accession] = []
    dico_per_dom[accession].append(file)

# Build a dictionary accession -> [output files] (only 1 element in the array)
dico_output = {}
for file in output_files:
    _buflist = file.split("/")[:-3]
    accession = _buflist.pop() 
    if  not accession in accessions:
        accessions.append(accession) 
    if not   accession in dico_output:
        dico_output[accession] = []
    dico_output[accession].append(file)
    
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
      
def getTaxid(accession_number=accession_number,input_file=organisms_file):
    df = pd.read_csv(organisms_file, sep='\t', header=0)
    taxid = df.loc[df['Assembly Accession'] == accession_number, 'Taxid'].values[0]    
    summarised_data["Taxid"] = taxid
    

for accession in accessions:
    print("\n\nProcessing genome "+accession)
    files_seq =  dico_per_seq[accession]
    files_dom = dico_per_dom[accession]
    file_output = dico_output[accession][0]  #ce tableau ne doit avoir qu'1 element
    
    print(".... Output file =        "+ file_output)
    print(".... Per sequence files = ", files_seq) 
    print(".... Per domain files =   ", files_dom)  
    noms_colonnes = ['SeqID']
    for domain in domains:
        noms_colonnes.append(domain+' Query')
        noms_colonnes.append(domain+' E-value')
        noms_colonnes.append(domain+' Score')
        noms_colonnes.append('Nb '+domain+' domains')
        noms_colonnes.append(domain+' domain start')
        noms_colonnes.append(domain+' domain end')

    data_list = []

    # All candidates must have a lest the 1st domain
    print(".... Processing 1st domain file "+files_seq[0])
    with open(files_seq[0]) as reader:
        domain = files_seq[0].split("/")[-1].split("_")[0]
        for line in reader:
            line_data = line.strip().split('\t')
            to_add = {'SeqID': line_data[0], domain+' Query': line_data[2], domain+' E-value': line_data[7], domain+' Score': line_data[8]}
            data_list.append(to_add)
        summarised_data = pd.DataFrame(data_list, columns=noms_colonnes)

    for file in       files_dom:
        print(".... Processing file " + file)
        domain = file.split("/")[-1].split("_")[0]
        print(".... Associated domain " + domain)
        process_domain_tabulated(domain,file)
        #print("Debug 2",summarised_data)  
        process_domain_summary(domain,file) 
        #print("Debug 3",summarised_data)  
 
    getTaxid()                
    summarised_data = summarised_data.fillna(0)                     
    summarised_data.to_csv(dico_output[accession][0], sep=';')
