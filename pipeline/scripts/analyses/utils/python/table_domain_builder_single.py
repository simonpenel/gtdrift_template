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
      
def getTaxid(accession_number=accession_number,input_file=organisms_file):
    df = pd.read_csv(organisms_file, sep='\t', header=0)
    taxid = df.loc[df['Assembly Accession'] == accession_number, 'Taxid'].values[0]    
    summarised_data["Taxid"] = taxid
    

domain = domain_per_sequence_tabulated_file.split("/")[-1].split("_")[0]  

noms_colonnes = ['SeqID']
noms_colonnes.append(domain+' Query')
noms_colonnes.append(domain+' E-value')
noms_colonnes.append(domain+' Score')
noms_colonnes.append('Nb '+domain+' domains')
noms_colonnes.append(domain+' domain start')
noms_colonnes.append(domain+' domain end')


data_list = []


print(".... Accession  " + accession_number )
print(".... Domain  " + domain )
# All candidates must have the domain
print(".... Processing file "+domain_per_sequence_tabulated_file)

with open(domain_per_sequence_tabulated_file) as reader:
    for line in reader:
        line_data = line.strip().split('\t')
        to_add = {'SeqID': line_data[0], domain+' Query': line_data[2], domain+' E-value': line_data[7], domain+' Score': line_data[8]}
        data_list.append(to_add)
    summarised_data = pd.DataFrame(data_list, columns=noms_colonnes)


process_domain_tabulated(domain,domain_per_domain_summary_file)
process_domain_summary(domain,domain_per_domain_summary_file)  
 
 
getTaxid()                
summarised_data = summarised_data.fillna(0)    
print("Output file = "+output_file)                 
summarised_data.to_csv(output_file, sep=';')
