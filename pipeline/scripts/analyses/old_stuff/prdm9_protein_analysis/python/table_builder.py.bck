import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Reads hmm_search processed files and domains processed files and creates an overview table in the csv format')

parser.add_argument('-a', '--accession', type=str, required=True, help='Organism genome NCBI accession number')
parser.add_argument('-o', '--output', type=str, required=True, help='Processed file path')

args = parser.parse_args()

accession_number = args.accession
noms_colonnes = ['SeqID', 'SET Query', 'SET E-value', 'SET Score', 'Nb SET domains', 'SET domain start', 'SET domain end',
                'KRAB Query', 'KRAB E-value', 'KRAB Score', 'Nb KRAB domains', 'KRAB domain start', 'KRAB domain end',
                'SSXRD Query', 'SSXRD E-value', 'SSXRD Score', 'Nb SSXRD domains', 'SSXRD domain start', 'SSXRD domain end',
                'ZF Query', 'ZF E-value', 'ZF Score', 'Nb ZF domains', 'ZF domain start', 'ZF domain end']


data_list = []

# All candidates must have a SET domain
with open(f"results/{accession_number}/hmm_search/tbl/SET_processed") as reader:
    for line in reader:
        line_data = line.strip().split('\t')
        to_add = {'SeqID': line_data[0], 'SET Query': line_data[2], 'SET E-value': line_data[7], 'SET Score': line_data[8]}
        data_list.append(to_add)

summarised_data = pd.DataFrame(data_list, columns=noms_colonnes)


def processed_data(domain, accession_number=accession_number):
    '''
    Lit les fichiers résultat de hmm_search après mise en forme (1 fichier pour chaque domain protéique) et saisit les valeurs d'intérêt (E-value, Score) dans un data frame
    '''
    with open(f"results/{accession_number}/hmm_search/domtbl/{domain}_domains_processed") as reader:
        for line in reader.readlines():
            line_data = line.split('\t')
            seq_id = line_data[0]
            if seq_id in summarised_data['SeqID'].values:
                summarised_data.loc[summarised_data['SeqID'] == seq_id, f"{domain} Query"] = line_data[3]
                summarised_data.loc[summarised_data['SeqID'] == seq_id, f"{domain} E-value"] = line_data[6]
                summarised_data.loc[summarised_data['SeqID'] == seq_id, f"{domain} Score"] = line_data[7]


def processed_domains(domain, accession_number=accession_number):
    '''
    Récupère les informations importantes (nombre de domaines identifiés, position) dans les fichiers résultat de hmm_search --domtblout et les saisit dans un dataframe
    '''
    with open(f"results/{accession_number}/hmm_search/domtbl/{domain}_domains_summary") as reader:
        for line in reader.readlines()[1:]:
            line_data = line.split('\t')
            nb_domains = line_data[9]
            seq_id = line_data[0]
            if seq_id in summarised_data['SeqID'].values:
                summarised_data.loc[summarised_data['SeqID'] == seq_id, f"Nb {domain} domains"] = nb_domains
                summarised_data.loc[summarised_data['SeqID'] == seq_id, f"{domain} domain start"] = line_data[17]
                summarised_data.loc[summarised_data['SeqID'] == seq_id, f"{domain} domain end"]= line_data[18]

def getTaxid(accession_number=accession_number):
    df = pd.read_csv('data/resources/organisms_data', sep='\t', header=0)
    taxid = df.loc[df['Assembly Accession'] == accession_number, 'Taxid'].values[0]
    summarised_data["Taxid"] = taxid


domain1 = 'KRAB'
processed_data(domain1)
processed_domains(domain1)
domain2 = 'SSXRD'
processed_data(domain2)
processed_domains(domain2)
domain3 = 'ZF'
processed_data(domain3)
processed_domains(domain3)
domain4 = 'SET'
processed_domains(domain4)
getTaxid()

summarised_data = summarised_data.fillna(0)                                   
summarised_data.to_csv(args.output, sep=';')
