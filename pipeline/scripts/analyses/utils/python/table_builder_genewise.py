import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Reads hmm_search processed files and domains processed files and creates an overview table in the csv format')

parser.add_argument('-a', '--accession', type=str, required=True, help='Organism genome NCBI accession number')
parser.add_argument('-e', '--error', type=str, required=True, help='Error check file')
parser.add_argument('-o', '--output', type=str, required=True, help='Processed file path')

args = parser.parse_args()

accession_number = args.accession
noms_colonnes = ['SeqID', 'Chromosome', 'Chr Start', 'Chr End', 'Strand', 'Protein Length', 'ProtRefID', 'Pseudogene (Genewise)', 
                'Stop/Shift Positions', 'SET Query', 'SET E-value', 'SET Score', 'Nb SET domains', 'SET non-truncated', 'SET domain start',
                'SET domain end', 'SET Intron', 'SET Stop/Frameshift', 'KRAB Query', 'KRAB E-value', 'KRAB Score', 'Nb KRAB domains',
                'KRAB non-truncated', 'KRAB domain start', 'KRAB domain end', 'KRAB Intron', 'KRAB Stop/Frameshift',
                'SSXRD Query', 'SSXRD E-value', 'SSXRD Score', 'Nb SSXRD domains', 'SSXRD non-truncated', 'SSXRD domain start',
                'SSXRD domain end', 'SSXRD Intron', 'SSXRD Stop/Frameshift', 'ZF Query', 'ZF E-value', 'ZF Score',
                'Nb ZF domains', 'ZF non-truncated', 'ZF domain start', 'ZF domain end', 'ZF Intron', 'ZF Stop/Frameshift']


data_list = []

# All candidates must have a SET domain
with open(f"results/{accession_number}/Step4_Hmm/tbl/SET_processed") as reader:
    for line in reader:
        line_data = line.strip().split('\t')
        to_add = {'SeqID': line_data[0], 'SET Query': line_data[2], 'SET E-value': line_data[7], 'SET Score': line_data[8]}
        data_list.append(to_add)
if args.error == 'results/'+ args.accession + '/error_check/errcheck.error_detected':
    to_add = {'SeqID': f'{accession_number}-ERROR', 'SET Query': 0, 'SET E-value': 0, 'SET Score': 0}
    data_list.append(to_add)

summarised_data = pd.DataFrame(data_list, columns=noms_colonnes)


def processed_data(domain, accession_number=accession_number):
    '''
    Lit les fichiers résultat de hmm_search après mise en forme (1 fichier pour chaque domain protéique) et saisit les valeurs d'intérêt (E-value, Score) dans un data frame
    '''
    with open(f"results/{accession_number}/Step4_Hmm/domtbl/{domain}_domains_processed") as reader:
        for line in reader.readlines():
            line_data = line.split('\t')
            seq_id = line_data[0]
            if seq_id in summarised_data['SeqID'].values:
                summarised_data.loc[summarised_data['SeqID'] == seq_id, f"{domain} Query"] = line_data[3]
                summarised_data.loc[summarised_data['SeqID'] == seq_id, f"{domain} E-value"] = line_data[6]
                summarised_data.loc[summarised_data['SeqID'] == seq_id, f"{domain} Score"] = line_data[7]
    if args.error == 'results/'+ args.accession + '/error_check/errcheck.error_detected':
        summarised_data[f"{domain} Query"] = 0
        summarised_data[f"{domain} E-value"] = 0
        summarised_data[f"{domain} Score"] = 0


def processed_domains(domain, accession_number=accession_number):
    '''
    Récupère les informations importantes (nombre de domaines identifiés, position) dans les fichiers résultat de hmm_search --domtblout et les saisit dans un dataframe
    '''
    with open(f"results/{accession_number}/Step4_Hmm/domtbl/{domain}_domains_summary") as reader:
        for line in reader.readlines()[1:]:
            line_data = line.split('\t')
            nb_domains = line_data[10]
            seq_id = line_data[0]
            pseudo = line_data[22].split(',')
            print("debug")
            print(pseudo)
            shift = 0
            intron = 0
            if seq_id in summarised_data['SeqID'].values:
                summarised_data.loc[summarised_data['SeqID'] == seq_id, f"Nb {domain} domains"] = nb_domains
                summarised_data.loc[summarised_data['SeqID'] == seq_id, f"Chromosome"]= pseudo[0]
                summarised_data.loc[summarised_data['SeqID'] == seq_id, f"Chr Start"]= pseudo[4].split('-')[0]
                summarised_data.loc[summarised_data['SeqID'] == seq_id, f"Chr End"]= pseudo[4].split('-')[1]
                summarised_data.loc[summarised_data['SeqID'] == seq_id, f"Strand"]= pseudo[5]
                summarised_data.loc[summarised_data['SeqID'] == seq_id, f"Protein Length"]= pseudo[3]
                summarised_data.loc[summarised_data['SeqID'] == seq_id, f"ProtRefID"]= pseudo[6]
                summarised_data.loc[summarised_data['SeqID'] == seq_id, f"Pseudogene (Genewise)"] = pseudo[7]
                summarised_data.loc[summarised_data['SeqID'] == seq_id, f"Stop/Shift Positions"] = pseudo[1]

## ZF has to be treated differently because all domains must be accounted for
                if domain == 'ZF' and not summarised_data.loc[summarised_data['SeqID'] == seq_id, f"{domain} domain start"].isnull().iloc[0]:
                    summarised_data.loc[summarised_data['SeqID'] == seq_id, f"{domain} domain end"]= line_data[18]
                
                else:
                    summarised_data.loc[summarised_data['SeqID'] == seq_id, f"{domain} domain start"]= line_data[17]
                    summarised_data.loc[summarised_data['SeqID'] == seq_id, f"{domain} domain end"]= line_data[18]

                posint = pseudo[2].split(';')
                if posint == ['NA']:
                    posint = []
                posshift = pseudo[1].split(';')
                if posshift == ['NA']:
                    posshift = []
                    
                #print(posint)
                if posint != [''] :
                    for i in range(len(posint)):
                        #print(i,posint[i])
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
    
    if args.error == 'results/'+ args.accession + '/error_check/errcheck.error_detected':
        summarised_data[f"Nb {domain} domains"] = 0
        summarised_data[f"Chromosome"] = 0
        summarised_data[f"Chr Start"] = 0
        summarised_data[f"Chr End"] = 0
        summarised_data[f"Strand"] = 0
        summarised_data[f"Protein Length"] = 0
        summarised_data[f"ProtRefID"] = 0
        summarised_data[f"Pseudogene (Genewise)"] = 0
        summarised_data[f"Stop/Shift Positions"] = 0
        summarised_data[f"{domain} domain start"] = 0
        summarised_data[f"{domain} domain end"] = 0
        summarised_data[f"{domain} Intron"] = 0
        summarised_data[f"{domain} Stop/Frameshift"] = 0
        summarised_data[f"{domain} non-truncated"] = 0

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

if args.error == 'results/'+ args.accession + '/error_check/errcheck.error_detected':
        summarised_data[f"SeqID"] = f'{accession_number}-ERROR'

summarised_data = summarised_data.fillna(0)

summarised_data.to_csv(args.output, sep=';')
