import pandas as pd


#parser = argparse.ArgumentParser(description='Reads hmm_search processed files and domains processed files and creates an overview table in the csv format')

#parser.add_argument('-i', '--input_dir', type=str, required=True, help='Input dir path')
#parser.add_argument('-a', '--accession', type=str, required=True, help='Organism genome NCBI accession number')
#parser.add_argument('-o', '--output', type=str, required=True, help='Output file path')
#parser.add_argument('-s', '--set', type=str, required=True, help='Tabulated SET domain hmmsearch results file path')
#args = parser.parse_args()

accession_number = snakemake.params.accession
#input_dir = args.input_dir

SET_per_sequence_tabulated_file = snakemake.input.SET_per_sequence_tabulated

SET_per_domain_tabulated_file = snakemake.input.SET_per_domain_tabulated
KRAB_per_domain_tabulated_file = snakemake.input.KRAB_per_domain_tabulated
SSXRD_per_domain_tabulated_file = snakemake.input.SSXRD_per_domain_tabulated
ZF_per_domain_tabulated_file = snakemake.input.ZF_per_domain_tabulated


SET_per_domain_summary_file = snakemake.input.SET_per_domain_summary
KRAB_per_domain_summary_file = snakemake.input.KRAB_per_domain_summary
SSXRD_per_domain_summary_file = snakemake.input.SSXRD_per_domain_summary
ZF_per_domain_summary_file = snakemake.input.ZF_per_domain_summary

#SET_hmm_tab = args.s
noms_colonnes = ['SeqID', 'SET Query', 'SET E-value', 'SET Score', 'Nb SET domains', 'SET domain start', 'SET domain end',
                'KRAB Query', 'KRAB E-value', 'KRAB Score', 'Nb KRAB domains', 'KRAB domain start', 'KRAB domain end',
                'SSXRD Query', 'SSXRD E-value', 'SSXRD Score', 'Nb SSXRD domains', 'SSXRD domain start', 'SSXRD domain end',
                'ZF Query', 'ZF E-value', 'ZF Score', 'Nb ZF domains', 'ZF domain start', 'ZF domain end']


data_list = []

# All candidates must have a SET domain
#with open(f"{input_dir}/{accession_number}/analyses/prdm9_prot/hmm_search/tbl/SET_tabulated") as reader:
with open(SET_per_sequence_tabulated_file) as reader:
    for line in reader:
        line_data = line.strip().split('\t')
        to_add = {'SeqID': line_data[0], 'SET Query': line_data[2], 'SET E-value': line_data[7], 'SET Score': line_data[8]}
        data_list.append(to_add)

summarised_data = pd.DataFrame(data_list, columns=noms_colonnes)


def process_domain_tabulated(domain, domain_tabulated_file, accession_number=accession_number):
    '''
    Lit les fichiers résultat de hmm_search après mise en forme (1 fichier pour chaque domain protéique) et saisit les valeurs d'intérêt (E-value, Score) dans un data frame
    '''
#    with open(f"{input_dir}/{accession_number}/analyses/prdm9_prot/hmm_search/domtbl/{domain}_domains_tabulated") as reader:
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
#    with open(f"{input_dir}/{accession_number}/analyses/prdm9_prot/hmm_search/domtbl/{domain}_domains_summary") as reader:
    with open(domain_summary_file) as reader:    
        for line in reader.readlines()[1:]:
            line_data = line.split('\t')
            nb_domains = line_data[9]
            seq_id = line_data[0]
            if seq_id in summarised_data['SeqID'].values:
                summarised_data.loc[summarised_data['SeqID'] == seq_id, f"Nb {domain} domains"] = nb_domains
                summarised_data.loc[summarised_data['SeqID'] == seq_id, f"{domain} domain start"] = line_data[17]
                summarised_data.loc[summarised_data['SeqID'] == seq_id, f"{domain} domain end"]= line_data[18]

#def getTaxid(accession_number=accession_number,input_dir=input_dir):
#    df = pd.read_csv(input_dir+'/../organisms_data', sep='\t', header=0)
#    taxid = df.loc[df['Assembly Accession'] == accession_number, 'Taxid'].values[0]    
#    print("DEBUG "+accession_number+": "+str(taxid))
#    summarised_data["Taxid"] = taxid


#domain1 = 'KRAB'
process_domain_tabulated("KRAB",KRAB_per_domain_tabulated_file )
process_domain_summary("KRAB",KRAB_per_domain_summary_file )

process_domain_tabulated("SET",SET_per_domain_tabulated_file )
process_domain_summary("SET",SET_per_domain_summary_file )

process_domain_tabulated("SSXRD",SSXRD_per_domain_tabulated_file )
process_domain_summary("SSXRD",SSXRD_per_domain_summary_file )

process_domain_tabulated("ZF",ZF_per_domain_tabulated_file )
process_domain_summary("ZF",ZF_per_domain_summary_file )

#processed_domains(domain1)
#domain2 = 'SSXRD'
#processed_data(domain2)
#processed_domains(domain2)
#domain3 = 'ZF'
#processed_data(domain3)
#processed_domains(domain3)
#domain4 = 'SET'
#processed_domains(domain4)
# a remettre getTaxid()

summarised_data = summarised_data.fillna(0)                                   
summarised_data.to_csv(args.output, sep=';')
