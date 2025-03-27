import pandas as pd

#accession_number = snakemake.params.accession

#organisms_file = snakemake.input.organisms_file

#domain_per_sequence_tabulated_file = snakemake.input.domain_per_sequence_tabulated

#domain_per_domain_summary_file = snakemake.input.domain_per_domain_summary
domain_per_domain_summary_file = "KRAB_domains_summary"
#output_file = snakemake.output[0]
output_file = "toto.csv"
#
# def process_domain_tabulated(domain, domain_tabulated_file, accession_number=accession_number):
#     '''
#     Lit les fichiers résultat de hmm_search après mise en forme (1 fichier pour chaque domain protéique) et saisit les valeurs d'intérêt (E-value, Score) dans un data frame
#     '''
#     with open(domain_tabulated_file) as reader:
#         for line in reader.readlines():
#             line_data = line.split('\t')
#             seq_id = line_data[0]
#             if seq_id in summarised_data['SeqID'].values:
#                 summarised_data.loc[summarised_data['SeqID'] == seq_id, f"{domain} Query"] = line_data[3]
#                 summarised_data.loc[summarised_data['SeqID'] == seq_id, f"{domain} E-value"] = line_data[6]
#                 summarised_data.loc[summarised_data['SeqID'] == seq_id, f"{domain} Score"] = line_data[7]

def process_domain(domain, domain_summary_file, accession_number="toto"):
    '''
    Récupère les informations importantes (nombre de domaines identifiés, position) dans les fichiers résultat de hmm_search --domtblout et les saisit dans un dataframe
    '''
    dico = {}
    with open(domain_summary_file) as reader:
        for line in reader.readlines()[1:]:
            line_data = line.split('\t')
            # nb_domains = line_data[9]
            seq_id = line_data[0]
            hmm_id = line_data[3]
            prot_len = int(line_data[2])
            hmm_len = int(line_data[5])
            hmm_range = [int(line_data[15]),int(line_data[16])]
            prot_range = [int(line_data[19]),int(line_data[20])]
            print(hmm_id,hmm_len,hmm_range,hmm_range[1]-hmm_range[0]+1)
            print(seq_id,prot_len,prot_range,prot_range[1]-prot_range[0]+1)
            if seq_id in dico :
                _val = dico[seq_id]
                hmm = _val[0]
                sequence = _val[1]
                for i in range(hmm_range[0],hmm_range[1]) :
                    hmm[i-1] += 1
                for i in range(prot_range[0],prot_range[1]) :
                    sequence[i-1] += 1
                #print(hmm)
                dico[seq_id] = [hmm,sequence]

            else :
                hmm = [0] * hmm_len
                for i in range(hmm_range[0],hmm_range[1]) :
                    hmm[i-1] = 1
                #print(hmm)
                sequence  = [0] * prot_len
                for i in range(prot_range[0],prot_range[1]) :
                    sequence[i-1] = 1
                #print(sequence)
                dico[seq_id] = [hmm,sequence]

    for seq_id in  dico:
        print(seq_id)
        _val = dico[seq_id]
        print(_val)
        hmm = _val[0]
        couv_hmm = 0
        for i in hmm:
            if i > 0 :
                couv_hmm += 1
        score_hmm = couv_hmm / len(hmm)
        print("Couverture hmm "+str(couv_hmm))
        print("Score hmm "+str(score_hmm))

        prot = _val[0]
        couv_prot = 0
        for i in prot:
            if i > 0 :
                couv_prot += 1
        score_prot = couv_prot / len(hmm)
        print("Couverture prot "+str(couv_prot))
        print("Score prot "+str(score_prot))

            # if seq_id in summarised_data['SeqID'].values:
            #     summarised_data.loc[summarised_data['SeqID'] == seq_id, f"Nb {domain} domains"] = nb_domains
            #     summarised_data.loc[summarised_data['SeqID'] == seq_id, f"{domain} domain start"] = line_data[17]
            #     summarised_data.loc[summarised_data['SeqID'] == seq_id, f"{domain} domain end"]= line_data[18]
#
# def getTaxid(accession_number=accession_number,input_file=organisms_file):
#     df = pd.read_csv(organisms_file, sep='\t', header=0)
#     taxid = df.loc[df['Assembly Accession'] == accession_number, 'Taxid'].values[0]
#     summarised_data["Taxid"] = taxid
#
# def getSpecies(accession_number=accession_number,input_file=organisms_file):
#     df = pd.read_csv(organisms_file, sep='\t', header=0)
#     species = df.loc[df['Assembly Accession'] == accession_number, 'Species Name'].values[0]
#     summarised_data["Species"] = species

#domain = domain_per_sequence_tabulated_file.split("/")[-1].split("_")[0]

# noms_colonnes = ['SeqID']
# noms_colonnes.append('Assembly')
# noms_colonnes.append(domain+' Query')
# noms_colonnes.append(domain+' E-value')
# noms_colonnes.append(domain+' Score')
# noms_colonnes.append('Nb '+domain+' domains')
# noms_colonnes.append(domain+' domain start')
# noms_colonnes.append(domain+' domain end')


data_list = []
noms_colonnes = ['SeqID']

# print(".... Accession  " + accession_number )
# print(".... Domain  " + domain )
# # All candidates must have the domain
# print(".... Processing file "+domain_per_sequence_tabulated_file)
#
# with open(domain_per_sequence_tabulated_file) as reader:
#     for line in reader:
#         line_data = line.strip().split('\t')
#         to_add = {'SeqID': line_data[0], 'Assembly':accession_number, domain+' Query': line_data[2], domain+' E-value': line_data[7], domain+' Score': line_data[8]}
#         data_list.append(to_add)
#     summarised_data = pd.DataFrame(data_list, columns=noms_colonnes)

summarised_data = pd.DataFrame(data_list, columns=noms_colonnes)


#process_domain_tabulated(domain,domain_per_domain_summary_file)
process_domain("KRAB",domain_per_domain_summary_file)


#getTaxid()
#getSpecies()
summarised_data = summarised_data.fillna(0)
print("Output file = "+output_file)
summarised_data.to_csv(output_file, sep=';')
