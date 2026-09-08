import pandas as pd

import warnings
warnings.filterwarnings("error")

accession_number = snakemake.params.accession

to_be_merged = snakemake.params.to_be_merged

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

def display_info(line) :
        seqname = line[0]
        hmm_length =line[5]
        global_evalue =line[6]
        global_score =line[7]
        num_hit = line[9]
        nb_hits =  line[10]
        hit_cvalue = line[11]
        hit_score = line[13]
        start_in_hmm = line[15]
        end_in_hmm = line[16]
        start_in_prot = line[17]
        end_in_prot = line[18]        
        print(seqname +" "+num_hit+"/"+nb_hits)
        print("\thmm_length "+hmm_length)
        print("\tglobal evalue "+global_evalue)
        print("\tglobal score "+global_score)        
        print("\thit cvalue "+hit_cvalue)
        print("\thit score "+hit_score)
        print("\tstart prot "+start_in_prot)
        print("\tend prot "+end_in_prot)                  
        print("\tstart hmm "+start_in_hmm)
        print("\tend hmm "+end_in_hmm) 

def get_hmm_info(line) :
        seqname = line[0]
        hmm_length =line[5]
        global_evalue =line[6]
        global_score =line[7]
        num_hit = line[9]
        nb_hits =  line[10]
        hit_cvalue = line[11]
        hit_score = line[13]
        start_in_hmm = line[15]
        end_in_hmm = line[16]
        start_in_prot = line[17]
        end_in_prot = line[18]
        info__hit = {
            "seqname" : seqname,
            "hmm_length" : hmm_length,  
            "global_evalue" : global_evalue,
            "global_score" : global_score,
            "num_hit" : num_hit,
            "nb_hits" : nb_hits,
            "hit_cvalue" : hit_cvalue,
            "hit_score" : hit_score,
            "start_in_hmm" : start_in_hmm,
            "end_in_hmm" : end_in_hmm,
            "start_in_prot": start_in_prot,
            "end_in_prot" : end_in_prot,
            "merged" : False          
        }
        return info__hit     

def process_domain_merge(domain, domain_summary_file, accession_number=accession_number):
    '''
    Merge les domaines
    '''
    with open(domain_summary_file) as reader:
        sequence_domains = {}    
        for line in reader.readlines()[1:]:
            line_data = line.split('\t')
            seq_id = line_data[0]

            nb_domains = int(line_data[9])
            if nb_domains == 1 :
                sequence_domains[seq_id] = []
                sequence_domains[seq_id].append(line_data)
            else :
                if not seq_id in sequence_domains :
                    sequence_domains[seq_id] = []
                sequence_domains[seq_id].append(line_data)
        for sequence in sequence_domains:
            # if len(sequence_domains[sequence]) == 1 :# pour l'affichage
            #     break 
            print("\n\nSEQUENCE "+ sequence)
            index_hmm_line = 1
            merged_domains = []
            for sequence_hmm_line in sequence_domains[sequence] :
                # print("Line "+str(index_hmm_line))
                # print(sequence_hmm_line)
                display_info(sequence_hmm_line)
                hmm_info = get_hmm_info(sequence_hmm_line)
                # print(hmm_info)
                if index_hmm_line == 1 :
                    current_hmm_info = hmm_info
                    merged_hmm_info = hmm_info
                else :
                    print("Check protein position : is curent hit start "+hmm_info["start_in_prot"] + " > previous hit end " +  current_hmm_info["end_in_prot"] + " ?")
                    if int(hmm_info["start_in_prot"]) < int(current_hmm_info["end_in_prot"]):
                        print("superposed hit in protein : new hmm domain")
                        merged_domains.append(merged_hmm_info)
                        merged_hmm_info = hmm_info
                    else :
                        print("current hit in protein is compatible with previous hit")
                        print("Check hmm position : is curent hit start " + hmm_info["start_in_hmm"] + " > previous hit end " +  current_hmm_info["end_in_hmm"] + " ?")
                        if int(hmm_info["start_in_hmm"]) > int(current_hmm_info["end_in_hmm"]):
                            print("splited hmm domain")
                            merged_hmm_info["end_in_hmm"] = hmm_info["end_in_hmm"]
                            merged_hmm_info["end_in_prot"] = hmm_info["end_in_prot"]
                            merged_hmm_info["hit_score"] = str(float(merged_hmm_info["hit_score"]) + float(hmm_info["hit_score"]))
                            merged_hmm_info["merged"] = True
                        else : 
                            print("superposed hit in hmm : new hmm domain")
                            merged_domains.append(merged_hmm_info)
                            merged_hmm_info = hmm_info

                    current_hmm_info = hmm_info
                index_hmm_line += 1
            merged_domains.append(merged_hmm_info)    
            print("Merged domains :")
            print("-------------- :")
            for buf in merged_domains:
                print(buf)
            if  to_be_merged == "merged" :
                print("Combine all domains (" + domain + ")")
                combined_domain = merged_domains.pop(0)
                for merged_domain in merged_domains:
                    print("current")
                    print(combined_domain)
                    print("new ")
                    print(merged_domain)                    
                    if int(merged_domain["start_in_prot"]) >= int(combined_domain["end_in_prot"]) - 2:
                        print("ok")
                        combined_domain["end_in_prot"] = merged_domain["end_in_prot"]
                    else :
                        sys.exit("problem")
                print("Combined domain :")
                print(combined_domain)
                if sequence in summarised_data['SeqID'].values:
                    summarised_data.loc[summarised_data['SeqID'] == sequence, f"Nb {domain} hits"] = len(sequence_domains[sequence])
                    summarised_data.loc[summarised_data['SeqID'] == sequence, f"Nb {domain} domains (after merging domains)"] = 1
                    summarised_data.loc[summarised_data['SeqID'] == sequence, f"{domain} domain start"] = int(combined_domain["start_in_prot"])
                    summarised_data.loc[summarised_data['SeqID'] == sequence, f"{domain} domain end"]= int(combined_domain["end_in_prot"])

            else :
                print("Select best domain (" + domain + ")")
                # Select the domain with the highest score
                best_domain = merged_domains[0]
                score_max = best_domain["hit_score"]
                if len(merged_domains) > 1 :
                    for merged_domain in merged_domains:
                        if merged_domain["hit_score"] > score_max :
                            best_domain = merged_domain
                            score_max = merged_domain["hit_score"]
                print("Best domain :")
                print(best_domain)
                if sequence in summarised_data['SeqID'].values:
                    print("debug add "+sequence)
                    summarised_data.loc[summarised_data['SeqID'] == sequence, f"Nb {domain} hits"] = len(sequence_domains[sequence])
                    summarised_data.loc[summarised_data['SeqID'] == sequence, f"Nb {domain} domains (after merging splited hits)"] = len(merged_domains)
                    summarised_data.loc[summarised_data['SeqID'] == sequence, f"{domain} domain start"] = int(best_domain["start_in_prot"])
                    summarised_data.loc[summarised_data['SeqID'] == sequence, f"{domain} domain end"]= int(best_domain["end_in_prot"])
                    print(summarised_data["Nb "+domain+" hits"])


def process_domain_summary(domain, domain_summary_file, accession_number=accession_number):
    '''
    Récupère les informations importantes (nombre de domaines identifiés, position) dans les fichiers résultat de hmm_search --domtblout et les saisit dans un dataframe
    '''
    with open(domain_summary_file) as reader:    
        for line in reader.readlines()[1:]:
            line_data = line.split('\t')
            nb_domains = int(line_data[9])
            #if nb_domains > 1 :
                # print("warning")
                # sys.exit(1)
            seq_id = line_data[0]
            if seq_id in summarised_data['SeqID'].values:
                summarised_data.loc[summarised_data['SeqID'] == seq_id, f"Nb {domain} domains"] = nb_domains
                summarised_data.loc[summarised_data['SeqID'] == seq_id, f"{domain} domain start"] = int(line_data[17])
                summarised_data.loc[summarised_data['SeqID'] == seq_id, f"{domain} domain end"]= int(line_data[18])
      
      
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
noms_colonnes.append('Nb '+domain+' hits')
#noms_colonnes.append('Nb '+domain+' domains')
if  to_be_merged == "merged" :
    noms_colonnes.append('Nb '+domain+' domains (after merging domains)')
else :
    noms_colonnes.append('Nb '+domain+' domains (after merging splited hits)')
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
    summarised_data = summarised_data.astype({domain+' HMM cov. pos.': "string"})
    summarised_data = summarised_data.astype({domain+' domain start': "Int32"})
    summarised_data = summarised_data.astype({domain+' domain end': "Int32"})
    if  to_be_merged == "merged" :
        summarised_data = summarised_data.astype({ 'Nb '+domain+' domains (after merging domains)': "Int32"})
    else :
        summarised_data = summarised_data.astype({ 'Nb '+domain+' domains (after merging splited hits)': "Int32"})
    summarised_data = summarised_data.astype({'Nb '+ domain+' hits': "Int32"})    

print(".... Processing file "+domain_per_domain_summary_file)
process_domain_tabulated(domain,domain_per_domain_summary_file)
process_domain_merge(domain,domain_per_domain_summary_file)  
# process_domain_summary(domain,domain_per_domain_summary_file)  
process_hmm_cov(domain,domain_per_domain_summary_file)   
 
getTaxid()                
getSpecies()                

# summarised_data = summarised_data.fillna(0)    
summarised_data[domain+' HMM cov. pos.'] = summarised_data[domain+' HMM cov. pos.'].fillna(value="0")
summarised_data[domain+' domain start'] = summarised_data[domain+' domain start'].fillna(value=0)
summarised_data[domain+' domain end'] = summarised_data[domain+' domain end'].fillna(value=0)

#summarised_data['Nb '+domain+' domains'] = summarised_data['Nb '+domain+' domains'].fillna(value=0)

print("Output file = "+output_file)                 
summarised_data.to_csv(output_file, sep=';')

