import pandas as pd

import warnings
warnings.filterwarnings("error")

accession_number = snakemake.params.accession

to_be_merged = snakemake.params.to_be_merged

organisms_file = snakemake.input.organisms_file

domain_per_sequence_tabulated_file = snakemake.input.domain_per_sequence_tabulated

domain_per_domain_summary_file = snakemake.input.domain_per_domain_summary

output_file = snakemake.output[0]


def get_coverage(segments):
    max = 0
    for segment in segments :
        if segment[1] > max :
            max = segment[1]
    sequence = [0] * (max + 1)
    for segment in segments:
        index =  segment[0]
        while index <= segment[1] :
            sequence[index] = 1
            index += 1
    print(sequence)
    couverture = 0
    index = 0
    while index <= max :
        if sequence[index] == 1 :
            couverture +=1
        index +=1
    return couverture


def get_positions(segments):
    max = 0
    for segment in segments :
        if segment[1] > max :
            max = segment[1]
    sequence = [0] * (max + 1)
    for segment in segments:
        index =  segment[0]
        while index <= segment[1] :
            sequence[index] = 1
            index += 1
    print(sequence)
    positions = []
    flag = 0
    index  = 0
    position = []
    flag = sequence[index]
    if sequence[index] == 1 :
        position.add(index)
    index = 1
    while index <= max :
        if sequence[index] != flag :
            if flag == 1 :
                position.append(index)
                positions.append(position)
                position = []
                flag = sequence[index]
            else :
                position.append(index)
                flag = sequence[index]
        index +=1
    if flag == 1 :
        index -= 1
        position.append(index)
        positions.append(position)
    print("POSITIONS:")
    print(positions)
    return positions

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

def display_info(line) :
        seqname = line[0]
        hmm_length = int(line[5])
        global_evalue =float(line[6])
        global_score = float(line[7])
        num_hit = int(line[9])
        nb_hits =  int(line[10])
        hit_cvalue = float(line[11])
        hit_score = float(line[13])
        start_in_hmm = int(line[15])
        end_in_hmm = int(line[16])
        start_in_prot = int(line[17])
        end_in_prot = int(line[18])        
        print(seqname + " " + str(num_hit) + "/" + str(nb_hits))
        print("\thmm_length " + str(hmm_length))
        print("\tglobal evalue " + str(global_evalue))
        print("\tglobal score " + str(global_score))        
        print("\thit cvalue " + str(hit_cvalue))
        print("\thit score " + str(hit_score))
        print("\tstart prot " + str(start_in_prot))
        print("\tend prot " + str(end_in_prot))                  
        print("\tstart hmm " + str(start_in_hmm))
        print("\tend hmm " + str(end_in_hmm)) 

def get_hmm_info(line) :
        seqname = line[0]
        hmm_length = int(line[5])
        global_evalue = float(line[6])
        global_score = float(line[7])
        num_hit = int(line[9])
        nb_hits =  int(line[10])
        hit_cvalue = float(line[11])
        hit_score = float(line[13])
        start_in_hmm = int(line[15])
        end_in_hmm = int(line[16])
        start_in_prot = int(line[17])
        end_in_prot = int(line[18])
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
            "merged" : False,
            "segments_prot" : [],
            "segments_hmm" : [],         
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
            print("\n\nSEQUENCE "+ sequence)
            index_hmm_line = 1
            merged_domains = []
            for sequence_hmm_line in sequence_domains[sequence] :
                display_info(sequence_hmm_line)
                hmm_info = get_hmm_info(sequence_hmm_line)
                if index_hmm_line == 1 :
                    hmm_info["segments_prot"].append([hmm_info["start_in_prot"], hmm_info["end_in_prot"]])
                    hmm_info["segments_hmm"].append([hmm_info["start_in_hmm"], hmm_info["end_in_hmm"]])
                    current_hmm_info = hmm_info
                    merged_hmm_info = hmm_info
                    evalue_min = float(hmm_info["hit_cvalue"])
                else :
                    print("Check protein position : is curent hit start " + str(hmm_info["start_in_prot"]) + " > previous hit end " +  str(current_hmm_info["end_in_prot"]) + " ?")
                    if hmm_info["start_in_prot"] < current_hmm_info["end_in_prot"]:
                        print("superposed hit in protein : new hmm domain")
                        merged_domains.append(merged_hmm_info)
                        merged_hmm_info = hmm_info
                        evalue_min = float(hmm_info["hit_cvalue"])
                    else :
                        print("current hit in protein is compatible with previous hit")
                        print("Check hmm position : is curent hit start " + str(hmm_info["start_in_hmm"]) + " > previous hit end " +  str(current_hmm_info["end_in_hmm"]) + " ?")
                        if hmm_info["start_in_hmm"] > current_hmm_info["end_in_hmm"]:
                            print("splited hmm domain")
                            merged_hmm_info["end_in_hmm"] = hmm_info["end_in_hmm"]
                            merged_hmm_info["end_in_prot"] = hmm_info["end_in_prot"]
                            merged_hmm_info["hit_score"] = round((float(merged_hmm_info["hit_score"]) + float(hmm_info["hit_score"])),2)
                            merged_hmm_info["merged"] = True
                            merged_hmm_info["segments_prot"].append([hmm_info["start_in_prot"],hmm_info["end_in_prot"]])
                            merged_hmm_info["segments_hmm"].append([hmm_info["start_in_hmm"],hmm_info["end_in_hmm"]])
                            if hmm_info["hit_cvalue"] < evalue_min:
                                evalue_min = hmm_info["hit_cvalue"]
                                merged_hmm_info["hit_cvalue"] = evalue_min
                        else : 
                            print("superposed hit in hmm : new hmm domain")
                            hmm_info["segments_prot"].append([hmm_info["start_in_prot"],hmm_info["end_in_prot"]])
                            hmm_info["segments_hmm"].append([hmm_info["start_in_hmm"],hmm_info["end_in_hmm"]])
                            merged_domains.append(merged_hmm_info)
                            merged_hmm_info = hmm_info
                            evalue_min = hmm_info["hit_cvalue"]

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
                segments_hmm = []
                segments_hmm += combined_domain["segments_hmm"]
                segments_prot = []
                segments_prot += combined_domain["segments_prot"]
                for merged_domain in merged_domains:                
                    combined_domain["end_in_prot"] = merged_domain["end_in_prot"]
                    segments_hmm += merged_domain["segments_hmm"]
                    segments_prot += merged_domain["segments_prot"]

                combined_domain["segments_hmm"] = segments_hmm
                combined_domain["segments_prot"] = segments_prot
                print(combined_domain)
                # Calcul longueur cumullee en proteine    
                cumulated_length = get_coverage(combined_domain["segments_prot"])
                print("Cumulated length : "+str(cumulated_length))
                # Calcul couverture hmm  
                coverage = get_coverage(combined_domain["segments_hmm"])
                coverage = round(coverage /  combined_domain["hmm_length"],2)
                print("Coverage : "+str(coverage))
                # Calcul des positions
                positions = get_positions(combined_domain["segments_hmm"])
                if sequence in summarised_data['SeqID'].values:
                    summarised_data.loc[summarised_data['SeqID'] == sequence, f"Nb {domain} hits"] = len(sequence_domains[sequence])
                    summarised_data.loc[summarised_data['SeqID'] == sequence, f"Nb {domain} domains"] = len(merged_domains)
                    summarised_data.loc[summarised_data['SeqID'] == sequence, f"{domain} domain start"] = int(combined_domain["start_in_prot"])
                    summarised_data.loc[summarised_data['SeqID'] == sequence, f"{domain} domain end"] = int(combined_domain["end_in_prot"])
                    summarised_data.loc[summarised_data['SeqID'] == sequence, f"{domain} Score"] = combined_domain["global_score"]
                    summarised_data.loc[summarised_data['SeqID'] == sequence, f"{domain} position"] = str(positions)
                    summarised_data.loc[summarised_data['SeqID'] == sequence, f"{domain} Length"] = combined_domain["hmm_length"]
                    summarised_data.loc[summarised_data['SeqID'] == sequence, f"{domain} Prot Length"] = cumulated_length
                    summarised_data.loc[summarised_data['SeqID'] == sequence, f"{domain} coverage"] = coverage

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
                # Calcul longueur cumullee en proteine    
                cumulated_length = get_coverage(best_domain["segments_prot"])
                print("Cumulated length : "+str(cumulated_length))
                # Calcul couverture hmm  
                coverage = get_coverage(best_domain["segments_hmm"])
                coverage = round(coverage /  best_domain["hmm_length"],2)
                print("Coverage : "+str(coverage))                
                if sequence in summarised_data['SeqID'].values:
                    summarised_data.loc[summarised_data['SeqID'] == sequence, f"Nb {domain} hits"] = len(sequence_domains[sequence])
                    summarised_data.loc[summarised_data['SeqID'] == sequence, f"Nb {domain} domains"] = len(merged_domains)
                    summarised_data.loc[summarised_data['SeqID'] == sequence, f"{domain} domain start"] = best_domain["start_in_prot"]
                    summarised_data.loc[summarised_data['SeqID'] == sequence, f"{domain} domain end"] = best_domain["end_in_prot"]
                    summarised_data.loc[summarised_data['SeqID'] == sequence, f"{domain} Score"] = best_domain["hit_score"]
                    summarised_data.loc[summarised_data['SeqID'] == sequence, f"{domain} position"] = str(best_domain["segments_hmm"])
                    summarised_data.loc[summarised_data['SeqID'] == sequence, f"{domain} Length"] = best_domain["hmm_length"]
                    summarised_data.loc[summarised_data['SeqID'] == sequence, f"{domain} Prot Length"] = cumulated_length
                    summarised_data.loc[summarised_data['SeqID'] == sequence, f"{domain} coverage"] = coverage
                    

def process_domain_summary(domain, domain_summary_file, accession_number=accession_number):
    '''
    Récupère les informations importantes (nombre de domaines identifiés, position) dans les fichiers résultat de hmm_search --domtblout et les saisit dans un dataframe
    '''
    with open(domain_summary_file) as reader:    
        for line in reader.readlines()[1:]:
            line_data = line.split('\t')
            nb_domains = int(line_data[9])
            seq_id = line_data[0]
            if seq_id in summarised_data['SeqID'].values:
                summarised_data.loc[summarised_data['SeqID'] == seq_id, f"Nb {domain} domains"] = nb_domains
                summarised_data.loc[summarised_data['SeqID'] == seq_id, f"{domain} domain start"] = int(line_data[17])
                summarised_data.loc[summarised_data['SeqID'] == seq_id, f"{domain} domain end"]= int(line_data[18])

      
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
noms_colonnes.append('Nb '+domain+' domains')
noms_colonnes.append(domain+' Prot Length')
noms_colonnes.append(domain+' Length')
noms_colonnes.append(domain+' domain start')
noms_colonnes.append(domain+' domain end')
noms_colonnes.append(domain+' position')
noms_colonnes.append(domain+' coverage')
data_list = []


print(".... Accession  " + accession_number )
print(".... Domain  " + domain )
# All candidates must have the domain
print(".... Processing file "+domain_per_sequence_tabulated_file)

with open(domain_per_sequence_tabulated_file) as reader:
    for line in reader:
        line_data = line.strip().split('\t')
        # TO DO A QOI SERT SCORE ER EVALU ICI
        to_add = {'SeqID': line_data[0], 'Assembly':accession_number, domain+' Query': line_data[2], domain+' E-value': float(line_data[7]), domain+' Score': float(line_data[8])}
        data_list.append(to_add)
    summarised_data = pd.DataFrame(data_list, columns=noms_colonnes)
    summarised_data = summarised_data.astype({domain+' domain start': "Int32"})
    summarised_data = summarised_data.astype({domain+' domain end': "Int32"})
    summarised_data = summarised_data.astype({ 'Nb '+domain+' domains': "Int32"})    
    summarised_data = summarised_data.astype({'Nb '+ domain+' hits': "Int32"})  
    summarised_data = summarised_data.astype({domain+' position': "str"})    

print(".... Processing file "+domain_per_domain_summary_file)
process_domain_tabulated(domain,domain_per_domain_summary_file)
process_domain_merge(domain,domain_per_domain_summary_file)  

getTaxid()                
getSpecies()                

summarised_data = summarised_data.astype({domain+' Prot Length': "Int32"})   
summarised_data = summarised_data.astype({domain+' Length': "Int32"}) 
summarised_data['Nb '+domain+' domains'] = summarised_data['Nb '+domain+' domains'].fillna(value=0)
summarised_data['Nb '+domain+' hits'] = summarised_data['Nb '+domain+' hits'].fillna(value=0)
print("Output file = "+output_file)                 
#summarised_data.to_csv(output_file, sep=';',index=False,na_rep="N/A")
summarised_data.to_csv(output_file, sep=';',na_rep="NA") # On garde l'index car il est utilise par le script python suivant
