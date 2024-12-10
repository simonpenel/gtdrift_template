import sys

# This script is used to generate the json and query used in json files of other scripts
# python3 generate_json_and_query.py data/resources/organisms_list conf_assembly_list conf_query 

curated_only = False
if len(sys.argv) > 4:
	if sys.argv[4] == "curated":
	    curated_only = True
	    print("Curated only")
dico = {}
with open(sys.argv[1]) as reader:
    data = reader.readlines()
    for line in data[1:]:
        line_data = line.strip().split('\t')
        print(line_data[-1]+" "+ line_data[3])
        full_ac = line_data[2]
        url = line_data[-1]
        is_annoted = line_data[3]
        acc = full_ac.split(".")[0][4:]
        print ("url "+url+" cur onnly "+str(curated_only)+ " is-ano " + is_annoted)      
        if ((url != 'None') and (( curated_only == False ) or ( is_annoted == "True" ))) :
            print("OK")
            if acc in dico:
                dico[acc].append(full_ac)
            else :
                dico[acc] = [full_ac]
unic_dico = {}
for key in dico:
    val  = dico[key]
    if len(val) == 1 :
        unic_dico[key] = val[0]
    if len(val) == 2 :
        test = 0
        for candidat in val:
            if candidat[0:3] == "GCF" :
                unic_dico[key] = candidat
                test += 1
        if test != 1 :
            sys.exit("we want one and only one  GCF")
    if len(val) > 2 :
        print(val)
        sys.exit("too many assemblies")
accessions = list(unic_dico.values())
with  open(sys.argv[2]+".col", 'w') as writer_col,  open(sys.argv[3]+".col", 'w') as writer2_col:
    
    writer_col.write("  \"assembly_list\": [\n")
    writer_col.write('"'+accessions[0]+'" ')
    writer2_col.write('('+accessions[0].split(".")[0]+') ')
    for line in accessions[1:]:
            writer_col.write(', \n"'+line+'" ')
            writer2_col.write('\n OR ('+line.split(".")[0]+') ')
    writer_col.write("\n]")
