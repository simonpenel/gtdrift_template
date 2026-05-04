import sys

# This script is used to generate the json and query used in json files of other scripts
# python3 generate_json_and_query.py data/resources/organisms_list conf_assembly_list conf_query 

CURATED = []
UNCURATED = []
ALL = []
with open(sys.argv[1]) as reader:
    data = reader.readlines()
    for line in data[1:]:
        #print(line)
        line_data = line.strip().split('\t')
        if line_data[-1] != 'None' and line_data[3] == 'True': # if there is an existing URL and genome is curated
                CURATED.append(line_data[2])
                ALL.append(line_data[2])                
        elif line_data[-1] != 'None':
                UNCURATED.append(line_data[2])
                ALL.append(line_data[2])

writer = open(sys.argv[2], 'w')
writer.write("{\n\
    \"storagetype\": \"irods\",\n\
    \"assembly_list\": [\n")
ac=ALL[0]
writer.write("\""+ac+"\"")   
for ac in ALL[1:]:
    writer.write(",\n\""+ac+"\"")
writer.write("\n]\n}\n")
writer.close()

writer = open("curated_only_"+sys.argv[2], 'w')
writer.write("{\n\
    \"storagetype\": \"irods\",\n\
    \"assembly_list\": [\n")
ac=CURATED[0]
writer.write("\""+ac+"\"")   
for ac in CURATED[1:]:
    writer.write(",\n\""+ac+"\"")
writer.write("\n]\n}\n")
