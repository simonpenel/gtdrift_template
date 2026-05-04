import sys

# This script is used to generate the json and query used in json files of other scripts
# python3 generate_json_and_query.py data/resources/organisms_list conf_assembly_list conf_query 

CURATED = []
UNCURATED = []
ALL = []
maxln = int(sys.argv[2])
nl = 0
index = 1
writer = open(sys.argv[1]+"."+str(maxln)+"."+str(index), 'w')
with open(sys.argv[1]) as reader:
    data = reader.readlines()
    header = data[0]
    writer.write(header)
    for line in data[1:]:
        nl = nl + 1
        if nl > maxln:
            writer.close()
            index = index + 1
            nl = 1
            writer = open(sys.argv[1]+"."+str(maxln)+"."+str(index), 'w')
            writer.write(header)
        writer.write(line)
writer.close()            

