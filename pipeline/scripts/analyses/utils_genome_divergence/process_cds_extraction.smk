import json

# Function to load JSON files
def load_json(file_path):
    with open(file_path, 'r') as file:
        return json.load(file)

# Assign environment variables
globals().update(load_json("../environment_path.json"))

configfile: "analyse.json"
configfile: "assemblies_annoted.json"
configfile: "assemblies_not_annoted.json"
ACCESSNB_ANNOTED = config["assembly_list_annoted"]
ACCESSNB_NOT_ANNOTED = config["assembly_list_not_annoted"]
nb_not_annoted = len(ACCESSNB_NOT_ANNOTED )
print(len(ACCESSNB_NOT_ANNOTED ))
annoted_set = set(ACCESSNB_ANNOTED)
ACCESSNB_NOT_ANNOTED = [
    x for x in ACCESSNB_NOT_ANNOTED
    if x not in annoted_set
]
print(len(ACCESSNB_NOT_ANNOTED ))
if len(ACCESSNB_NOT_ANNOTED ) < nb_not_annoted :
    print("Acessions have been removed from not annoated accessions")


#BUSCO_DIR = config["busco_dir"]

# Name of global results directory.
# ---------------------------------
# The directory is located in pathGTDriftGlobalResults
GLOBAL_RESULTS = config["global_analyse_dir_name"]

# Name of genome specific results directory.
# ------------------------------------------
# The directory is located in genome_assembly/{accession}/analyses/
GENOME_RESULTS = config["analyse_dir_name"]

rule all:
    """
    Get extracted sequences for BUSCOs
    """
    input:
        pair_list = pathGTDriftGlobalResults + GLOBAL_RESULTS + "busco_full.fa"


rule concatenate_all_buscos:
    """
    Concatenate BUSCO files of all species
    """
    input:
        busco_prot = expand(pathGTDriftGlobalResults + GLOBAL_RESULTS 
        + "extracted_buscos/{accession}_all_buscos.fa", accession=ACCESSNB_ANNOTED),
        busco_dna = expand(pathGTDriftGlobalResults + GLOBAL_RESULTS 
        + "extracted_buscos/{accession}_all_buscos.fa", accession=ACCESSNB_NOT_ANNOTED)
    output:
        busco_cat = pathGTDriftGlobalResults + GLOBAL_RESULTS + "busco_full.fa"
    shell:
        """
        cat {input.busco_prot} {input.busco_dna} > {output}
        """