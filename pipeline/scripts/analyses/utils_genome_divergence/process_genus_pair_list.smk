import json
import glob
import os


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


# Name of global results directory.
# ---------------------------------
# The directory is located in pathGTDriftGlobalResults
GLOBAL_RESULTS = config["global_analyse_dir_name"]

# Name of genome specific results directory.
# ------------------------------------------
# The directory is located in genome_assembly/{accession}/analyses/
GENOME_RESULTS = config["analyse_dir_name"]

#BUSCO_DIR = config["busco_dir"]

# Function to load JSON files
def load_json(file_path):
    with open(file_path, 'r') as file:
        return json.load(file)

# Assign environment variables
globals().update(load_json("../environment_path.json"))

rule all:
    input:
        fam =  pathGTDriftGlobalResults + GLOBAL_RESULTS  + "list_pairs_genus",
        singles =  pathGTDriftGlobalResults + GLOBAL_RESULTS + "list_orphans_genus",
        busco =  pathGTDriftGlobalResults + GLOBAL_RESULTS + "busco_pairs_genus"

rule make_blast_db:
    input:
        busco =  pathGTDriftGlobalResults + GLOBAL_RESULTS + "busco_full.fa"
    output:
        pathGTDriftGlobalResults + GLOBAL_RESULTS + "busco_full.fa.ndb"

    shell:
        """
        makeblastdb -in {input.busco} -dbtype nucl -parse_seqids
        """

rule get_passed_busco:
    """
    Gets the list of species which passed BUSCO with at least one single or multi copy gene
    """
    input:
        busco =  pathGTDriftGlobalResults + GLOBAL_RESULTS + "busco_full.fa",
        db =  pathGTDriftGlobalResults + GLOBAL_RESULTS + "busco_full.fa.ndb"
    output:
         pathGTDriftGlobalResults + GLOBAL_RESULTS + "full_species_list"
    shell:
        """
        grep '>' {input.busco} |awk -F'\t' '{{ print $1 }}' |awk -F'-' '{{ print $2 }}' |uniq > {output}
        """

rule get_genus_list:
    input:
        #tax_data = pathGTDriftResource + "ncbi_dataset_eukaryota.taxonomy",
        tax_data = "gca_taxonomy.tsv",
        species =  pathGTDriftGlobalResults + GLOBAL_RESULTS + "full_species_list"
    output:
        fam =  pathGTDriftGlobalResults + GLOBAL_RESULTS + "list_pairs_genus",
        singles =  pathGTDriftGlobalResults + GLOBAL_RESULTS + "list_orphans_genus"
    shell:
        """
        python3 ../utils/python/generate_pairs.py -i {input.species} -t {input.tax_data} -l "Genus" -o {output.fam} -s {output.singles}
        """

rule create_pair_list_genus:
    input:
        pairs =  pathGTDriftGlobalResults + GLOBAL_RESULTS + "list_pairs_genus",
        busco =  pathGTDriftGlobalResults + GLOBAL_RESULTS + "busco_full.fa"
    output:
         pathGTDriftGlobalResults + GLOBAL_RESULTS + "busco_pairs_genus"
    shell:
        """
        python3  ../utils/python/create_busco_pairs.py -i {input.pairs} -b {input.busco} -o {output}
        """