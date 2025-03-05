import json
import os
import pandas as pd
import glob

# Configuration
# -------------
configfile: "analyse.json"
configfile: "assemblies.json"


# List of assemblies
# ------------------
ACCESSNB = config["assembly_list"]

# Name of global results directory.
# The directory is located in pathGTDriftGlobalResults
# ----------------------------------------------------
GLOBAL_RESULTS = config["analyse_dir_name"]

# Name of genome specific results directory.
# The directory is located in genome_assembly/{accession}/analyses/
# -----------------------------------------------------------------
GENOME_RESULTS = config["analyse_dir_name"]

# Function to load JSON files
def load_json(file_path):
    with open(file_path, 'r') as file:
        return json.load(file)

# Assign environment variables
globals().update(load_json("../environment_path.json"))

# Rule all
rule all:
    input:
        tyrosine_output = pathGTDriftGlobalResults 
        + GLOBAL_RESULTS + "all_SET_tyrosines.csv"


# Rule to analyze extracted sequences using the SET_tyrosines.py script
rule analyze_prdm9_candidates:
    input:
        fasta_file = pathGTDriftData
            + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS
#            + "selected_candidates.fasta",
            + "candidates_SET.fasta",
    output:
        csv_output = pathGTDriftData
            + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS
            + "SET_tyrosines.csv",
    shell:
        """
        python3 ../utils/python/SET_tyrosines.py {input.fasta_file} {output.csv_output}
        """

# Rule to combine all SET_tyrosines CSV files into one
rule combine_set_tyrosines:
    input:        
        expand(pathGTDriftData 
        + "genome_assembly/{accession}/analyses/" 
        + GENOME_RESULTS + "SET_tyrosines.csv", accession=ACCESSNB)
    output:
        pathGTDriftGlobalResults + GLOBAL_RESULTS + "all_SET_tyrosines.csv"
    script:
        "../utils/python/merge_tyrosines.py"             



