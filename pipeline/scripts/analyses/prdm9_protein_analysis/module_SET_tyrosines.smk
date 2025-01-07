import json
import os
import glob
import pandas as pd

configfile: "config.json"
ACCESSNB = config["assembly_list"]

# Function to load JSON files
def load_json(file_path):
    with open(file_path, 'r') as file:
        return json.load(file)

# Assign environment variables
globals().update(load_json("../environment_path.json"))

# Rule to analyze extracted sequences using the SET_tyrosines.py script
rule analyze_prdm9_candidates:
    input:
        fasta_file = pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/candidates_prdm9.fasta"
    output:
        csv_output = pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/SET_tyrosines.csv"
    shell:
        """
        python3 python/SET_tyrosines.py {input.fasta_file} {output.csv_output}
        """

# Rule to combine all SET_tyrosines CSV files into one
rule combine_set_tyrosines:
    input:
        expand(pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/SET_tyrosines.csv", accession=ACCESSNB)
    output:
        combined_csv = pathGTDriftGlobalResults + GLOBAL_RESULTS + "table_results/SET_tyrosines.csv"
    run:
        # Combine all CSV files into one DataFrame
        dfs = [pd.read_csv(file) for file in input]
        combined_df = pd.concat(dfs, ignore_index=True)
        combined_df.to_csv(output.combined_csv, index=False)

