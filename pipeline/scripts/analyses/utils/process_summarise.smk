## Main snakemake for domain analysis on protein data
## Date : Decembre 2024
## Authors :


# Import python modules
# ---------------------
import os
import json

# Function to load JSON files
# ---------------------------
def load_json(file_path):
    with open(file_path, "r") as file:
        return json.load(file)

# Assign environment variables
# ----------------------------
globals().update(load_json("../environment_path.json"))


# Configuration
# -------------
configfile: "analyse.json"
configfile: "assemblies.json"

# List of assemblies
# ------------------
ANALYSES = config["analyse_list"]

# List of assemblies
# ------------------
ACCESSNB = config["assembly_list"]

GLOBAL_RESULTS = config["analyse_dir_name"]

GENOME_RESULTS = config["analyse_dir_name"]

rule all:
    input:
    # Joining  all the analyses for each accession
        resume=expand(
            pathGTDriftData
            + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS
            + "all_analyses_summary.txt",
            accession=ACCESSNB),
rule concat:
    """
    Join analyses
    """
    input:  
        # Analyse summaries for each analyse
        expand(
            pathGTDriftData
            + "genome_assembly/{{accession}}/analyses/{analyse}/"
            + "whole_summary.csv",
            analyse=ANALYSES),

    output:
            pathGTDriftData
            + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS
            + "all_analyses_summary.txt",
    script:
        "../utils/python/concatenate.py"



