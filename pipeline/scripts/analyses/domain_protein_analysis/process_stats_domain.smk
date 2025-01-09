## Main snakemake for prdm9 analysis on protein data
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
configfile: "config.json"


ACCESSNB = config["assembly_list"]
DOMAIN = ["KRAB"]

GLOBAL_RESULTS = "domain_protein_analysis/"
GENOME_RESULTS = "domain_protein_analysis/"
# Rules
# -----

# --------------------------------------------------------------
# all : inputs define the files generated at the end of process. 
# --------------------------------------------------------------
rule all:
    """
    Get the prdm9 stats on protein data
    """
    input:
        # summary of prdm9 hmmsearch on the proteome for each assembly.
        # -------------------------------------------------------------
        stats_domain=expand(
            pathGTDriftData
            + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + 
            + "summary_hmmsearch_prdm9_{accession}.csv",
            accession=ACCESSNB,
        ),
        
        

# Modules snakemake
# -----------------

# -----------------------------------------------------------
# module_stats_domain.smk:
# general statistics on domain search in protein data.
# outputs : stats_domain
# -----------------------------------------------------------
include: "module_stats_domain.smk"



