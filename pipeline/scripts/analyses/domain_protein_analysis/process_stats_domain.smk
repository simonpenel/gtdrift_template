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
configfile: "config.json"

ACCESSNB = config["assembly_list"]

DOMAIN_HMM_DIR  = config["domain_hmm_dir"]

REFERENCE_DOMAIN_HMM = config["reference_domain_hmm"]

REFERENCE = config["reference"]

DOMAINS = config["domains"]

GLOBAL_RESULTS = config["global_results_dir"]

GENOME_RESULTS = config["genome_results_dir"]

if config["mode"] == "guix":
    RUNCMD = "guix shell hmmer -- "
else:
    RUNCMD = ""
    

# Rules
# -----

# --------------------------------------------------------------
# all : inputs define the files generated at the end of process. 
# --------------------------------------------------------------
rule all:
    """
    Get the domain stats on protein data
    """
    input:
        # Candidates after pararlog checks
        candidates=expand(
            pathGTDriftData
            + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS
            + "candidates.fasta",
            accession=ACCESSNB,
        ),
        
        

# Modules snakemake
# -----------------

# -----------------------------------------------------------
# module_stats_domain.smk:
# general statistics on domain search in protein data.
# outputs : stats_domain
# -----------------------------------------------------------
include: "../utils/module_stats_domain.smk"

include: "../utils/module_check_paralogs.smk"

