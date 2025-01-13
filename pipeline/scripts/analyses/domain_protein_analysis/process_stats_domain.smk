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
        # summary of domain hmmsearch on the proteome for each assembly.
        # -------------------------------------------------------------
        stats_domain=expand(
            pathGTDriftData
            + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS
            + "summary_hmmsearch_{accession}.csv",
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



