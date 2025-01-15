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

# List of assemblies
ACCESSNB = config["assembly_list"]

# Name of the directory containing fasta alignment files for each domain. 
# The directory is located in pathGTDriftResource/ref_align/
# -----------------------------------------------------------------------
DOMAIN_HMM_DIR  = config["domain_hmm_dir"]


# Path of the HMM profile used to check for paralogs. 
# The path is related  to pathGTDriftResource/
# ---------------------------------------------------
REFERENCE_DOMAIN_HMM = config["reference_domain_hmm"]

# The reference domain used to filter out paralogs
# ------------------------------------------------
REFERENCE = config["reference"]

# List of domains. Preliminary candidates are selected only if the first domain is present.
# -----------------------------------------------------------------------------------------
DOMAINS = config["domains"]

# Name of global results directory.
# The directory is located in pathGTDriftGlobalResults
# ----------------------------------------------------
GLOBAL_RESULTS = config["global_results_dir"]

# Name of genome specific results directory.
# The directory is located in genome_assembly/{accession}/analyses/
# -----------------------------------------------------------------
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
    Get the candidates in fasta format
    """
    input:
        # Candidates after pararlog checks
        candidates=expand(
            pathGTDriftData
           + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS
            + "candidates_{domain}.fasta",
            accession=ACCESSNB,domain=DOMAINS
        ),
        # Hmm profiles for paralogy checking
        hmm_db = expand(pathGTDriftResource + "hmm_build/paralogy_check/" + DOMAIN_HMM_DIR +"profil_{domain}.hmm",domain=DOMAINS),

               
        

# Modules snakemake
# -----------------

include: "../utils/module_stats_domain.smk"
include: "../utils/module_check_paralogs.smk"
include: "../utils/module_build_paralogs.smk"
