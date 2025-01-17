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
# ------------------
ACCESSNB = config["assembly_list"]

# Name of the subdirectory used in several places for hmm calculations:
# - in pathGTDriftResource/ref_align/:
#   it contains the fasta alignment files for each domain. 
# - in pathGTDriftResource/hmm_build/search/
#   it contains the hmm profile for each domain.
# - in pathGTDriftResource/ref_align_for_paralogy_check/
#   it contains a directory for each domain, wich contains several 
#   alignnments of the domain, which will be used to build a hmm database
#   for each domain.
# - in pathGTDriftResource/hmm_build/paralogy_check/
#   it contains the hmm database for each domain for paralogy checks. 
# -----------------------------------------------------------------------
DOMAIN_HMM_DIR  = config["domain_hmm_dir"]

# The reference domain used to filter out paralogs during paralogy checks
# -----------------------------------------------------------------------
REFERENCE = config["reference"]

# List of domains to be processed. 
# ------------'-------------------
DOMAINS = config["domains"]

# List of domains to be processed. 
# ------------'-------------------
DOMAINS_SIMPLE = config["domains_simple"]

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
    
def get_reference(wildcards):
    fname = DOMAINS.get(wildcards.domain, "")
    return fname 

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
        # --------------------------------
        candidates=expand(
            pathGTDriftData
           + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS
            + "candidates_{domain}.fasta",
            accession=ACCESSNB,domain=DOMAINS
        ),
        # Hmm profiles for paralogy checking
        # ----------------------------------
        hmm_db = expand(pathGTDriftResource + "hmm_build/paralogy_check/" + DOMAIN_HMM_DIR +"profil_{domain}.hmm",domain=DOMAINS),
        # Candidates without pararlog checks
        # --------------------------------
        candidates_simple=expand(
            pathGTDriftData
           + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS
            + "candidates_1_{domain}.fasta",
            accession=ACCESSNB,domain=DOMAINS_SIMPLE
        ),
        summary_simple=expand(
            pathGTDriftData
            + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS +"summary_hmmsearch_{accession}_{domain}.csv",
            accession=ACCESSNB,domain=DOMAINS_SIMPLE)
               
        

# Modules snakemake
# -----------------

include: "../utils/module_stats_domain.smk"
include: "../utils/module_check_paralogs.smk"
include: "../utils/module_build_paralogs.smk"
