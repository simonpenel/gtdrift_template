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
DOMAIN = ["KRAB", "SET", "SSXRD", "ZF"]
GLOBAL_RESULTS = "prm9_protein_analysis/"

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
        stats_prdm9=expand(
            pathGTDriftData
            + "genome_assembly/{accession}/analyses/prdm9_prot/"
            + "summary_hmmsearch_prdm9_{accession}.csv",
            accession=ACCESSNB,
        ),
        
        # global statistics on KRAB domains
        # ---------------------------------
        global_krab=pathGTDriftGlobalResults + GLOBAL_RESULTS + "table_results/krab_data.csv",
        
        # global statistics on KRAB and ZF domains
        # ----------------------------------------
        global_krabzf=pathGTDriftGlobalResults
        + GLOBAL_RESULTS + "table_results/krabzf_data.csv",
        
        # global statistics on ZF domains count
        # -------------------------------------
        global_zf_count=pathGTDriftGlobalResults 
        + GLOBAL_RESULTS + "table_results/zf_count.csv",
        
        # global statistics on prdm9
        # --------------------------
        global_table_prdm9=pathGTDriftGlobalResults
        + GLOBAL_RESULTS + "table_results/table_prdm9.csv",
        
        # summary of prdm9 candidates
        # ---------------------------
        global_PRDM9_candidates=pathGTDriftGlobalResults
        + GLOBAL_RESULTS + "table_results/global_prdm9_candidates.csv",
        
        # statistics on zinc finger diversity
        # -----------------------------------
        global_zincfinger_diversity=pathGTDriftGlobalResults
        + GLOBAL_RESULTS + "table_results/zinc_finger.csv",
        
        # statistics on SET tyrosines
        # ---------------------------
        global_SET_tyrosines_stats=pathGTDriftGlobalResults
        + GLOBAL_RESULTS + "table_results/SET_tyrosines.csv",


# Modules snakemake
# -----------------

# -----------------------------------------------------------
# module_stats_prdm9.smk:
# general statistics on prdm9 search in protein data.
# outputs : stats_prdm9
# -----------------------------------------------------------
include: "module_stats_prdm9.smk"

# -----------------------------------------------------------
# module_stats_zinfinger.smk:
# statistics on zinc fingers in protein data, summary of data
# outputs : global_PRDM9_candidates, global_zincfinger_diversity, global_table_prdm9
# -----------------------------------------------------------
include: "module_stats_zincfinger.smk"

# -----------------------------------------------------------
# module_SET_tyrosines.smk:
# statistics on SET tyrosines
# outputs : global_SET_tyrosines_stats
# -----------------------------------------------------------
include: "module_SET_tyrosines.smk"

# -----------------------------------------------------------
# module_concatenate_results.smk:
# concatenate
# outputs : global_krab,global_krabzf,global_zf_count
# -----------------------------------------------------------
include: "module_concatenate_results.smk"
